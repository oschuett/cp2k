!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module for equivariant PAO-ML based on PyTorch.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_model
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp,&
                                              sp
   USE message_passing,                 ONLY: mp_para_env_type
   USE pao_types,                       ONLY: pao_env_type,&
                                              pao_model_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE torch_api,                       ONLY: &
        torch_dict_create, torch_dict_get, torch_dict_insert, torch_dict_release, torch_dict_type, &
        torch_model_eval, torch_model_freeze, torch_model_get_attr, torch_model_load
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_model_load, pao_model_predict, pao_model_type

CONTAINS

! **************************************************************************************************
!> \brief Loads a PAO-ML model.
!> \param model ...
!> \param filename ...
!> \param atomic_kind_set ...
! **************************************************************************************************
   SUBROUTINE pao_model_load(model, filename, atomic_kind_set)
      TYPE(pao_model_type), INTENT(INOUT)                :: model
      CHARACTER(LEN=default_path_length), INTENT(IN)     :: filename
      TYPE(atomic_kind_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: atomic_kind_set

      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: feature_kind_names
      INTEGER                                            :: ikind, jkind

      WRITE (*, *) "Loading pao model from : ", TRIM(filename)
      CALL torch_model_load(model%torch_model, filename)

      ! Read model attributes.
      CALL torch_model_get_attr(model%torch_model, "pao_model_version", model%version)
      CALL torch_model_get_attr(model%torch_model, "kind_name", model%kind_name)
      CALL torch_model_get_attr(model%torch_model, "prim_basis_name", model%prim_basis_name)
      CALL torch_model_get_attr(model%torch_model, "pao_basis_size", model%pao_basis_size)
      CALL torch_model_get_attr(model%torch_model, "feature_kind_names", feature_kind_names)
      CALL torch_model_get_attr(model%torch_model, "num_neighbors", model%num_neighbors)
      CALL torch_model_get_attr(model%torch_model, "cutoff", model%cutoff)

      ! Freeze model after all attributes have been read.
      CALL torch_model_freeze(model%torch_model)

      ! For each feature kind name lookup its corresponding atomic kind number.
      ALLOCATE (model%feature_kinds(SIZE(feature_kind_names)))
      model%feature_kinds(:) = -1
      DO ikind = 1, SIZE(feature_kind_names)
         DO jkind = 1, SIZE(atomic_kind_set)
            IF (TRIM(atomic_kind_set(jkind)%name) == TRIM(feature_kind_names(ikind))) THEN
               model%feature_kinds(ikind) = jkind
            END IF
         END DO
         IF (model%feature_kinds(ikind) < 0) THEN
            WRITE (*, *) "Model supports feature kind '", TRIM(feature_kind_names(ikind)), "' that is not present in subsys."
         END IF
      END DO

      ! Check for missing kinds.
      DO jkind = 1, SIZE(atomic_kind_set)
         IF (ALL(model%feature_kinds /= atomic_kind_set(jkind)%kind_number)) THEN
            WRITE (*, *) "Model lacks feature kind '", TRIM(atomic_kind_set(jkind)%name), "' that is present in subsys."
         END IF
      END DO

      !TODO check compatibility
      WRITE (*, *) "pao_model_version: ", model%version

      !      IF (pao%iw > 0) &
      !         WRITE (pao%iw, "(A,I3,A,E10.1,1X,E10.1,1X,E10.1)") " PAO|ML| Descriptor for kind: "// &

   END SUBROUTINE pao_model_load

! **************************************************************************************************
!> \brief Fills pao%matrix_X based on machine learning predictions
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_model_predict(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_model_predict'

      INTEGER                                            :: acol, arow, handle, iatom, ikind, jatom, &
                                                            jkind, jneighbor, natoms
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: neighbors_indices
      REAL(dp), DIMENSION(3)                             :: relpos
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X
      REAL(sp), ALLOCATABLE, DIMENSION(:, :)             :: neighbors_features, neighbors_relpos
      REAL(sp), DIMENSION(:, :), POINTER                 :: predicted_xblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pao_model_type), POINTER                      :: model
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(torch_dict_type)                              :: model_inputs, model_outputs

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      cell=cell, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      natom=natoms)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,qs_env,particle_set,qs_kind_set,cell,natoms) &
!$OMP PRIVATE(iter,arow,acol,iatom,jatom,ikind,jkind,jneighbor,relpos,block_X) &
!$OMP PRIVATE(neighbors_indices, neighbors_relpos,neighbors_features) &
!$OMP PRIVATE(model, model_inputs,model_outputs, predicted_xblock)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         IF (SIZE(block_X) == 0) CYCLE ! pao disabled for iatom

         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         model => pao%models(ikind)
         CPASSERT(model%version > 0)

         ! Find neighbors.
         ! TODO support abitrary num_neighbors
         CPASSERT(model%num_neighbors == natoms - 1)
         ALLOCATE (neighbors_indices(model%num_neighbors))
         DO jneighbor = 1, model%num_neighbors
            IF (jneighbor < iatom) THEN
               neighbors_indices(jneighbor) = jneighbor
            ELSE
               neighbors_indices(jneighbor) = jneighbor + 1 ! skip over central iatom
            END IF
         END DO

         ! Compute neighbors relative positions.
         ALLOCATE (neighbors_relpos(3, model%num_neighbors))
         DO jneighbor = 1, model%num_neighbors
            jatom = neighbors_indices(jneighbor)
            ! TODO pbc
            relpos = particle_set(jatom)%r - particle_set(iatom)%r
            neighbors_relpos(:, jneighbor) = REAL(angstrom*relpos, kind=sp)
         END DO

         ! Compute neighbors features.
         ALLOCATE (neighbors_features(SIZE(model%feature_kinds), model%num_neighbors))
         neighbors_features(:, :) = 0.0_sp
         DO jneighbor = 1, model%num_neighbors
            jatom = neighbors_indices(jneighbor)
            jkind = particle_set(jatom)%atomic_kind%kind_number
            WHERE (model%feature_kinds == jkind) neighbors_features(:, jneighbor) = 1.0_sp
         END DO

         ! Inference.
         model_inputs = torch_dict_type()
         CALL torch_dict_create(model_inputs)
         CALL torch_dict_insert(model_inputs, "neighbors_relpos", neighbors_relpos)
         CALL torch_dict_insert(model_inputs, "neighbors_features", neighbors_features)
         model_outputs = torch_dict_type()
         CALL torch_dict_create(model_outputs)
         CALL torch_model_eval(model%torch_model, model_inputs, model_outputs)

         ! Copy predicted XBlock.
         NULLIFY (predicted_xblock)
         CALL torch_dict_get(model_outputs, "xblock", predicted_xblock)
         block_X = RESHAPE(predicted_xblock, (/SIZE(block_X), 1/))
         WRITE (*, *) "Predicted xblock for atom: ", iatom

         ! Clean up.
         CALL torch_dict_release(model_inputs)
         CALL torch_dict_release(model_outputs)
         DEALLOCATE (predicted_xblock, neighbors_indices, neighbors_relpos, neighbors_features)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_model_predict

END MODULE pao_model

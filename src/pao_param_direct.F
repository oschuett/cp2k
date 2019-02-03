!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Direct parametrization
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_direct
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_get_block_p, dbcsr_p_type, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_release, dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE pao_param_methods,               ONLY: pao_calc_grad_lnv_wrt_AB
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE pao_potentials,                  ONLY: pao_guess_initial_potential
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_direct'

   PUBLIC :: pao_param_init_direct, pao_param_finalize_direct, pao_calc_AB_direct
   PUBLIC :: pao_param_count_direct, pao_param_initguess_direct

CONTAINS

! **************************************************************************************************
!> \brief Initialize direct parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_direct(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      IF (pao%precondition) &
         CPABORT("PAO preconditioning not supported for selected parametrization.")

   END SUBROUTINE pao_param_init_direct

! **************************************************************************************************
!> \brief Finalize direct parametrization
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_direct()

      ! Nothing todo.

   END SUBROUTINE pao_param_finalize_direct

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_direct(qs_env, ikind, nparams)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_count_direct'

      INTEGER                                            :: pao_basis_size, pri_basis_size
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), &
                       basis_set=basis_set, &
                       pao_basis_size=pao_basis_size)
      pri_basis_size = basis_set%nsgf

      nparams = pao_basis_size*pri_basis_size

   END SUBROUTINE pao_param_count_direct

! **************************************************************************************************
!> \brief Fills matrix_X with an initial guess
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_direct(pao, qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initguess_direct'

      INTEGER                                            :: acol, arow, handle, iatom, m, n
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X, H, A, V0, block_H0, H_evecs, &
                                                            block_N, block_N_inv
      REAL(dp), DIMENSION(:), POINTER                    :: H_evals
      TYPE(dbcsr_iterator_type)                          :: iter
      LOGICAL                                            :: found

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,blk_sizes_pri,blk_sizes_pao) &
!$OMP PRIVATE(iter,arow,acol,iatom,n,m,i,block_X)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)

         CALL dbcsr_get_block_p(matrix=pao%matrix_H0, row=iatom, col=iatom, block=block_H0, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_diag, row=iatom, col=iatom, block=block_N, found=found)
         ! TODO: matrix_N_inv has wrong distribution,  pao%matrix_N_inv_diag does not yet exist
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=iatom, col=iatom, block=block_N_inv, found=found)
         !CALL dbcsr_get_block_p(matrix=pao%matrix_N_inv_diag, row=iatom, col=iatom, block=block_N_inv, found=found)
         CPASSERT(ASSOCIATED(block_H0) .AND. ASSOCIATED(block_N) .AND. ASSOCIATED(block_N_inv))

         n = blk_sizes_pri(iatom) ! size of primary basis
         m = blk_sizes_pao(iatom) ! size of pao basis

         ALLOCATE (V0(n, n))
         CALL pao_guess_initial_potential(qs_env, iatom, V0)

         ! construct H
         ALLOCATE (H(n, n))
         H = MATMUL(MATMUL(block_N, block_H0 + V0), block_N) ! transform into orthonormal basis
         !H = block_H0 + V0

         ! diagonalize H
         ALLOCATE (H_evecs(n, n), H_evals(n))
         H_evecs = H
         CALL diamat_all(H_evecs, H_evals)

         ! use first m eigenvectors as initial guess
         ALLOCATE (A(n, m))
         A = MATMUL(block_N_inv, H_evecs(:, 1:m))
         !A = H_evecs(:, 1:m)

         WRITE (*,*) "OLE atom:", iatom, " init guess ", A(:, 1)

         block_X = RESHAPE(A, (/n*m,1/))
         DEALLOCATE (H, V0, A, H_evecs, H_evals)

      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_param_initguess_direct

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_direct(pao, qs_env, ls_scf_env, gradient)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_AB_direct'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, m, n
      LOGICAL                                            :: found
      REAL(dp)                                           :: denom, nom, w
      REAL(dp), DIMENSION(:), POINTER                    :: AA_evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: AA, AA_evecs, AA_inv, block_A, block_B, &
                                                            block_G, block_Ma, block_Mb, block_X, &
                                                            D1, G, M1, M2, M3, block_N
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_Ma, matrix_Mb
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct

      CALL timeset(routineN, handle)
      ls_mstruct => ls_scf_env%ls_mstruct

      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)
      ENDIF

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,matrix_Ma,matrix_Mb) &
!$OMP PRIVATE(iter,arow,acol,iatom,found,n,m,w,i,j,k) &
!$OMP PRIVATE(block_X,block_A,block_B,AA,AA_evals,AA_evecs,AA_inv) &
!$OMP PRIVATE(block_Ma, block_Mb, block_G)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_A, row=iatom, col=iatom, block=block_A, found=found)
         CPASSERT(ASSOCIATED(block_A))
         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_B, row=iatom, col=iatom, block=block_B, found=found)
         CPASSERT(ASSOCIATED(block_B))
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_diag, row=iatom, col=iatom, block=block_N, found=found)
         CPASSERT(ASSOCIATED(block_N))


         n = SIZE(block_A, 1) ! size of primary basis
         m = SIZE(block_A, 2) ! size of pao basis
         block_A = RESHAPE(block_X, (/n, m/))

         !TODO: rename AA -> AS_invA
         ALLOCATE (AA(m, m))
         !AA = MATMUL(TRANSPOSE(block_A), block_A)
         !AA = MATMUL(MATMUL(TRANSPOSE(block_A), block_s), block_A)
         AA = MATMUL(MATMUL(MATMUL(TRANSPOSE(block_A), block_N), block_N), block_A)

         ! diagonalize AA
         ALLOCATE (AA_evecs(m, m), AA_evals(m))
         AA_evecs(:, :) = AA
         CALL diamat_all(AA_evecs, AA_evals)

         ! build AA_inv
         ALLOCATE (AA_inv(m, m))
         AA_inv(:, :) = 0.0_dp
         DO k = 1, m
            w = 1.0_dp / AA_evals(k)
            DO i = 1, m
            DO j = 1, m
               AA_inv(i, j) = AA_inv(i, j)+w*AA_evecs(i, k)*AA_evecs(j, k)
            ENDDO
            ENDDO
         ENDDO

         !B = 1/S * A * 1/(A^T 1/S A)
         block_B = MATMUL(MATMUL(MATMUL(block_N, block_N), block_A), AA_inv)

         WRITE (*,*) "OLE atom:", iatom, " block_A(:,1)", block_A(:,1)
         WRITE (*,*) "OLE atom:", iatom, " block_B(:,1)", block_B(:,1)

         WRITE (*,*) "OLE atom:", iatom, " A*B^T", MATMUL(block_A, TRANSPOSE(block_B))
         WRITE (*,*) "OLE atom:", iatom, " A^T*B", MATMUL(TRANSPOSE(block_A), block_B)

         ! TURNING POINT (if calc grad) ------------------------------------------
         IF (gradient) THEN
            CALL dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
            CPASSERT(ASSOCIATED(block_G))
            CALL dbcsr_get_block_p(matrix=matrix_Ma, row=iatom, col=iatom, block=block_Ma, found=found)
            CALL dbcsr_get_block_p(matrix=matrix_Mb, row=iatom, col=iatom, block=block_Mb, found=found)
            ! don't check ASSOCIATED(block_M), it might have been filtered out.

            ALLOCATE (G(n, m))
            G(:, :) = 0.0_dp
            IF (ASSOCIATED(block_Ma)) THEN
               G = block_Ma
            ENDIF

            IF (ASSOCIATED(block_Mb)) THEN

               G = G+MATMUL(block_Mb, AA_inv)

               ALLOCATE (D1(m, m), M1(m, m), M2(m, m), M3(m, m))
               ! apparently this term is really negelitable :-/

               ! calculate derivatives dAA_inv/ dAA
               DO i = 1, m
               DO j = 1, m
                  denom = AA_evals(i)-AA_evals(j)
                  IF (i == j) THEN
                     D1(i, i) = -1.0_dp/AA_evals(i)**2 ! diagonal elements
                  ELSE IF (ABS(denom) > 1e-10_dp) THEN
                     nom = -1.0_dp/AA_evals(i)**2+1.0_dp/AA_evals(j)**2
                     D1(i, j) = nom/denom
                  ELSE
                     D1(i, j) = 1.0_dp !TODO: ! limit according to L'Hospital's rule
                  ENDIF
               ENDDO
               ENDDO

               M1 = MATMUL(TRANSPOSE(block_A), block_Mb) !TODO: or the transposed?
               M2 = M1*D1 ! Hadamard product
               M3 = MATMUL(MATMUL(AA_evecs, M2), TRANSPOSE(AA_evecs))

               ! calculate derivatives dAA / dA
               G = G+2.0_dp*MATMUL(block_A, M3)
               DEALLOCATE (D1, M1, M2, M3)
            ENDIF

            block_G = RESHAPE(G, (/n*m, 1/))
            DEALLOCATE (G)
         ENDIF

         DEALLOCATE (AA, AA_evecs, AA_evals, AA_inv)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      IF (gradient) THEN
         CALL dbcsr_release(matrix_Ma)
         CALL dbcsr_release(matrix_Mb)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pao_calc_AB_direct

END MODULE pao_param_direct

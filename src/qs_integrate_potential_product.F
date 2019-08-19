!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Build up the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      Joost VandeVondele (02.2002)
!>            1) rewrote collocate_pgf for increased accuracy and speed
!>            2) collocate_core hack for PGI compiler
!>            3) added multiple grid feature
!>            4) new way to go over the grid
!>      Joost VandeVondele (05.2002)
!>            1) prelim. introduction of the real space grid type
!>      JGH [30.08.02] multigrid arrays independent from potential
!>      JGH [17.07.03] distributed real space code
!>      JGH [23.11.03] refactoring and new loop ordering
!>      JGH [04.12.03] OpneMP parallelization of main loops
!>      Joost VandeVondele (12.2003)
!>           1) modified to compute tau
!>      Joost removed incremental build feature
!>      Joost introduced map consistent
!>      Rewrote grid integration/collocation routines, [Joost VandeVondele,03.2007]
!>      JGH [26.06.15] modification to allow for k-points
!> \author Matthias Krack (03.04.2001)
! **************************************************************************************************
MODULE qs_integrate_potential_product
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE gaussian_gridlevels,             ONLY: gridlevel_info_type
   USE grid_ibase_mp,                   ONLY: grid_integrate_mp
   USE grid_ibase_ref,                  ONLY: grid_integrate_ref
   USE input_constants,                 ONLY: do_admm_exch_scaling_merlot
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_p_type,&
                                              rs_grid_release,&
                                              rs_grid_retain
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE task_list_methods,               ONLY: rs_distribute_matrix
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   INTEGER :: debug_count = 0

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_product'

! *** Public subroutines ***
! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

   PUBLIC :: integrate_v_rspace

CONTAINS

! **************************************************************************************************
!> \brief computes matrix elements corresponding to a given potential
!> \param v_rspace ...
!> \param hmat ...
!> \param hmat_kp ...
!> \param pmat ...
!> \param pmat_kp ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param force_adm whether force of in aux. dens. matrix is calculated
!> \param ispin ...
!> \param compute_tau ...
!> \param gapw ...
!> \param basis_type ...
!> \param pw_env_external ...
!> \param task_list_external ...
!> \par History
!>      IAB (29-Apr-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!>      Some refactoring, get priorities for options correct (JGH, 04.2014)
!>      Added options to allow for k-points
!>      For a smooth transition we allow for old and new (vector) matrices (hmat, pmat) (JGH, 06.2015)
!> \note
!>     integrates a given potential (or other object on a real
!>     space grid) = v_rspace using a multi grid technique (mgrid_*)
!>     over the basis set producing a number for every element of h
!>     (should have the same sparsity structure of S)
!>     additional screening is available using the magnitude of the
!>     elements in p (? I'm not sure this is a very good idea)
!>     this argument is optional
!>     derivatives of these matrix elements with respect to the ionic
!>     coordinates can be computed as well
! **************************************************************************************************
   SUBROUTINE integrate_v_rspace(v_rspace, hmat, hmat_kp, pmat, pmat_kp, &
                                 qs_env, calculate_forces, force_adm, ispin, &
                                 compute_tau, gapw, basis_type, pw_env_external, task_list_external)

      TYPE(pw_p_type)                                    :: v_rspace
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hmat
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: hmat_kp
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: pmat
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: pmat_kp
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL, INTENT(IN), OPTIONAL                      :: force_adm
      INTEGER, INTENT(IN), OPTIONAL                      :: ispin
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau, gapw
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: basis_type
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: my_basis_type
      INTEGER                                            :: code_base, handle, i, iatom, &
                                                            igrid_level, img, natom, nimages, &
                                                            offs_dv
      INTEGER(KIND=int_8), DIMENSION(:), POINTER         :: atom_pair_recv, atom_pair_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      LOGICAL :: distributed_grids, do_kp, h_duplicated, map_consistent, my_compute_tau, &
         my_force_adm, my_gapw, p_duplicated, pab_required, scatter, use_virial
      REAL(KIND=dp)                                      :: admm_scal_fac, eps_gvg_rspace
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: posat
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap, dhmat, htemp
      TYPE(dbcsr_type), POINTER                          :: href
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! we test here if the provided operator matrices are consistent
      CPASSERT(PRESENT(hmat) .OR. PRESENT(hmat_kp))
      do_kp = .FALSE.
      IF (PRESENT(hmat_kp)) do_kp = .TRUE.
      IF (PRESENT(pmat)) THEN
         CPASSERT(PRESENT(hmat))
      ELSE IF (PRESENT(pmat_kp)) THEN
         CPASSERT(PRESENT(hmat_kp))
      END IF

      NULLIFY (pw_env, rs_descs, admm_env)

      offs_dv = 0

      ! this routine works in two modes:
      ! normal mode : <a| V | b>
      ! tau mode    : < nabla a| V | nabla b>
      my_compute_tau = .FALSE.
      IF (PRESENT(compute_tau)) my_compute_tau = compute_tau

      my_force_adm = .FALSE.
      IF (PRESENT(force_adm)) my_force_adm = force_adm

      ! this sets the basis set to be used. GAPW(==soft basis) overwrites basis_type
      ! default is "ORB"
      my_gapw = .FALSE.
      IF (PRESENT(gapw)) my_gapw = gapw
      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      ! get the task lists
      ! task lists have to be in sync with basis sets
      ! there is an option to provide the task list from outside (not through qs_env)
      ! outside option has highest priority
      SELECT CASE (my_basis_type)
      CASE ("ORB")
         CALL get_qs_env(qs_env=qs_env, &
                         task_list=task_list, &
                         task_list_soft=task_list_soft)
      CASE ("AUX_FIT")
         CALL get_qs_env(qs_env=qs_env, &
                         task_list_aux_fit=task_list, &
                         task_list_soft=task_list_soft)
      END SELECT
      IF (my_gapw) task_list => task_list_soft
      IF (PRESENT(task_list_external)) task_list => task_list_external
      CPASSERT(ASSOCIATED(task_list))

      ! the information on the grids is provided through pw_env
      ! pw_env has to be the parent env for the potential grid (input)
      ! there is an option to provide an external grid
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      IF (PRESENT(pw_env_external)) pw_env => pw_env_external

      ! get all the general information on the system we are working on
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      force=force, &
                      virial=virial)

      admm_scal_fac = 1.0_dp
      IF (my_force_adm) THEN
         CALL get_qs_env(qs_env=qs_env, admm_env=admm_env)
         ! Calculate bare scaling of force according to Merlot, 1. IF: ADMMP, 2. IF: ADMMS,
         IF ((.NOT. admm_env%charge_constrain) .AND. &
             (admm_env%scaling_model == do_admm_exch_scaling_merlot)) THEN
            admm_scal_fac = admm_env%gsi(ispin)**2
         ELSE IF (admm_env%charge_constrain .AND. &
                  (admm_env%scaling_model == do_admm_exch_scaling_merlot)) THEN
            admm_scal_fac = (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)
         END IF
      END IF

      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
      END DO

      ! assign from pw_env
      gridlevel_info => pw_env%gridlevel_info
      cube_info => pw_env%cube_info

      ! transform the potential on the rs_multigrids
      CALL potential_pw2rs(rs_v, v_rspace, pw_env)

      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         CPASSERT(do_kp)
      END IF
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      ! iatom - ikind mapping
      natom = SIZE(particle_set)
      ALLOCATE (kind_of(natom))
      DO iatom = 1, natom
         kind_of(iatom) = particle_set(iatom)%atomic_kind%kind_number
      END DO
      ! atomic positions
      ALLOCATE (posat(3, natom))
      DO iatom = 1, natom
         posat(:, iatom) = pbc(particle_set(iatom)%r, cell)
      END DO
      IF (calculate_forces) THEN
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)
      END IF

      map_consistent = dft_control%qs_control%map_consistent
      IF (map_consistent) THEN
         ! needs to be consistent with rho_rspace
         eps_gvg_rspace = dft_control%qs_control%eps_rho_rspace
      ELSE
         eps_gvg_rspace = dft_control%qs_control%eps_gvg_rspace
      ENDIF

      pab_required = (PRESENT(pmat) .OR. PRESENT(pmat_kp)) &
                     .AND. (calculate_forces .OR. .NOT. map_consistent)

      distributed_grids = .FALSE.
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         IF (rs_v(igrid_level)%rs_grid%desc%distributed) THEN
            distributed_grids = .TRUE.
         ENDIF
      ENDDO

      ! initialize the working hmat structures
      h_duplicated = .FALSE.
      ALLOCATE (dhmat(nimages))
      IF (do_kp) THEN
         DO img = 1, nimages
            dhmat(img)%matrix => hmat_kp(img)%matrix
         END DO
      ELSE
         dhmat(1)%matrix => hmat%matrix
      END IF
      IF (distributed_grids) THEN
         h_duplicated = .TRUE.
         href => dhmat(1)%matrix
         DO img = 1, nimages
            NULLIFY (dhmat(img)%matrix)
            ALLOCATE (dhmat(img)%matrix)
            CALL dbcsr_create(dhmat(img)%matrix, template=href, name='LocalH')
         END DO
      END IF

      p_duplicated = .FALSE.
      IF (pab_required) THEN
         ! initialize the working pmat structures
         ALLOCATE (deltap(nimages))
         IF (do_kp) THEN
            DO img = 1, nimages
               deltap(img)%matrix => pmat_kp(img)%matrix
            END DO
         ELSE
            deltap(1)%matrix => pmat%matrix
         END IF
         IF (distributed_grids) THEN
            p_duplicated = .TRUE.
            DO img = 1, nimages
               NULLIFY (deltap(img)%matrix)
               ALLOCATE (deltap(img)%matrix)
            END DO
            IF (do_kp) THEN
               DO img = 1, nimages
                  CALL dbcsr_copy(deltap(img)%matrix, pmat_kp(img)%matrix, name="LocalP")
               END DO
            ELSE
               CALL dbcsr_copy(deltap(1)%matrix, pmat%matrix, name="LocalP")
            END IF
         END IF
      END IF

      atom_pair_send => task_list%atom_pair_send
      atom_pair_recv => task_list%atom_pair_recv
      IF (distributed_grids .AND. pab_required) THEN
         CALL rs_distribute_matrix(rs_descs, deltap, atom_pair_send, atom_pair_recv, &
                                   natom, nimages, scatter=.TRUE.)
      ENDIF

      code_base = 0
      SELECT CASE (code_base)
      CASE (0)
         CALL grid_integrate_ref(rs_v, dhmat, deltap, task_list, cube_info, gridlevel_info, cell, &
                                 eps_gvg_rspace, admm_scal_fac, &
                                 my_compute_tau, pab_required, map_consistent, atom_of_kind, kind_of, posat, &
                                 calculate_forces, force, use_virial, virial)
      CASE (1)
         CALL grid_integrate_mp(rs_v, dhmat, deltap, task_list, cube_info, gridlevel_info, cell, &
                                eps_gvg_rspace, admm_scal_fac, &
                                my_compute_tau, pab_required, map_consistent, atom_of_kind, kind_of, posat, &
                                calculate_forces, force, use_virial, virial)
      CASE DEFAULT
         CPABORT("NA")
      END SELECT

      IF (h_duplicated) THEN
         ! Reconstruct H matrix if using distributed RS grids
         ! note send and recv direction reversed WRT collocate
         scatter = .FALSE.
         IF (do_kp) THEN
            CALL rs_distribute_matrix(rs_descs, dhmat, atom_pair_recv, atom_pair_send, &
                                      natom, nimages, scatter, hmats=hmat_kp)
         ELSE
            ALLOCATE (htemp(1))
            htemp(1)%matrix => hmat%matrix

            CALL rs_distribute_matrix(rs_descs, dhmat, atom_pair_recv, atom_pair_send, &
                                      natom, nimages, scatter, hmats=htemp)

            NULLIFY (htemp(1)%matrix)
            DEALLOCATE (htemp)
         END IF
         CALL dbcsr_deallocate_matrix_set(dhmat)
      ELSE
         DO img = 1, nimages
            NULLIFY (dhmat(img)%matrix)
         END DO
         DEALLOCATE (dhmat)
      END IF

      IF (pab_required) THEN
         IF (p_duplicated) THEN
            CALL dbcsr_deallocate_matrix_set(deltap)
         ELSE
            DO img = 1, nimages
               NULLIFY (deltap(img)%matrix)
            END DO
            DEALLOCATE (deltap)
         END IF
      END IF

      ! Release work storage
      IF (ASSOCIATED(rs_v)) THEN
         DO i = 1, SIZE(rs_v)
            CALL rs_grid_release(rs_v(i)%rs_grid)
         END DO
      END IF

      DEALLOCATE (kind_of, posat)
      IF (calculate_forces) THEN
         DEALLOCATE (atom_of_kind)
      END IF

      CALL timestop(handle)

   END SUBROUTINE integrate_v_rspace

END MODULE qs_integrate_potential_product

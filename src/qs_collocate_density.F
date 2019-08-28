!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      - rewrote collocate for increased accuracy and speed
!>      - introduced the PGI hack for increased speed with that compiler
!>        (22.02.02)
!>      - Added Multiple Grid feature
!>      - new way to get over the grid (01.03.02)
!>      - removed timing calls since they were getting expensive
!>      - Updated with the new QS data structures (09.04.02,MK)
!>      - introduction of the real space grid type ( prelim. version JVdV 05.02)
!>      - parallel FFT (JGH 22.05.02)
!>      - multigrid arrays independent from density (JGH 30.08.02)
!>      - old density stored in g space (JGH 30.08.02)
!>      - distributed real space code (JGH 17.07.03)
!>      - refactoring and new loop ordering (JGH 23.11.03)
!>      - OpenMP parallelization (JGH 03.12.03)
!>      - Modified to compute tau (Joost 12.03)
!>      - removed the incremental density rebuild (Joost 01.04)
!>      - introduced realspace multigridding (Joost 02.04)
!>      - introduced map_consistent (Joost 02.04)
!>      - Addition of the subroutine calculate_atomic_charge_density (TdK, 08.05)
!>      - rewrite of the collocate/integrate kernels (Joost VandeVondele, 03.07)
!> \author Matthias Krack (03.04.2001)
!>      1) Joost VandeVondele (01.2002)
!>      Thomas D. Kuehne (04.08.2005)
! **************************************************************************************************
MODULE qs_collocate_density
   USE ao_util,                         ONLY: exp_radius_single
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE cube_utils,                      ONLY: cube_info_type
   USE d3_poly,                         ONLY: poly_cp2k2d3
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type
   USE gauss_colloc,                    ONLY: collocGauss
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                              gridlevel_info_type
   USE grid_base_mp,                    ONLY: collocate_pgf_product_rspace,&
                                              grid_collocate_mp
   USE grid_base_ref,                   ONLY: grid_base_type,&
                                              grid_collocate_ref
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE lri_environment_types,           ONLY: lri_kind_type
   USE mathconstants,                   ONLY: pi,&
                                              twopi
   USE memory_utilities,                ONLY: reallocate
   USE orbital_pointers,                ONLY: coset,&
                                              indco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_integrate_function,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type,&
                                              pw_pools_create_pws,&
                                              pw_pools_give_back_pws
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_modify_pab_block,             ONLY: &
        FUNC_AB, FUNC_ADBmDAB, FUNC_ARB, FUNC_ARDBmDARB, FUNC_DABpADB, FUNC_DADB, FUNC_DER, &
        FUNC_DX, FUNC_DXDX, FUNC_DXDY, FUNC_DY, FUNC_DYDY, FUNC_DYDZ, FUNC_DZ, FUNC_DZDX, FUNC_DZDZ
   USE realspace_grid_types,            ONLY: realspace_grid_desc_p_type,&
                                              realspace_grid_p_type,&
                                              realspace_grid_type,&
                                              rs2pw,&
                                              rs_grid_release,&
                                              rs_grid_retain,&
                                              rs_grid_zero,&
                                              rs_pw_transfer
   USE rs_pw_interface,                 ONLY: density_rs2pw
   USE task_list_methods,               ONLY: rs_distribute_matrix
   USE task_list_types,                 ONLY: task_list_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_collocate_density'
! *** Public subroutines ***

   PUBLIC :: calculate_ppl_grid, &
             calculate_rho_core, &
             calculate_rho_elec, &
             calculate_lri_rho_elec, &
             calculate_rho_single_gaussian, &
             calculate_rho_metal, &
             calculate_rho_resp_single, &
             calculate_rho_resp_all, &
             calculate_wavefunction, &
             collocate_pgf_rspace, &
             collocate_pgf_product_gspace, &
             collocate_pgf_product_rspace, &
             calculate_rho_nlcc

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief computes the density of the non-linear core correction on the grid
!> \param rho_nlcc ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_rho_nlcc(rho_nlcc, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_nlcc
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_nlcc', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, iexp_nlcc, ikind, &
                                                            ithread, j, lmax_global, n, natom, nc, &
                                                            nexp_nlcc, ni, npme, nthread
      INTEGER(KIND=int_8)                                :: subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, mylmax, nct_nlcc
      LOGICAL                                            :: nlcc
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_nlcc
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval_nlcc, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, pab, particle_set, atomic_kind_set, qs_kind_set, &
               atom_list, pw_env, rs_rho, auxbas_pw_pool, cores, mylmax)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      cube_info = pw_env%cube_info(1)
      ! be careful in parallel nsmax is choosen with multigrid in mind!
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      !Find the maximum la_max over all of the atomic_kind_set
      !We are using FUNC_AB so no need to cater for different ga_gb_functions
      lmax_global = 0
      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)
         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, nct_nlcc=mylmax, nlcc_present=nlcc)
         IF (.NOT. nlcc) CYCLE
         lmax_global = MAX(MAXVAL(mylmax), lmax_global)
      END DO
      !lmax_global set according to ni = 2*nc-2 below
      lmax_global = 2*lmax_global-2

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)

         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, nlcc_present=nlcc, nexp_nlcc=nexp_nlcc, &
                            alpha_nlcc=alpha_nlcc, nct_nlcc=nct_nlcc, cval_nlcc=cval_nlcc)

         IF (.NOT. nlcc) CYCLE

         DO iexp_nlcc = 1, nexp_nlcc

            alpha = alpha_nlcc(iexp_nlcc)
            nc = nct_nlcc(iexp_nlcc)

            ni = ncoset(2*nc-2)
            ALLOCATE (pab(ni, 1))
            pab = 0._dp

            nthread = 1
            ithread = 0

            CALL reallocate(cores, 1, natom)
            npme = 0
            cores = 0

            ! prepare core function
            DO j = 1, nc
               SELECT CASE (j)
               CASE (1)
                  pab(1, 1) = cval_nlcc(1, iexp_nlcc)
               CASE (2)
                  n = coset(2, 0, 0)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
                  n = coset(0, 2, 0)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
                  n = coset(0, 0, 2)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
               CASE (3)
                  n = coset(4, 0, 0)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 4, 0)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 0, 4)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(2, 2, 0)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(2, 0, 2)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 2, 2)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
               CASE (4)
                  n = coset(6, 0, 0)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 6, 0)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 0, 6)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(4, 2, 0)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(4, 0, 2)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 4, 0)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 0, 4)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 4, 2)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 2, 4)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 2, 2)
                  pab(n, 1) = 6._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
               CASE DEFAULT
                  CPABORT("")
               END SELECT
            END DO
            IF (dft_control%nspins == 2) pab = pab*0.5_dp

            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
                  ! replicated realspace grid, split the atoms up between procs
                  IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
                     npme = npme+1
                     cores(npme) = iatom
                  ENDIF
               ELSE
                  npme = npme+1
                  cores(npme) = iatom
               ENDIF
            END DO

            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               subpatch_pattern = 0
               ni = 2*nc-2
               radius = exp_radius_single(0, ni, pab(:, 1), 1.0_dp/(2.0_dp*alpha**2), eps_rho_rspace, 1.0_dp)
               CALL collocate_pgf_rspace(ni, 1.0_dp/(2.0_dp*alpha**2), 0, ra, pab(:, 1), rs_rho, cell, radius, 1.0_dp)
            END DO

            DEALLOCATE (pab)

         END DO

      END DO

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF

      CALL rs_pw_transfer(rs_rho, rho_nlcc%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_nlcc

! **************************************************************************************************
!> \brief computes the local pseudopotential (without erf term) on the grid
!> \param vppl ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calculate_ppl_grid(vppl, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: vppl
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_ppl_grid', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, ikind, ithread, &
                                                            j, lmax_global, lppl, n, natom, ni, &
                                                            npme, nthread
      INTEGER(KIND=int_8)                                :: subpatch_pattern
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cexp_ppl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, pab, atomic_kind_set, qs_kind_set, particle_set, &
               atom_list, pw_env, rs_rho, auxbas_pw_pool, cores)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      cube_info = pw_env%cube_info(1)
      ! be careful in parallel nsmax is choosen with multigrid in mind!
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      !Find maximum la_max over all of the atomic kind set where la_max=ni & ni=2*lppl-2
      !Thus need to find the max of lppl over atomic_kind_set and use to define lmax_global
      !We are using FUNC_AB so no need to cater for different ga_gb_functions
      lmax_global = 0
      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)
         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, nexp_ppl=lppl)
         lmax_global = MAX(lppl, lmax_global)
      END DO
      lmax_global = 2*lmax_global-2

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)

         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, alpha_ppl=alpha, nexp_ppl=lppl, cexp_ppl=cexp_ppl)

         IF (lppl <= 0) CYCLE

         ni = ncoset(2*lppl-2)
         ALLOCATE (pab(ni, 1))
         pab = 0._dp

         nthread = 1
         ithread = 0

         CALL reallocate(cores, 1, natom)
         npme = 0
         cores = 0

         ! prepare core function
         DO j = 1, lppl
            SELECT CASE (j)
            CASE (1)
               pab(1, 1) = cexp_ppl(1)
            CASE (2)
               n = coset(2, 0, 0)
               pab(n, 1) = cexp_ppl(2)
               n = coset(0, 2, 0)
               pab(n, 1) = cexp_ppl(2)
               n = coset(0, 0, 2)
               pab(n, 1) = cexp_ppl(2)
            CASE (3)
               n = coset(4, 0, 0)
               pab(n, 1) = cexp_ppl(3)
               n = coset(0, 4, 0)
               pab(n, 1) = cexp_ppl(3)
               n = coset(0, 0, 4)
               pab(n, 1) = cexp_ppl(3)
               n = coset(2, 2, 0)
               pab(n, 1) = 2._dp*cexp_ppl(3)
               n = coset(2, 0, 2)
               pab(n, 1) = 2._dp*cexp_ppl(3)
               n = coset(0, 2, 2)
               pab(n, 1) = 2._dp*cexp_ppl(3)
            CASE (4)
               n = coset(6, 0, 0)
               pab(n, 1) = cexp_ppl(4)
               n = coset(0, 6, 0)
               pab(n, 1) = cexp_ppl(4)
               n = coset(0, 0, 6)
               pab(n, 1) = cexp_ppl(4)
               n = coset(4, 2, 0)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(4, 0, 2)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 4, 0)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 0, 4)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(0, 4, 2)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(0, 2, 4)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 2, 2)
               pab(n, 1) = 6._dp*cexp_ppl(4)
            CASE DEFAULT
               CPABORT("")
            END SELECT
         END DO

         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
                  npme = npme+1
                  cores(npme) = iatom
               ENDIF
            ELSE
               npme = npme+1
               cores(npme) = iatom
            ENDIF
         END DO

         IF (npme .GT. 0) THEN
            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               subpatch_pattern = 0
               ni = 2*lppl-2
               radius = exp_radius_single(0, ni, pab(:, 1), alpha, eps_rho_rspace, 1.0_dp)
               CALL collocate_pgf_rspace(ni, alpha, 0, ra, pab(:, 1), rs_rho, cell, radius, 1.0_dp)
            END DO
         ENDIF

         DEALLOCATE (pab)

      END DO

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF

      CALL rs_pw_transfer(rs_rho, vppl%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL timestop(handle)

   END SUBROUTINE calculate_ppl_grid

! **************************************************************************************************
!> \brief Collocates the fitted lri density on a grid.
!> \param lri_rho_g ...
!> \param lri_rho_r ...
!> \param qs_env ...
!> \param lri_coef ...
!> \param total_rho ...
!> \param basis_type ...
!> \param exact_1c_terms ...
!> \param pmat replicated block diagonal density matrix (optional)
!> \param atomlist list of atoms to be included (optional)
!> \par History
!>      04.2013
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_lri_rho_elec(lri_rho_g, lri_rho_r, qs_env, &
                                     lri_coef, total_rho, basis_type, exact_1c_terms, pmat, atomlist)

      TYPE(pw_p_type), INTENT(INOUT)                     :: lri_rho_g, lri_rho_r
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: lri_coef
      REAL(KIND=dp), INTENT(OUT)                         :: total_rho
      CHARACTER(len=*), INTENT(IN)                       :: basis_type
      LOGICAL, INTENT(IN)                                :: exact_1c_terms
      TYPE(dbcsr_type), OPTIONAL                         :: pmat
      INTEGER, DIMENSION(:), OPTIONAL                    :: atomlist

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_lri_rho_elec', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, dir, group_size, handle, iatom, igrid_level, ikind, ipgf, iset, jpgf, &
         jset, lmax_global, m1, maxco, maxpgf, maxset, maxsgf, maxsgf_set, my_pos, na1, natom, &
         nb1, ncoa, ncob, nseta, offset, sgfa, sgfb
      INTEGER, DIMENSION(3)                              :: lb, location, tp, ub
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, la_max, la_min, mylmax, &
                                                            npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, map_consistent
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: map_it
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: map_it2
      REAL(KIND=dp)                                      :: eps_rho_rspace, rab2, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: aci
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, pab, sphi_a, work, zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: lri_basis_set, orb_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: mgrid_gspace, mgrid_rspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(realspace_grid_type), POINTER                 :: rs_grid

      NULLIFY (aci, atomic_kind_set, qs_kind_set, atom_list, cell, cube_info, &
               dft_control, first_sgfa, gridlevel_info, la_max, la_min, lri_basis_set, &
               mgrid_gspace, mgrid_rspace, npgfa, nsgfa, pab, &
               particle_set, pw_env, pw_pools, rs_grid, rs_rho, sphi_a, &
               work, zeta, mylmax)

      CALL timeset(routineN, handle)

      IF (exact_1c_terms) THEN
         CPASSERT(PRESENT(pmat))
      END IF

      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, &
                      cell=cell, particle_set=particle_set, &
                      pw_env=pw_env, &
                      dft_control=dft_control)

      cube_info => pw_env%cube_info
      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      map_consistent = dft_control%qs_control%map_consistent
      gridlevel_info => pw_env%gridlevel_info

      ! *** set up the pw multi-grids *** !
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env=pw_env, rs_grids=rs_rho, pw_pools=pw_pools)

      CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                               use_data=REALDATA3D, &
                               in_space=REALSPACE)

      CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                               use_data=COMPLEXDATA1D, &
                               in_space=RECIPROCALSPACE)

      ! *** set up the rs multi-grids *** !
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
      END DO

      !take maxco from the LRI basis set!
      CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                           maxco=maxco, basis_type=basis_type)

      ALLOCATE (pab(maxco, 1))
      offset = 0
      my_pos = mgrid_rspace(1)%pw%pw_grid%para%my_pos
      group_size = mgrid_rspace(1)%pw%pw_grid%para%group_size

      !Find lmax_global across atomic_kind_set
      lmax_global = 0
      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=lri_basis_set, basis_type=basis_type)
         CALL get_gto_basis_set(gto_basis_set=lri_basis_set, lmax=mylmax)
         lmax_global = MAX(MAXVAL(mylmax), lmax_global)
      END DO

      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=lri_basis_set, basis_type=basis_type)

         !Take the lri basis set here!
         CALL get_gto_basis_set(gto_basis_set=lri_basis_set, lmax=la_max, &
                                lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa, &
                                sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)

         DO iatom = 1, natom
            atom_a = atom_list(iatom)
            IF (PRESENT(ATOMLIST)) THEN
               IF (atomlist(atom_a) == 0) CYCLE
            END IF
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            aci => lri_coef(ikind)%acoef(iatom, :)

            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE (map_it(m1))
            DO iset = 1, nseta
               ! collocate this set locally?
               map_it = .FALSE.
               DO ipgf = 1, npgfa(iset)
                  igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                  rs_grid => rs_rho(igrid_level)%rs_grid
                  IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                     DO dir = 1, 3
                        ! bounds of local grid (i.e. removing the 'wings'), if periodic
                        tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                        tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                        IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                           lb(dir) = rs_grid%lb_local(dir)+rs_grid%desc%border
                           ub(dir) = rs_grid%ub_local(dir)-rs_grid%desc%border
                        ELSE
                           lb(dir) = rs_grid%lb_local(dir)
                           ub(dir) = rs_grid%ub_local(dir)
                        ENDIF
                        ! distributed grid, only map if it is local to the grid
                        location(dir) = tp(dir)+rs_grid%desc%lb(dir)
                     ENDDO
                     IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                         lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                         lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                        map_it(ipgf) = .TRUE.
                     ENDIF
                  ELSE
                     ! not distributed, just a round-robin distribution over the full set of CPUs
                     IF (MODULO(offset, group_size) == my_pos) map_it(ipgf) = .TRUE.
                  ENDIF
               END DO
               offset = offset+1

               IF (ANY(map_it(1:npgfa(iset)))) THEN
                  sgfa = first_sgfa(1, iset)
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  m1 = sgfa+nsgfa(iset)-1
                  ALLOCATE (work(nsgfa(iset), 1))
                  work(1:nsgfa(iset), 1) = aci(sgfa:m1)
                  pab = 0._dp

                  CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), 1.0_dp, lri_basis_set%sphi(1, sgfa), &
                             SIZE(lri_basis_set%sphi, 1), work(1, 1), SIZE(work, 1), 0.0_dp, pab(1, 1), &
                             SIZE(pab, 1))

                  DO ipgf = 1, npgfa(iset)
                     na1 = (ipgf-1)*ncoset(la_max(iset))
                     igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                     rs_grid => rs_rho(igrid_level)%rs_grid
                     IF (map_it(ipgf)) THEN
                        CALL collocate_pgf_product_rspace(la_max=la_max(iset), &
                                                          zeta=zeta(ipgf, iset), &
                                                          la_min=la_min(iset), &
                                                          lb_max=0, zetb=0.0_dp, lb_min=0, &
                                                          ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                                          rab2=0.0_dp, scale=1._dp, &
                                                          pab=pab, o1=na1, o2=0, &
                                                          rsgrid=rs_grid, cell=cell, &
                                                          cube_info=cube_info(igrid_level), &
                                                          eps_rho_rspace=eps_rho_rspace, &
                                                          ga_gb_function=FUNC_AB, &
                                                          map_consistent=map_consistent, &
                                                          lmax_global=lmax_global)
                     ENDIF
                  END DO
                  DEALLOCATE (work)
               END IF
            END DO
            DEALLOCATE (map_it)
         END DO
      END DO

      DEALLOCATE (pab)

      ! process the one-center terms
      IF (exact_1c_terms) THEN
         ! find maximum numbers
         maxset = 0
         maxpgf = 0
         lmax_global = 0
         offset = 0
         DO ikind = 1, SIZE(qs_kind_set)
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind=qs_kind, &
                             basis_set=orb_basis_set, basis_type="ORB")
            IF (.NOT. ASSOCIATED(orb_basis_set)) CYCLE
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                   npgf=npgfa, nset=nseta, lmax=mylmax)
            maxset = MAX(nseta, maxset)
            maxpgf = MAX(MAXVAL(npgfa), maxpgf)
            lmax_global = MAX(MAXVAL(mylmax), lmax_global)
         END DO
         CALL get_qs_kind_set(qs_kind_set=qs_kind_set, &
                              maxco=maxco, &
                              maxsgf=maxsgf, &
                              maxsgf_set=maxsgf_set, &
                              basis_type="ORB")
         ALLOCATE (pab(maxco, maxco), work(maxco, maxsgf_set))

         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type="ORB")
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, lmax=la_max, &
                                   lmin=la_min, zet=zeta, nset=nseta, npgf=npgfa, &
                                   sphi=sphi_a, first_sgf=first_sgfa, nsgf_set=nsgfa)
            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               rab(:) = 0.0_dp
               rab2 = 0.0_dp
               CALL dbcsr_get_block_p(matrix=pmat, row=atom_a, col=atom_a, BLOCK=p_block, found=found)
               m1 = MAXVAL(npgfa(1:nseta))
               ALLOCATE (map_it2(m1, m1))
               DO iset = 1, nseta
                  DO jset = 1, nseta
                     ! processor mappint
                     map_it2 = .FALSE.
                     DO ipgf = 1, npgfa(iset)
                        DO jpgf = 1, npgfa(jset)
                           zetp = zeta(ipgf, iset)+zeta(jpgf, jset)
                           igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                           rs_grid => rs_rho(igrid_level)%rs_grid
                           IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                              DO dir = 1, 3
                                 ! bounds of local grid (i.e. removing the 'wings'), if periodic
                                 tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                                 tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                                 IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                                    lb(dir) = rs_grid%lb_local(dir)+rs_grid%desc%border
                                    ub(dir) = rs_grid%ub_local(dir)-rs_grid%desc%border
                                 ELSE
                                    lb(dir) = rs_grid%lb_local(dir)
                                    ub(dir) = rs_grid%ub_local(dir)
                                 ENDIF
                                 ! distributed grid, only map if it is local to the grid
                                 location(dir) = tp(dir)+rs_grid%desc%lb(dir)
                              ENDDO
                              IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                                  lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                                  lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                                 map_it2(ipgf, jpgf) = .TRUE.
                              ENDIF
                           ELSE
                              ! not distributed, just a round-robin distribution over the full set of CPUs
                              IF (MODULO(offset, group_size) == my_pos) map_it2(ipgf, jpgf) = .TRUE.
                           ENDIF
                        END DO
                     END DO
                     offset = offset+1
                     !
                     IF (ANY(map_it2(1:npgfa(iset), 1:npgfa(jset)))) THEN
                        ncoa = npgfa(iset)*ncoset(la_max(iset))
                        sgfa = first_sgfa(1, iset)
                        ncob = npgfa(jset)*ncoset(la_max(jset))
                        sgfb = first_sgfa(1, jset)
                        ! decontract density block
                        CALL dgemm("N", "N", ncoa, nsgfa(jset), nsgfa(iset), &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   p_block(sgfa, sgfb), SIZE(p_block, 1), &
                                   0.0_dp, work(1, 1), maxco)
                        CALL dgemm("N", "T", ncoa, ncob, nsgfa(jset), &
                                   1.0_dp, work(1, 1), maxco, &
                                   sphi_a(1, sgfb), SIZE(sphi_a, 1), &
                                   0.0_dp, pab(1, 1), maxco)
                        DO ipgf = 1, npgfa(iset)
                           DO jpgf = 1, npgfa(jset)
                              zetp = zeta(ipgf, iset)+zeta(jpgf, jset)
                              igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                              rs_grid => rs_rho(igrid_level)%rs_grid

                              na1 = (ipgf-1)*ncoset(la_max(iset))
                              nb1 = (jpgf-1)*ncoset(la_max(jset))

                              IF (map_it2(ipgf, jpgf)) THEN
                                 CALL collocate_pgf_product_rspace( &
                                    la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                    la_max(jset), zeta(jpgf, jset), la_min(jset), &
                                    ra, rab, rab2, 1.0_dp, pab, na1, nb1, &
                                    rs_grid, cell, cube_info(igrid_level), &
                                    eps_rho_rspace, ga_gb_function=FUNC_AB, &
                                    map_consistent=map_consistent, lmax_global=lmax_global)
                              END IF
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
               DEALLOCATE (map_it2)
               !
            END DO
         END DO
         DEALLOCATE (pab, work)
      END IF

      CALL pw_zero(lri_rho_g%pw)
      CALL pw_zero(lri_rho_r%pw)

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_zero(mgrid_rspace(igrid_level)%pw)
         CALL rs_pw_transfer(rs=rs_rho(igrid_level)%rs_grid, &
                             pw=mgrid_rspace(igrid_level)%pw, dir=rs2pw)
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid)
      END DO

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_zero(mgrid_gspace(igrid_level)%pw)
         CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                          mgrid_gspace(igrid_level)%pw)
         CALL pw_axpy(mgrid_gspace(igrid_level)%pw, lri_rho_g%pw)
      END DO
      CALL pw_transfer(lri_rho_g%pw, lri_rho_r%pw)
      total_rho = pw_integrate_function(lri_rho_r%pw, isign=-1)

      ! *** give back the multi-grids *** !
      CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)

      CALL timestop(handle)

   END SUBROUTINE calculate_lri_rho_elec

! **************************************************************************************************
!> \brief computes the density of the core charges on the grid
!> \param rho_core ...
!> \param total_rho ...
!> \param qs_env ...
!> \param only_nopaw ...
! **************************************************************************************************
   SUBROUTINE calculate_rho_core(rho_core, total_rho, qs_env, only_nopaw)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_core
      REAL(KIND=dp), INTENT(OUT)                         :: total_rho
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: only_nopaw

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_core', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, ikind, j, natom, &
                                                            npme
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: my_only_nopaw, paw_atom
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(1)                        :: cab
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type)                               :: cube_info
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, atomic_kind_set, qs_kind_set, particle_set, &
               atom_list, pw_env, rs_rho, auxbas_pw_pool, cores)

      my_only_nopaw = .FALSE.
      IF (PRESENT(only_nopaw)) my_only_nopaw = only_nopaw

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      cube_info = pw_env%cube_info(1)
      ! be careful in parallel nsmax is choosen with multigrid in mind!
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom, &
                          alpha_core_charge=alpha, ccore_charge=cab(1))

         IF (my_only_nopaw .AND. paw_atom) CYCLE
         IF (alpha == 0.0_dp .OR. cab(1) == 0.0_dp) CYCLE

         CALL reallocate(cores, 1, natom)
         npme = 0
         cores = 0

         DO iatom = 1, natom
            IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
                  npme = npme+1
                  cores(npme) = iatom
               ENDIF
            ELSE
               npme = npme+1
               cores(npme) = iatom
            ENDIF
         END DO

         IF (npme .GT. 0) THEN
            DO j = 1, npme
               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               radius = exp_radius_single(0, 0, cab, alpha, eps_rho_rspace, -1.0_dp)
               CALL collocate_pgf_rspace(0, alpha, 0, ra, cab, rs_rho, cell, radius, -1.0_dp)
            END DO
         ENDIF

      END DO

      IF (ASSOCIATED(cores)) THEN
         DEALLOCATE (cores)
      END IF

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      total_rho = pw_integrate_function(rhoc_r%pw, isign=-1)

      CALL pw_transfer(rhoc_r%pw, rho_core%pw)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_core

! **************************************************************************************************
!> \brief collocate a single Gaussian on the grid
!> \param rho_gb charge density generated by a single gaussian
!> \param qs_env qs environment
!> \param iatom_in atom index
!> \par History
!>        12.2011 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_single_gaussian(rho_gb, qs_env, iatom_in)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_gb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom_in

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_single_gaussian', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, npme
      REAL(KIND=dp)                                      :: eps_rho_rspace, eta, radius
      REAL(KIND=dp), DIMENSION(1)                        :: cab
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pw_env, rs_rho, auxbas_pw_pool)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      cab(1) = 1.0_dp
      iatom = iatom_in

      npme = 0

      IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
         IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
            npme = npme+1
         ENDIF
      ELSE
         npme = npme+1
      ENDIF

      IF (npme .GT. 0) THEN
         atom_a = qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
         ra(:) = pbc(qs_env%qmmm_env_qm%image_charge_pot%particles_all(atom_a)%r, cell)
         eta = qs_env%qmmm_env_qm%image_charge_pot%eta
         radius = exp_radius_single(0, 0, cab, eta, eps_rho_rspace, 1.0_dp)
         CALL collocate_pgf_rspace(0, eta, 0, ra, cab, rs_rho, cell, radius, 1.0_dp)
      ENDIF

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL pw_transfer(rhoc_r%pw, rho_gb%pw)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_single_gaussian

! **************************************************************************************************
!> \brief computes the image charge density on the grid (including coeffcients)
!> \param rho_metal image charge density
!> \param coeff expansion coefficients of the image charge density, i.e.
!>        rho_metal=sum_a c_a*g_a
!> \param total_rho_metal total induced image charge density
!> \param qs_env qs environment
!> \par History
!>        01.2012 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_metal(rho_metal, coeff, total_rho_metal, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_metal
      REAL(KIND=dp), DIMENSION(:), POINTER               :: coeff
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: total_rho_metal
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_metal', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, handle, iatom, j, natom, npme
      INTEGER, DIMENSION(:), POINTER                     :: cores
      REAL(KIND=dp)                                      :: eps_rho_rspace, eta, radius
      REAL(KIND=dp), DIMENSION(1)                        :: cab
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, dft_control, pw_env, rs_rho, auxbas_pw_pool, cores)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      cab(1) = 1.0_dp

      natom = SIZE(qs_env%qmmm_env_qm%image_charge_pot%image_mm_list)

      CALL reallocate(cores, 1, natom)
      npme = 0
      cores = 0

      DO iatom = 1, natom
         IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
               npme = npme+1
               cores(npme) = iatom
            ENDIF
         ELSE
            npme = npme+1
            cores(npme) = iatom
         ENDIF
      ENDDO

      IF (npme .GT. 0) THEN
         DO j = 1, npme
            iatom = cores(j)
            atom_a = qs_env%qmmm_env_qm%image_charge_pot%image_mm_list(iatom)
            ra(:) = pbc(qs_env%qmmm_env_qm%image_charge_pot%particles_all(atom_a)%r, cell)
            eta = qs_env%qmmm_env_qm%image_charge_pot%eta
            radius = exp_radius_single(0, 0, cab, eta, eps_rho_rspace, coeff(iatom))
            CALL collocate_pgf_rspace(0, eta, 0, ra, cab, rs_rho, cell, radius, coeff(iatom))
         ENDDO
      ENDIF

      DEALLOCATE (cores)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      IF (PRESENT(total_rho_metal)) &
         !minus sign: account for the fact that rho_metal has opposite sign
         total_rho_metal = pw_integrate_function(rhoc_r%pw, isign=-1)

      CALL pw_transfer(rhoc_r%pw, rho_metal%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_metal

! **************************************************************************************************
!> \brief collocate a single Gaussian on the grid for periodic RESP fitting
!> \param rho_gb charge density generated by a single gaussian
!> \param qs_env qs environment
!> \param eta width of single Gaussian
!> \param iatom_in atom index
!> \par History
!>        06.2012 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_resp_single(rho_gb, qs_env, eta, iatom_in)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_gb
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: eta
      INTEGER, INTENT(IN)                                :: iatom_in

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_resp_single', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, npme
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(1)                        :: cab
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pw_env, rs_rho, auxbas_pw_pool, particle_set)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      cab(1) = 1.0_dp
      iatom = iatom_in

      npme = 0

      IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
         IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
            npme = npme+1
         ENDIF
      ELSE
         npme = npme+1
      ENDIF

      IF (npme .GT. 0) THEN
         ra(:) = pbc(particle_set(iatom)%r, cell)
         radius = exp_radius_single(0, 0, cab, eta, eps_rho_rspace, 1.0_dp)
         CALL collocate_pgf_rspace(0, eta, 0, ra, cab, rs_rho, cell, radius, 1.0_dp)
      ENDIF

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL pw_transfer(rhoc_r%pw, rho_gb%pw)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_resp_single

! **************************************************************************************************
!> \brief computes the RESP charge density on a grid based on the RESP charges
!> \param rho_resp RESP charge density
!> \param coeff RESP charges, take care of normalization factor
!>        (eta/pi)**1.5 later
!> \param natom number of atoms
!> \param eta width of single Gaussian
!> \param qs_env qs environment
!> \par History
!>        01.2012 created
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE calculate_rho_resp_all(rho_resp, coeff, natom, eta, qs_env)

      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_resp
      REAL(KIND=dp), DIMENSION(:), POINTER               :: coeff
      INTEGER, INTENT(IN)                                :: natom
      REAL(KIND=dp), INTENT(IN)                          :: eta
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_resp_all', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, j, npme
      INTEGER, DIMENSION(:), POINTER                     :: cores
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(1)                        :: cab
      REAL(KIND=dp), DIMENSION(3)                        :: ra
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoc_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(realspace_grid_type), POINTER                 :: rs_rho

      CALL timeset(routineN, handle)

      NULLIFY (cell, cores, dft_control, pw_env, rs_rho, auxbas_pw_pool, &
               particle_set)

      CALL get_qs_env(qs_env=qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_rs_grid=rs_rho, &
                      auxbas_pw_pool=auxbas_pw_pool)
      CALL rs_grid_retain(rs_rho)
      CALL rs_grid_zero(rs_rho)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      cab(1) = 1.0_dp

      CALL reallocate(cores, 1, natom)
      npme = 0
      cores = 0

      DO iatom = 1, natom
         IF (rs_rho%desc%parallel .AND. .NOT. rs_rho%desc%distributed) THEN
            IF (MODULO(iatom, rs_rho%desc%group_size) == rs_rho%desc%my_pos) THEN
               npme = npme+1
               cores(npme) = iatom
            ENDIF
         ELSE
            npme = npme+1
            cores(npme) = iatom
         ENDIF
      ENDDO

      IF (npme .GT. 0) THEN
         DO j = 1, npme
            iatom = cores(j)
            ra(:) = pbc(particle_set(iatom)%r, cell)
            radius = exp_radius_single(0, 0, cab, eta, eps_rho_rspace, coeff(iatom))
            CALL collocate_pgf_rspace(0, eta, 0, ra, cab, rs_rho, cell, radius, coeff(iatom))
         ENDDO
      ENDIF

      DEALLOCATE (cores)

      CALL pw_pool_create_pw(auxbas_pw_pool, rhoc_r%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)

      CALL rs_pw_transfer(rs_rho, rhoc_r%pw, rs2pw)
      CALL rs_grid_release(rs_rho)

      CALL pw_transfer(rhoc_r%pw, rho_resp%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoc_r%pw)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_resp_all

! **************************************************************************************************
!> \brief maps a given wavefunction on the grid
!> \param mo_vectors ...
!> \param ivector ...
!> \param rho ...
!> \param rho_gspace ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param cell ...
!> \param dft_control ...
!> \param particle_set ...
!> \param pw_env ...
!> \param basis_type ...
!> \param external_vector ...
!> \par History
!>      08.2002 created [Joost VandeVondele]
!>      03.2006 made independent of qs_env [Joost VandeVondele]
!> \note
!>      modified calculate_rho_elec, should write the wavefunction represented by
!>      it's presumably dominated by the FFT and the rs->pw and back routines
!>
!>    BUGS ???
!>      does it take correctly the periodic images of the basis functions into account
!>      i.e. is only correct if the basis functions are localised enough to be just in 1 cell ?
! **************************************************************************************************
   SUBROUTINE calculate_wavefunction(mo_vectors, ivector, rho, rho_gspace, &
                                     atomic_kind_set, qs_kind_set, cell, dft_control, particle_set, &
                                     pw_env, basis_type, external_vector)

      TYPE(cp_fm_type), POINTER                          :: mo_vectors
      INTEGER                                            :: ivector
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho, rho_gspace
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: external_vector

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_wavefunction', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER :: dir, group, group_size, handle, i, iatom, igrid_level, ikind, ipgf, iset, &
         lmax_global, maxco, maxsgf_set, my_pos, na1, na2, nao, natom, ncoa, nseta, offset, sgfa
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: where_is_the_point
      INTEGER, DIMENSION(3)                              :: lb, location, tp, ub
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, mylmax, npgfa, nsgfa
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: local, map_it_here
      REAL(KIND=dp)                                      :: dab, eps_rho_rspace, rab2, scale, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvector
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab, sphi_a, work, zeta
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: mgrid_gspace, mgrid_rspace
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(REALSPACE_GRID_P_TYPE), DIMENSION(:), POINTER :: rs_rho

      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      CALL timeset(routineN, handle)

      NULLIFY (eigenvector, orb_basis_set, &
               pab, work, la_max, la_min, &
               npgfa, nsgfa, &
               sphi_a, zeta, first_sgfa, &
               rs_rho, pw_pools, mgrid_rspace, mgrid_gspace, mylmax)

      IF (PRESENT(external_vector)) THEN
         nao = SIZE(external_vector)
         ALLOCATE (eigenvector(nao))
         eigenvector = external_vector
      ELSE
         CALL cp_fm_get_info(matrix=mo_vectors, nrow_global=nao)
         ALLOCATE (eigenvector(nao))
         DO i = 1, nao
            CALL cp_fm_get_element(mo_vectors, i, ivector, eigenvector(i), local)
         ENDDO
      ENDIF

      ! *** set up the pw multi-grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_grids=rs_rho, pw_pools=pw_pools, &
                      cube_info=cube_info, gridlevel_info=gridlevel_info)

      CALL pw_pools_create_pws(pw_pools, mgrid_gspace, &
                               use_data=COMPLEXDATA1D, &
                               in_space=RECIPROCALSPACE)
      CALL pw_pools_create_pws(pw_pools, mgrid_rspace, &
                               use_data=REALDATA3D, &
                               in_space=REALSPACE)

      ! *** set up rs multi-grids
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
      END DO

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
!   *** Allocate work storage ***
      CALL get_atomic_kind_set(atomic_kind_set, natom=natom)
      CALL get_qs_kind_set(qs_kind_set, &
                           maxco=maxco, &
                           maxsgf_set=maxsgf_set, &
                           basis_type=my_basis_type)

      ALLOCATE (pab(maxco, 1))
      ALLOCATE (work(maxco, 1))

      offset = 0
      group = mgrid_rspace(1)%pw%pw_grid%para%group
      my_pos = mgrid_rspace(1)%pw%pw_grid%para%my_pos
      group_size = mgrid_rspace(1)%pw%pw_grid%para%group_size
      ALLOCATE (where_is_the_point(0:group_size-1))

      !loop over natom to find maximum value of la_max
      lmax_global = 0
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=my_basis_type)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, lmax=mylmax)
         lmax_global = MAX(MAXVAL(mylmax), lmax_global)
      END DO

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=my_basis_type)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                first_sgf=first_sgfa, &
                                lmax=la_max, &
                                lmin=la_min, &
                                npgf=npgfa, &
                                nset=nseta, &
                                nsgf_set=nsgfa, &
                                sphi=sphi_a, &
                                zet=zeta)
         ra(:) = pbc(particle_set(iatom)%r, cell)
         rab(:) = 0.0_dp
         rab2 = 0.0_dp
         dab = 0.0_dp

         DO iset = 1, nseta

            ncoa = npgfa(iset)*ncoset(la_max(iset))
            sgfa = first_sgfa(1, iset)

            DO i = 1, nsgfa(iset)
               work(i, 1) = eigenvector(offset+i)
            ENDDO

            CALL dgemm("N", "N", ncoa, 1, nsgfa(iset), &
                       1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                       work(1, 1), SIZE(work, 1), &
                       0.0_dp, pab(1, 1), SIZE(pab, 1))

            DO ipgf = 1, npgfa(iset)

               na1 = (ipgf-1)*ncoset(la_max(iset))+1
               na2 = ipgf*ncoset(la_max(iset))

               scale = 1.0_dp
               zetp = zeta(ipgf, iset)
               igrid_level = gaussian_gridlevel(gridlevel_info, zetp)

               map_it_here = .FALSE.

               IF (.NOT. ALL(rs_rho(igrid_level)%rs_grid%desc%perd == 1)) THEN
                  DO dir = 1, 3
                     ! bounds of local grid (i.e. removing the 'wings'), if periodic
                     tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_rho(igrid_level)%rs_grid%desc%npts(dir))
                     tp(dir) = MODULO(tp(dir), rs_rho(igrid_level)%rs_grid%desc%npts(dir))
                     IF (rs_rho(igrid_level)%rs_grid%desc%perd(dir) .NE. 1) THEN
                        lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local(dir)+rs_rho(igrid_level)%rs_grid%desc%border
                        ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local(dir)-rs_rho(igrid_level)%rs_grid%desc%border
                     ELSE
                        lb(dir) = rs_rho(igrid_level)%rs_grid%lb_local(dir)
                        ub(dir) = rs_rho(igrid_level)%rs_grid%ub_local(dir)
                     ENDIF
                     ! distributed grid, only map if it is local to the grid
                     location(dir) = tp(dir)+rs_rho(igrid_level)%rs_grid%desc%lb(dir)
                  ENDDO
                  IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                      lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                      lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                     map_it_here = .TRUE.
                  ENDIF
               ELSE
                  ! not distributed, just a round-robin distribution over the full set of CPUs
                  IF (MODULO(offset, group_size) == my_pos) map_it_here = .TRUE.
               ENDIF

               IF (map_it_here) THEN
                  CALL collocate_pgf_product_rspace( &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     0, 0.0_dp, 0, &
                     ra, rab, rab2, scale, pab, na1-1, 0, &
                     rs_rho(igrid_level)%rs_grid, cell, cube_info(igrid_level), &
                     eps_rho_rspace, map_consistent=.TRUE., ga_gb_function=FUNC_AB, &
                     lmax_global=lmax_global)
               END IF
            END DO

            offset = offset+nsgfa(iset)

         END DO

      END DO

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_pw_transfer(rs_rho(igrid_level)%rs_grid, &
                             mgrid_rspace(igrid_level)%pw, rs2pw)
         CALL rs_grid_release(rs_rho(igrid_level)%rs_grid)
      ENDDO

      CALL pw_zero(rho_gspace%pw)
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL pw_transfer(mgrid_rspace(igrid_level)%pw, &
                          mgrid_gspace(igrid_level)%pw)
         CALL pw_axpy(mgrid_gspace(igrid_level)%pw, rho_gspace%pw)
      END DO

      CALL pw_transfer(rho_gspace%pw, rho%pw)

      ! Release work storage
      DEALLOCATE (eigenvector)

      DEALLOCATE (pab)

      DEALLOCATE (work)

      ! give back the pw multi-grids
      CALL pw_pools_give_back_pws(pw_pools, mgrid_gspace)
      CALL pw_pools_give_back_pws(pw_pools, mgrid_rspace)

      CALL timestop(handle)

   END SUBROUTINE calculate_wavefunction

! **************************************************************************************************
!> \brief computes the density corresponding to a given density matrix on the grid
!> \param matrix_p ...
!> \param matrix_p_kp ...
!> \param rho ...
!> \param rho_gspace ...
!> \param total_rho ...
!> \param ks_env ...
!> \param soft_valid ...
!> \param compute_tau ...
!> \param compute_grad ...
!> \param basis_type ...
!> \param der_type ...
!> \param idir ...
!> \param task_list_external ...
!> \param pw_env_external ...
!> \par History
!>      IAB (15-Feb-2010): Added OpenMP parallelisation to task loop
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2010 on behalf of the HECToR project
!> \note
!>      both rho and rho_gspace contain the new rho
!>      (in real and g-space respectively)
! **************************************************************************************************
   SUBROUTINE calculate_rho_elec(matrix_p, matrix_p_kp, rho, rho_gspace, total_rho, &
                                 ks_env, soft_valid, compute_tau, compute_grad, &
                                 basis_type, der_type, idir, task_list_external, pw_env_external)

      TYPE(dbcsr_type), OPTIONAL, POINTER                :: matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_p_kp
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho, rho_gspace
      REAL(KIND=dp), INTENT(OUT)                         :: total_rho
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: soft_valid, compute_tau, compute_grad
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: basis_type
      INTEGER, INTENT(IN), OPTIONAL                      :: der_type, idir
      TYPE(task_list_type), OPTIONAL, POINTER            :: task_list_external
      TYPE(pw_env_type), OPTIONAL, POINTER               :: pw_env_external

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rho_elec', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER                                            :: code_base, ga_gb_function, handle, &
                                                            iatom, igrid_level, im, img, &
                                                            lmax_global, my_der_type, my_idir, n1, &
                                                            n2, natoms, nimages
      INTEGER(kind=int_8), DIMENSION(:), POINTER         :: atom_pair_recv, atom_pair_send
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      LOGICAL                                            :: distributed_rs_grids, do_kp, &
                                                            map_consistent, my_compute_grad, &
                                                            my_compute_tau, my_soft
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: posat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_base_type), DIMENSION(:), POINTER        :: rs_rho_base
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_desc_p_type), DIMENSION(:), &
         POINTER                                         :: rs_descs
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(task_list_type), POINTER                      :: task_list, task_list_soft

      CALL timeset(routineN, handle)

      CPASSERT(PRESENT(matrix_p) .OR. PRESENT(matrix_p_kp))
      do_kp = PRESENT(matrix_p_kp)

      ! by default, the full density is calculated
      my_soft = .FALSE.
      IF (PRESENT(soft_valid)) my_soft = soft_valid

      ! by default, do not compute the kinetic energy density (tau)
      ! if compute_tau, all grids referencing to rho are actually tau
      IF (PRESENT(compute_tau)) THEN
         my_compute_tau = compute_tau
      ELSE
         my_compute_tau = .FALSE.
      ENDIF

      IF (PRESENT(compute_grad)) THEN
         my_compute_grad = compute_grad
      ELSE
         my_compute_grad = .FALSE.
      ENDIF

      my_der_type = 0
      my_idir = 0
      IF (PRESENT(der_type)) THEN
         my_der_type = der_type
         ga_gb_function = FUNC_DER(my_der_type)
      ELSE IF (my_compute_tau) THEN
         ga_gb_function = FUNC_DADB
      ELSE IF (my_compute_grad) THEN
         ga_gb_function = FUNC_DABpADB
      ELSE
         ga_gb_function = FUNC_AB
      ENDIF
      IF (PRESENT(idir)) my_idir = idir

      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      CALL get_ks_env(ks_env, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set)

      ! get the task lists
      SELECT CASE (my_basis_type)
      CASE ("ORB")
         CALL get_ks_env(ks_env, &
                         task_list=task_list, &
                         task_list_soft=task_list_soft)
      CASE ("AUX_FIT")
         CALL get_ks_env(ks_env, &
                         task_list_aux_fit=task_list, &
                         task_list_soft=task_list_soft)
      END SELECT
      IF (PRESENT(task_list_external)) task_list => task_list_external
      ! the soft option will not work with task_list_external!!!!!!!!!
      IF (my_soft) task_list => task_list_soft
      CPASSERT(ASSOCIATED(task_list))

      nimages = dft_control%nimages
      CPASSERT(nimages == 1 .OR. do_kp)

      ! *** assign from pw_env
      CALL get_ks_env(ks_env, pw_env=pw_env)
      IF (PRESENT(pw_env_external)) pw_env => pw_env_external
      gridlevel_info => pw_env%gridlevel_info
      cube_info => pw_env%cube_info

      ! iatom - ikind mapping
      natoms = SIZE(particle_set)
      ALLOCATE (kind_of(natoms))
      DO iatom = 1, natoms
         kind_of(iatom) = particle_set(iatom)%atomic_kind%kind_number
      END DO
      ! atomic positions
      ALLOCATE (posat(3, natoms))
      DO iatom = 1, natoms
         posat(:, iatom) = pbc(particle_set(iatom)%r, cell)
      END DO

      ! Optionally increment lmax_global depending on which ga_gb_function
      ! will be used for collocation
      lmax_global = task_list%lmax_all
      SELECT CASE (ga_gb_function)
      CASE (FUNC_DADB)
         lmax_global = lmax_global+1
      CASE (FUNC_ADBmDAB)
         lmax_global = lmax_global+1
      CASE (FUNC_DABpADB)
         lmax_global = lmax_global+1
      CASE (FUNC_DX, FUNC_DY, FUNC_DZ)
         lmax_global = lmax_global+1
      CASE (FUNC_DXDY, FUNC_DYDZ, FUNC_DZDX)
         lmax_global = lmax_global+2
      CASE (FUNC_DXDX, FUNC_DYDY, FUNC_DZDZ)
         lmax_global = lmax_global+2
      CASE (FUNC_ARDBmDARB)
         lmax_global = lmax_global+2
      CASE (FUNC_ARB)
         lmax_global = lmax_global+1
      CASE (FUNC_AB)
         lmax_global = lmax_global
      CASE DEFAULT
         CPABORT("")
      END SELECT

      ! *** set up the pw multi-grids
      CPASSERT(ASSOCIATED(pw_env))
      CALL pw_env_get(pw_env, rs_descs=rs_descs, rs_grids=rs_rho)
      ! *** set up the rs multi-grids
      CPASSERT(ASSOCIATED(rs_rho))
      CPASSERT(ASSOCIATED(rs_descs))
      distributed_rs_grids = .FALSE.

      DO igrid_level = 1, gridlevel_info%ngrid_levels
         CALL rs_grid_retain(rs_rho(igrid_level)%rs_grid)
         CALL rs_grid_zero(rs_rho(igrid_level)%rs_grid)
         IF (rs_rho(igrid_level)%rs_grid%desc%distributed) THEN
            distributed_rs_grids = .TRUE.
         ENDIF
      END DO

      map_consistent = dft_control%qs_control%map_consistent

      !   *** Initialize working density matrix ***
      ! distributed rs grids require a matrix that will be changed
      ! whereas this is not the case for replicated grids
      NULLIFY (deltap)
      CALL dbcsr_allocate_matrix_set(deltap, nimages)
      IF (distributed_rs_grids) THEN
         DO img = 1, nimages
            ALLOCATE (deltap(img)%matrix)
         END DO
         ! this matrix has no strict sparsity pattern in parallel
         ! deltap%sparsity_id=-1
         IF (do_kp) THEN
            DO img = 1, nimages
               CALL dbcsr_copy(deltap(img)%matrix, matrix_p_kp(img)%matrix, &
                               name="DeltaP")
            END DO
         ELSE
            CALL dbcsr_copy(deltap(1)%matrix, matrix_p, name="DeltaP")
         END IF
      ELSE
         IF (do_kp) THEN
            DO img = 1, nimages
               deltap(img)%matrix => matrix_p_kp(img)%matrix
            END DO
         ELSE
            deltap(1)%matrix => matrix_p
         END IF
      ENDIF

      ! distribute the matrix
      IF (distributed_rs_grids) THEN
         atom_pair_send => task_list%atom_pair_send
         atom_pair_recv => task_list%atom_pair_recv
         CALL rs_distribute_matrix(rs_descs, deltap, atom_pair_send, atom_pair_recv, &
                                   task_list%natoms, task_list%nimages, scatter=.TRUE.)
      ENDIF

      ! load grid info to task_list
      CPASSERT(gridlevel_info%ngrid_levels == task_list%ngrid_levels)
      DO igrid_level = 1, gridlevel_info%ngrid_levels
         task_list%cinfo(igrid_level)%dh(1:3, 1:3) = cube_info(igrid_level)%dh(1:3, 1:3)
         task_list%cinfo(igrid_level)%dh_inv(1:3, 1:3) = cube_info(igrid_level)%dh_inv(1:3, 1:3)
         task_list%cinfo(igrid_level)%npts(1:3) = rs_rho(igrid_level)%rs_grid%desc%npts(1:3)
         task_list%cinfo(igrid_level)%lb_grid(1:3) = rs_rho(igrid_level)%rs_grid%desc%lb(1:3)- &
                                                     rs_rho(igrid_level)%rs_grid%lb_local(1:3)
         task_list%cinfo(igrid_level)%perd(1:3) = rs_rho(igrid_level)%rs_grid%desc%perd(1:3)
         task_list%cinfo(igrid_level)%drmin = cube_info(igrid_level)%drmin
         task_list%cinfo(igrid_level)%max_radius = cube_info(igrid_level)%max_radius
         task_list%cinfo(igrid_level)%max_rad_ga = cube_info(igrid_level)%max_rad_ga
         IF (ASSOCIATED(task_list%cinfo(igrid_level)%lb_cube)) &
            DEALLOCATE (task_list%cinfo(igrid_level)%lb_cube)
         IF (ASSOCIATED(task_list%cinfo(igrid_level)%ub_cube)) &
            DEALLOCATE (task_list%cinfo(igrid_level)%ub_cube)
         IF (ASSOCIATED(task_list%cinfo(igrid_level)%sphere_bounds)) &
            DEALLOCATE (task_list%cinfo(igrid_level)%sphere_bounds)
         IF (ASSOCIATED(cube_info(igrid_level)%lb_cube)) THEN
            n1 = SIZE(cube_info(igrid_level)%lb_cube, 1)
            n2 = SIZE(cube_info(igrid_level)%lb_cube, 2)
            ALLOCATE (task_list%cinfo(igrid_level)%lb_cube(n1, n2))
            task_list%cinfo(igrid_level)%lb_cube(:, :) = cube_info(igrid_level)%lb_cube(:, :)
         END IF
         IF (ASSOCIATED(cube_info(igrid_level)%ub_cube)) THEN
            n1 = SIZE(cube_info(igrid_level)%ub_cube, 1)
            n2 = SIZE(cube_info(igrid_level)%ub_cube, 2)
            ALLOCATE (task_list%cinfo(igrid_level)%ub_cube(n1, n2))
            task_list%cinfo(igrid_level)%ub_cube(:, :) = cube_info(igrid_level)%ub_cube(:, :)
         END IF
         IF (ASSOCIATED(cube_info(igrid_level)%sphere_bounds)) THEN
            n1 = MAXVAL(cube_info(igrid_level)%sphere_bounds_count)
            n2 = SIZE(cube_info(igrid_level)%sphere_bounds, 1)
            ALLOCATE (task_list%cinfo(igrid_level)%sphere_bounds(n1, n2))
            task_list%cinfo(igrid_level)%sphere_bounds = 0
            DO im = 1, n2
               n1 = cube_info(igrid_level)%sphere_bounds_count(im)
               task_list%cinfo(igrid_level)%sphere_bounds(1:n1, im) = cube_info(igrid_level)%sphere_bounds(im)%p(1:n1)
            END DO
         END IF
      END DO

      code_base = 0
      IF (.NOT. map_consistent) code_base = 1
      IF (.NOT. (ga_gb_function == FUNC_AB .OR. ga_gb_function == FUNC_DADB)) code_base = 1
      SELECT CASE (code_base)
      CASE (0)
         ALLOCATE (rs_rho_base(gridlevel_info%ngrid_levels))
         DO igrid_level = 1, gridlevel_info%ngrid_levels
            rs_rho_base(igrid_level)%grid => rs_rho(igrid_level)%rs_grid%r(:, :, :)
         END DO
         CALL grid_collocate_ref(rs_rho_base, deltap, task_list, &
                                 my_compute_tau, kind_of, posat, &
                                 lmax_global)
         DEALLOCATE (rs_rho_base)
      CASE (1)
         CALL grid_collocate_mp(rs_rho, deltap, task_list, cube_info, cell, &
                                ga_gb_function, my_idir, map_consistent, kind_of, posat, &
                                lmax_global)
      CASE DEFAULT
         CPABORT("NA")
      END SELECT

      ! Release work storage
      DEALLOCATE (kind_of, posat)

      IF (distributed_rs_grids) THEN
         CALL dbcsr_deallocate_matrix_set(deltap)
      ELSE
         DO img = 1, nimages
            NULLIFY (deltap(img)%matrix)
         END DO
         DEALLOCATE (deltap)
      ENDIF

      CALL density_rs2pw(pw_env, rs_rho, rho, rho_gspace)

      total_rho = pw_integrate_function(rho%pw, isign=-1)

      CALL timestop(handle)

   END SUBROUTINE calculate_rho_elec

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions in g-space
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rab2 ...
!> \param scale ...
!> \param pab ...
!> \param na ...
!> \param nb ...
!> \param eps_rho_gspace ...
!> \param gsq_max ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_gspace(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, rab2, scale, pab, na, nb, &
                                           eps_rho_gspace, gsq_max, pw)

      ! NOTE: this routine is much slower than collocate_pgf_product_rspace

      INTEGER, INTENT(IN)                                :: la_max
      REAL(dp), INTENT(IN)                               :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(dp), INTENT(IN)                               :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: ra, rab
      REAL(dp), INTENT(IN)                               :: rab2, scale
      REAL(dp), DIMENSION(:, :), POINTER                 :: pab
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(dp), INTENT(IN)                               :: eps_rho_gspace, gsq_max
      TYPE(pw_type), POINTER                             :: pw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'collocate_pgf_product_gspace', &
         routineP = moduleN//':'//routineN

      COMPLEX(dp), DIMENSION(3)                          :: phasefactor
      COMPLEX(dp), DIMENSION(:), POINTER                 :: rag, rbg
      COMPLEX(dp), DIMENSION(:, :, :, :), POINTER        :: cubeaxis
      INTEGER                                            :: ax, ay, az, bx, by, bz, handle, i, ico, &
                                                            ig, ig2, jco, jg, kg, la, lb, &
                                                            lb_cube_min, lb_grid, ub_cube_max, &
                                                            ub_grid
      INTEGER, DIMENSION(3)                              :: lb_cube, ub_cube
      REAL(dp)                                           :: f, fa, fb, pij, prefactor, rzetp, &
                                                            twozetp, zetp
      REAL(dp), DIMENSION(3)                             :: dg, expfactor, fap, fbp, rap, rbp, rp
      REAL(dp), DIMENSION(:), POINTER                    :: g

      CALL timeset(routineN, handle)

      dg(:) = twopi/(pw%pw_grid%npts(:)*pw%pw_grid%dr(:))

      zetp = zeta+zetb
      rzetp = 1.0_dp/zetp
      f = zetb*rzetp
      rap(:) = f*rab(:)
      rbp(:) = rap(:)-rab(:)
      rp(:) = ra(:)+rap(:)
      twozetp = 2.0_dp*zetp
      fap(:) = twozetp*rap(:)
      fbp(:) = twozetp*rbp(:)

      prefactor = scale*SQRT((pi*rzetp)**3)*EXP(-zeta*f*rab2)
      phasefactor(:) = EXP(CMPLX(0.0_dp, -rp(:)*dg(:), KIND=dp))
      expfactor(:) = EXP(-0.25*rzetp*dg(:)*dg(:))

      lb_cube(:) = pw%pw_grid%bounds(1, :)
      ub_cube(:) = pw%pw_grid%bounds(2, :)

      lb_cube_min = MINVAL(lb_cube(:))
      ub_cube_max = MAXVAL(ub_cube(:))

      NULLIFY (cubeaxis, g, rag, rbg)

      CALL reallocate(cubeaxis, lb_cube_min, ub_cube_max, 1, 3, 0, la_max, 0, lb_max)
      CALL reallocate(g, lb_cube_min, ub_cube_max)
      CALL reallocate(rag, lb_cube_min, ub_cube_max)
      CALL reallocate(rbg, lb_cube_min, ub_cube_max)

      lb_grid = LBOUND(pw%cc, 1)
      ub_grid = UBOUND(pw%cc, 1)

      DO i = 1, 3

         DO ig = lb_cube(i), ub_cube(i)
            ig2 = ig*ig
            cubeaxis(ig, i, 0, 0) = expfactor(i)**ig2*phasefactor(i)**ig
         END DO

         IF (la_max > 0) THEN
            DO ig = lb_cube(i), ub_cube(i)
               g(ig) = REAL(ig, dp)*dg(i)
               rag(ig) = CMPLX(fap(i), -g(ig), KIND=dp)
               cubeaxis(ig, i, 1, 0) = rag(ig)*cubeaxis(ig, i, 0, 0)
            END DO
            DO la = 2, la_max
               fa = REAL(la-1, dp)*twozetp
               DO ig = lb_cube(i), ub_cube(i)
                  cubeaxis(ig, i, la, 0) = rag(ig)*cubeaxis(ig, i, la-1, 0)+ &
                                           fa*cubeaxis(ig, i, la-2, 0)
               END DO
            END DO
            IF (lb_max > 0) THEN
               fa = twozetp
               DO ig = lb_cube(i), ub_cube(i)
                  rbg(ig) = CMPLX(fbp(i), -g(ig), KIND=dp)
                  cubeaxis(ig, i, 0, 1) = rbg(ig)*cubeaxis(ig, i, 0, 0)
                  cubeaxis(ig, i, 1, 1) = rbg(ig)*cubeaxis(ig, i, 1, 0)+ &
                                          fa*cubeaxis(ig, i, 0, 0)
               END DO
               DO lb = 2, lb_max
                  fb = REAL(lb-1, dp)*twozetp
                  DO ig = lb_cube(i), ub_cube(i)
                     cubeaxis(ig, i, 0, lb) = rbg(ig)*cubeaxis(ig, i, 0, lb-1)+ &
                                              fb*cubeaxis(ig, i, 0, lb-2)
                     cubeaxis(ig, i, 1, lb) = rbg(ig)*cubeaxis(ig, i, 1, lb-1)+ &
                                              fb*cubeaxis(ig, i, 1, lb-2)+ &
                                              fa*cubeaxis(ig, i, 0, lb-1)
                  END DO
               END DO
               DO la = 2, la_max
                  fa = REAL(la, dp)*twozetp
                  DO ig = lb_cube(i), ub_cube(i)
                     cubeaxis(ig, i, la, 1) = rbg(ig)*cubeaxis(ig, i, la, 0)+ &
                                              fa*cubeaxis(ig, i, la-1, 0)
                  END DO
                  DO lb = 2, lb_max
                     fb = REAL(lb-1, dp)*twozetp
                     DO ig = lb_cube(i), ub_cube(i)
                        cubeaxis(ig, i, la, lb) = rbg(ig)*cubeaxis(ig, i, la, lb-1)+ &
                                                  fb*cubeaxis(ig, i, la, lb-2)+ &
                                                  fa*cubeaxis(ig, i, la-1, lb-1)
                     END DO
                  END DO
               END DO
            END IF
         ELSE
            IF (lb_max > 0) THEN
               DO ig = lb_cube(i), ub_cube(i)
                  g(ig) = REAL(ig, dp)*dg(i)
                  rbg(ig) = CMPLX(fbp(i), -g(ig), KIND=dp)
                  cubeaxis(ig, i, 0, 1) = rbg(ig)*cubeaxis(ig, i, 0, 0)
               END DO
               DO lb = 2, lb_max
                  fb = REAL(lb-1, dp)*twozetp
                  DO ig = lb_cube(i), ub_cube(i)
                     cubeaxis(ig, i, 0, lb) = rbg(ig)*cubeaxis(ig, i, 0, lb-1)+ &
                                              fb*cubeaxis(ig, i, 0, lb-2)
                  END DO
               END DO
            END IF
         END IF

      END DO

      DO la = 0, la_max
         DO lb = 0, lb_max
            IF (la+lb == 0) CYCLE
            fa = (1.0_dp/twozetp)**(la+lb)
            DO i = 1, 3
               DO ig = lb_cube(i), ub_cube(i)
                  cubeaxis(ig, i, la, lb) = fa*cubeaxis(ig, i, la, lb)
               END DO
            END DO
         END DO
      END DO

      ! Add the current primitive Gaussian function product to grid

      DO ico = ncoset(la_min-1)+1, ncoset(la_max)

         ax = indco(1, ico)
         ay = indco(2, ico)
         az = indco(3, ico)

         DO jco = ncoset(lb_min-1)+1, ncoset(lb_max)

            pij = prefactor*pab(na+ico, nb+jco)

            IF (ABS(pij) < eps_rho_gspace) CYCLE

            bx = indco(1, jco)
            by = indco(2, jco)
            bz = indco(3, jco)

            DO i = lb_grid, ub_grid
               IF (pw%pw_grid%gsq(i) > gsq_max) CYCLE
               ig = pw%pw_grid%g_hat(1, i)
               jg = pw%pw_grid%g_hat(2, i)
               kg = pw%pw_grid%g_hat(3, i)
               pw%cc(i) = pw%cc(i)+pij*cubeaxis(ig, 1, ax, bx)* &
                          cubeaxis(jg, 2, ay, by)* &
                          cubeaxis(kg, 3, az, bz)
            END DO

         END DO

      END DO

      DEALLOCATE (cubeaxis)
      DEALLOCATE (g)
      DEALLOCATE (rag)
      DEALLOCATE (rbg)

      CALL timestop(handle)

   END SUBROUTINE collocate_pgf_product_gspace

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions
!>        this assumes the subpatch_collocate = T option
!>        i.e. no halos
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param ra ...
!> \param ca Vector of coefficients with length ncoset(la_max)
!> \param rsgrid ...
!> \param cell ...
!> \param radius ...
!> \param prefactor ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_rspace(la_max, zeta, la_min, ra, ca, rsgrid, cell, radius, prefactor)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra
      REAL(KIND=dp), DIMENSION(:)                        :: ca
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: radius, prefactor

      CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_rspace', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ico, l, lxa, lxyz, lya, lza, n1
      INTEGER, DIMENSION(2, 3)                           :: local_b
      INTEGER, DIMENSION(3)                              :: local_s, ng, periodic
      REAL(KIND=dp)                                      :: a
      REAL(KIND=dp), &
         DIMENSION((la_max+1)*(la_max+2)*(la_max+3)/6)   :: coef_xyz, poly_d3
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid

      a = MAXVAL(ABS(rsgrid%desc%dh))
      IF (radius .LT. a/2.0_dp) THEN
         RETURN
      END IF

      ng(:) = rsgrid%desc%npts(:)
      grid => rsgrid%r(:, :, :)

      n1 = ncoset(la_max)
      coef_xyz(1:n1) = 0.0_dp
      lxyz = 0
      DO lza = 0, la_max
      DO lya = 0, la_max-lza
      DO lxa = 0, la_max-lza-lya
         lxyz = lxyz+1
         l = lxa+lya+lza
         IF (l >= la_min) THEN
            ico = coset(lxa, lya, lza)
            coef_xyz(lxyz) = prefactor*ca(ico)
         END IF
      ENDDO
      ENDDO
      ENDDO

      periodic = 1 ! cell%perd
      CALL poly_cp2k2d3(coef_xyz, la_max, poly_d3)
      local_b(1, :) = rsgrid%lb_real-rsgrid%desc%lb
      local_b(2, :) = rsgrid%ub_real-rsgrid%desc%lb
      local_s = rsgrid%lb_real-rsgrid%lb_local
      CALL collocGauss(h=cell%hmat, h_inv=cell%h_inv, &
                       grid=grid, poly=poly_d3, alphai=zeta, posi=ra, max_r2=radius*radius, &
                       periodic=periodic, gdim=ng, local_bounds=local_b, local_shift=local_s)

   END SUBROUTINE collocate_pgf_rspace

END MODULE qs_collocate_density

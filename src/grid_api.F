!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Top-level API for the grid package.
! **************************************************************************************************
MODULE grid_api
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_get_info,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE gauss_colloc,                    ONLY: collocGauss
   USE grid_base_mp,                    ONLY: collocate_pgf_product_rspace,&
                                              grid_collocate_mp
   USE grid_base_ref,                   ONLY: grid_collocate_ref
   USE grid_ibase_mp,                   ONLY: grid_integrate_mp
   USE grid_ibase_ref,                  ONLY: grid_integrate_ref
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace,&
                                              integrate_pgf_rspace
   USE qs_modify_pab_block,             ONLY: &
        FUNC_AB, FUNC_ADBmDAB, FUNC_ARB, FUNC_ARDBmDARB, FUNC_DABpADB, FUNC_DADB, FUNC_DER, &
        FUNC_DX, FUNC_DXDX, FUNC_DXDY, FUNC_DY, FUNC_DYDY, FUNC_DYDZ, FUNC_DZ, FUNC_DZDX, FUNC_DZDZ
   USE realspace_grid_types,            ONLY: realspace_grid_p_type
   USE spmat_types,                     ONLY: spmat_allocate,&
                                              spmat_p_type,&
                                              spmat_put,&
                                              spmat_release,&
                                              spmat_type
   USE task_list_types,                 ONLY: grid_base_type,&
                                              task_list_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Possible backend to use for gpw
   INTEGER, PARAMETER :: grid_driver_openmp = 100
   INTEGER, PARAMETER :: grid_driver_reference = 200

   PUBLIC :: grid_driver_openmp, grid_driver_reference
   PUBLIC :: grid_collocate, grid_integrate

   ! from gauss_colloc
   PUBLIC :: collocGauss

   ! from grid_base_mp
   PUBLIC :: collocate_pgf_product_rspace

   ! from qs_integrate_potential_low
   PUBLIC :: integrate_pgf_product_rspace, integrate_pgf_rspace

   ! from qs_modify_pab_block
   PUBLIC :: FUNC_AB, FUNC_ADBmDAB, FUNC_ARB, FUNC_ARDBmDARB, FUNC_DABpADB, FUNC_DADB, FUNC_DER, &
             FUNC_DX, FUNC_DXDX, FUNC_DXDY, FUNC_DY, FUNC_DYDY, FUNC_DYDZ, FUNC_DZ, FUNC_DZDX, FUNC_DZDZ

! TODO: This is just to conveniently test both backend via the CI.
#if defined(_OPENMP)
   INTEGER, PARAMETER :: grid_driver_default = grid_driver_openmp
#else
   INTEGER, PARAMETER :: grid_driver_default = grid_driver_reference
#endif

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rs_rho ...
!> \param deltap ...
!> \param task_list ...
!> \param cube_info ...
!> \param cell ...
!> \param ga_gb_function ...
!> \param my_idir ...
!> \param map_consistent ...
!> \param compute_tau ...
!> \param kind_of ...
!> \param posat ...
!> \param lmax_global ...
! **************************************************************************************************
   SUBROUTINE grid_collocate(rs_rho, deltap, task_list, cube_info, cell, ga_gb_function, &
                             my_idir, map_consistent, compute_tau, kind_of, posat, lmax_global)

      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_rho
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(cell_type), POINTER                           :: cell
      INTEGER, INTENT(IN)                                :: ga_gb_function, my_idir
      LOGICAL, INTENT(IN)                                :: map_consistent, compute_tau
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      INTEGER, INTENT(IN)                                :: lmax_global

      CHARACTER(len=*), PARAMETER                        :: routineN = 'grid_collocate'

      INTEGER                                            :: driver, handle, igrid_level, img, &
                                                            my_pos, ngrid_levels, nimages, nrow
      TYPE(grid_base_type), DIMENSION(:), POINTER        :: rs_rho_base
      TYPE(spmat_p_type), ALLOCATABLE, DIMENSION(:)      :: spmats

      CALL timeset(routineN, handle)

      driver = grid_driver_default
      IF (.NOT. map_consistent .OR. .NOT. (ga_gb_function == FUNC_AB .OR. ga_gb_function == FUNC_DADB)) &
         driver = grid_driver_openmp

      SELECT CASE (driver)
      CASE (grid_driver_reference)
         ngrid_levels = SIZE(rs_rho)
         ALLOCATE (rs_rho_base(ngrid_levels))
         DO igrid_level = 1, ngrid_levels
            rs_rho_base(igrid_level)%grid => rs_rho(igrid_level)%rs_grid%r(:, :, :)
         END DO
         CALL dbcsr_get_info(deltap(1)%matrix, nblkrows_total=nrow)
         nimages = SIZE(deltap)
         ALLOCATE (spmats(nimages))
         DO img = 1, nimages
            CALL spmat_allocate(spmats(img)%spmat, nrow, nrow)
            CALL spmat_get_dbcsr(spmats(img)%spmat, deltap(img)%matrix)
         END DO
         my_pos = rs_rho(1)%rs_grid%desc%my_pos

         IF (.FALSE.) THEN
            CALL write_task_file(task_list, spmats, rs_rho_base, kind_of, posat, my_pos)
         END IF

         CALL grid_collocate_ref(rs_rho_base, spmats, task_list, compute_tau, kind_of, posat, lmax_global)

         DO img = 1, nimages
            CALL spmat_release(spmats(img)%spmat)
         END DO
         DEALLOCATE (spmats)
         DEALLOCATE (rs_rho_base)

      CASE (grid_driver_openmp)
         CALL grid_collocate_mp(rs_rho, deltap, task_list, cube_info, cell, ga_gb_function, &
                                my_idir, map_consistent, kind_of, posat, lmax_global)

      CASE DEFAULT
         CPABORT("Unknown grid driver.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE grid_collocate

! **************************************************************************************************
!> \brief ...
!> \param spmat ...
!> \param dbmat ...
! **************************************************************************************************
   SUBROUTINE spmat_get_dbcsr(spmat, dbmat)
      TYPE(spmat_type), INTENT(INOUT)                    :: spmat
      TYPE(dbcsr_type), INTENT(IN)                       :: dbmat

      INTEGER                                            :: blk, icol, irow
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dblock
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, dbmat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, dblock, blk)
         CALL spmat_put(spmat, irow, icol, dblock)
      ENDDO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE spmat_get_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param rs_v ...
!> \param dhmat ...
!> \param deltap ...
!> \param task_list ...
!> \param cube_info ...
!> \param cell ...
!> \param eps_gvg_rspace ...
!> \param admm_scal_fac ...
!> \param compute_tau ...
!> \param pab_required ...
!> \param map_consistent ...
!> \param kind_of ...
!> \param posat ...
!> \param calculate_forces ...
!> \param force ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE grid_integrate(rs_v, dhmat, deltap, task_list, cube_info, cell, eps_gvg_rspace, &
                             admm_scal_fac, compute_tau, pab_required, map_consistent, &
                             kind_of, posat, calculate_forces, force, use_virial, virial)

      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dhmat, deltap
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: eps_gvg_rspace, admm_scal_fac
      LOGICAL, INTENT(IN)                                :: compute_tau, pab_required, map_consistent
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: force
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(len=*), PARAMETER                        :: routineN = 'grid_integrate'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (grid_driver_default)
      CASE (grid_driver_reference)
         CALL grid_integrate_ref(rs_v, dhmat, deltap, task_list, cube_info, cell, eps_gvg_rspace, &
                                 admm_scal_fac, compute_tau, pab_required, map_consistent, &
                                 kind_of, posat, calculate_forces, force, use_virial, virial)

      CASE (grid_driver_openmp)
         CALL grid_integrate_mp(rs_v, dhmat, deltap, task_list, cube_info, cell, eps_gvg_rspace, &
                                admm_scal_fac, compute_tau, pab_required, map_consistent, &
                                kind_of, posat, calculate_forces, force, use_virial, virial)
      CASE DEFAULT
         CPABORT("Unknown grid driver.")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE grid_integrate

! **************************************************************************************************
!> \brief ...
!> \param task_list ...
!> \param spmats ...
!> \param rs_rho ...
!> \param kind_of ...
!> \param posat ...
!> \param mp ...
! **************************************************************************************************
   SUBROUTINE write_task_file(task_list, spmats, rs_rho, kind_of, posat, mp)
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(spmat_p_type), DIMENSION(:), INTENT(IN), &
         TARGET                                          :: spmats
      TYPE(grid_base_type), DIMENSION(:), INTENT(IN)     :: rs_rho
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      INTEGER, INTENT(IN)                                :: mp

      CHARACTER(LEN=default_string_length)               :: fname, proc
      INTEGER                                            :: ig, ir, iu, n1, n2, ng, nset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: irs
      INTEGER, DIMENSION(2, 3)                           :: gridbounds
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: mat
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      TYPE(spmat_type), POINTER                          :: spmat

      WRITE (proc, FMT="(I10)") mp
      fname = "tasklist_"//TRIM(ADJUSTL(proc))//".data"
      iu = 999
      OPEN (unit=iu, file=fname, position="rewind")
      !
      WRITE (iu, *) "S> TASK GENERAL"
      WRITE (iu, *) task_list%ntasks
      WRITE (iu, *) task_list%ngrid_levels
      WRITE (iu, *) task_list%nimages
      WRITE (iu, *) task_list%natoms
      WRITE (iu, *) task_list%maxset
      WRITE (iu, *) task_list%maxpgf
      WRITE (iu, *) task_list%lmax_all
      WRITE (iu, *) task_list%maxco
      WRITE (iu, *) task_list%maxsgf_set
      WRITE (iu, *) task_list%eps_rho_rspace
      WRITE (iu, *) "E> TASK GENERAL"
      !
      WRITE (iu, *) "S> TASK GRID"
      ng = SIZE(task_list%cinfo)
      WRITE (iu, *) ng
      DO ig = 1, ng
         WRITE (iu, *) task_list%cinfo(ig)%dh
         WRITE (iu, *) task_list%cinfo(ig)%dh_inv
         WRITE (iu, *) task_list%cinfo(ig)%perd
         WRITE (iu, *) task_list%cinfo(ig)%npts
         WRITE (iu, *) task_list%cinfo(ig)%lb_grid
         WRITE (iu, *) task_list%cinfo(ig)%drmin
         WRITE (iu, *) task_list%cinfo(ig)%max_radius
         WRITE (iu, *) task_list%cinfo(ig)%max_rad_ga
         IF (ASSOCIATED(task_list%cinfo(ig)%lb_cube)) THEN
            n1 = SIZE(task_list%cinfo(ig)%lb_cube, 1)
            n2 = SIZE(task_list%cinfo(ig)%lb_cube, 2)
            WRITE (iu, *) n1, n2
            WRITE (iu, *) task_list%cinfo(ig)%lb_cube
         ELSE
            n1 = 0
            n2 = 0
            WRITE (iu, *) n1, n2
         END IF
         IF (ASSOCIATED(task_list%cinfo(ig)%ub_cube)) THEN
            n1 = SIZE(task_list%cinfo(ig)%ub_cube, 1)
            n2 = SIZE(task_list%cinfo(ig)%ub_cube, 2)
            WRITE (iu, *) n1, n2
            WRITE (iu, *) task_list%cinfo(ig)%ub_cube
         ELSE
            n1 = 0
            n2 = 0
            WRITE (iu, *) n1, n2
         END IF
         IF (ASSOCIATED(task_list%cinfo(ig)%sphere_bounds)) THEN
            n1 = SIZE(task_list%cinfo(ig)%sphere_bounds, 1)
            n2 = SIZE(task_list%cinfo(ig)%sphere_bounds, 2)
            WRITE (iu, *) n1, n2
            WRITE (iu, *) task_list%cinfo(ig)%sphere_bounds
         ELSE
            n1 = 0
            n2 = 0
            WRITE (iu, *) n1, n2
         END IF
      END DO
      WRITE (iu, *) "E> TASK GRID"
      !
      WRITE (iu, *) "S> TASK BASIS"
      ng = SIZE(task_list%basis)
      WRITE (iu, *) ng
      DO ig = 1, ng
         IF (ASSOCIATED(task_list%basis(ig)%lmax)) THEN
            nset = SIZE(task_list%basis(ig)%lmax)
            WRITE (iu, *) nset
            WRITE (iu, *) task_list%basis(ig)%first_sgf
            WRITE (iu, *) task_list%basis(ig)%lmax
            WRITE (iu, *) task_list%basis(ig)%lmin
            WRITE (iu, *) task_list%basis(ig)%npgf
            WRITE (iu, *) task_list%basis(ig)%nsgf
            n1 = SIZE(task_list%basis(ig)%sphi, 1)
            n2 = SIZE(task_list%basis(ig)%sphi, 2)
            WRITE (iu, *) n1, n2
            WRITE (iu, *) task_list%basis(ig)%sphi
            n1 = SIZE(task_list%basis(ig)%zet, 1)
            n2 = SIZE(task_list%basis(ig)%zet, 2)
            WRITE (iu, *) n1, n2
            WRITE (iu, *) task_list%basis(ig)%zet
         ELSE
            nset = 0
            WRITE (iu, *) nset
         END IF
      END DO
      WRITE (iu, *) "E> TASK BASIS"
      !
      WRITE (iu, *) "S> TASK DIST"
      n1 = SIZE(task_list%dist_ab, 1)
      n2 = SIZE(task_list%dist_ab, 2)
      WRITE (iu, *) n1, n2
      WRITE (iu, *) task_list%dist_ab
      n1 = SIZE(task_list%taskstart, 1)
      n2 = SIZE(task_list%taskstart, 2)
      WRITE (iu, *) n1, n2
      WRITE (iu, *) task_list%taskstart
      WRITE (iu, *) task_list%taskstop
      n1 = SIZE(task_list%npairs, 1)
      WRITE (iu, *) n1
      WRITE (iu, *) task_list%npairs
      WRITE (iu, *) "E> TASK DIST"
      !
      CPABORT("TODO: serialize task_type routine not yet implemented")
      !WRITE (iu, *) "S> TASK LIST"
      !n1 = SIZE(task_list%tasks, 1)
      !n2 = SIZE(task_list%tasks, 2)
      !WRITE (iu, *) n1, n2
      !WRITE (iu, *) task_list%tasks
      !WRITE (iu, *) "E> TASK LIST"
      !
      ng = SIZE(spmats)
      WRITE (iu, *) "S> SPMAT"
      WRITE (iu, *) ng
      DO ig = 1, ng
         spmat => spmats(ig)%spmat
         WRITE (iu, *) spmat%nrow, spmat%ncol, spmat%nnz
         ALLOCATE (irs(spmat%nrow))
         irs = 0
         DO ir = 1, spmat%nrow
            IF (ALLOCATED(spmat%rind(ir)%col)) THEN
               irs(ir) = SIZE(spmat%rind(ir)%col)
            END IF
         END DO
         DO ir = 1, spmat%nrow
            WRITE (iu, *) irs(ir)
            IF (irs(ir) > 0) THEN
               WRITE (iu, *) spmat%rind(ir)%col(1:irs(ir))
               WRITE (iu, *) spmat%rind(ir)%blk(1:irs(ir))
            ELSE
               WRITE (iu, *)
               WRITE (iu, *)
            END IF
         END DO
         DEALLOCATE (irs)
         !data
         DO ir = 1, spmat%nnz
            mat => spmat%spdata(ir)%mat
            WRITE (iu, *) SIZE(mat, 1), SIZE(mat, 2)
            WRITE (iu, *) mat
         END DO
      END DO
      WRITE (iu, *) "E> SPMAT"
      !
      ng = SIZE(kind_of)
      WRITE (iu, *) "S> ATOMS"
      WRITE (iu, *) ng
      WRITE (iu, *) kind_of
      WRITE (iu, *) posat
      WRITE (iu, *) "E> ATOMS"
      !
      ng = SIZE(rs_rho)
      WRITE (iu, *) "S> RS GRID"
      WRITE (iu, *) ng
      DO ig = 1, ng
         grid => rs_rho(ig)%grid(:, :, :)
         gridbounds(1, 1) = LBOUND(grid, 1)
         gridbounds(2, 1) = UBOUND(grid, 1)
         gridbounds(1, 2) = LBOUND(grid, 2)
         gridbounds(2, 2) = UBOUND(grid, 2)
         gridbounds(1, 3) = LBOUND(grid, 3)
         gridbounds(2, 3) = UBOUND(grid, 3)
         WRITE (iu, *) gridbounds
      END DO
      WRITE (iu, *) "E> RS GRID"
      !
      CLOSE (unit=iu)

   END SUBROUTINE write_task_file

END MODULE grid_api

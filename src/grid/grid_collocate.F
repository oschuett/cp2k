!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE grid_collocate
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE realspace_grid_types,            ONLY: realspace_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_collocate'

   PUBLIC :: collocate_pgf_product

   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_AB = 100
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DADB = 200
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB_X = 301
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB_Y = 302
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ADBmDAB_Z = 303
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_XX = 411
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_XY = 412
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_XZ = 413
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_YX = 421
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_YY = 422
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_YZ = 423
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_ZX = 431
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_ZY = 432
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_ARDBmDARB_ZZ = 433
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB_X = 501
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB_Y = 502
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DABpADB_Z = 503
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DX = 601
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DY = 602
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZ = 603
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DXDY = 701
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DYDZ = 702
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZDX = 703
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DXDX = 801
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DYDY = 802
   INTEGER, PARAMETER, PUBLIC :: GRID_FUNC_DZDZ = 803

CONTAINS

! **************************************************************************************************
!> \brief low level collocation of primitive gaussian functions
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param scale ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param rsgrid ...
!> \param cell ...
!> \param cube_info ...
!> \param ga_gb_function ...
!> \param radius ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product(la_max, zeta, la_min, &
                                    lb_max, zetb, lb_min, &
                                    ra, rab, scale, pab, o1, o2, &
                                    rsgrid, cell, cube_info, &
                                    ga_gb_function, radius, &
                                    use_subpatch, subpatch_pattern)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rab
      REAL(KIND=dp), INTENT(IN)                          :: scale
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      TYPE(realspace_grid_type)                          :: rsgrid
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), INTENT(IN)                   :: cube_info
      INTEGER, INTENT(IN)                                :: ga_gb_function
      REAL(KIND=dp), INTENT(IN)                          :: radius
      LOGICAL, OPTIONAL                                  :: use_subpatch
      INTEGER(KIND=int_8), INTENT(IN), OPTIONAL          :: subpatch_pattern

      CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_product', &
         routineP = moduleN//':'//routineN

      INTEGER(KIND=int_8)                                :: my_subpatch_pattern
      INTEGER, DIMENSION(3), TARGET                      :: lb_grid, ngrid
      LOGICAL                                            :: my_use_subpatch, use_ortho
      LOGICAL(KIND=C_BOOL), DIMENSION(3), TARGET         :: periodic
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      INTERFACE
         SUBROUTINE grid_collocate_pgf_product_cpu(use_ortho, func, &
                                                   la_max, la_min, lb_max, lb_min, &
                                                   zeta, zetb, rscale, &
                                                   dh, dh_inv, ra, rab, npts, ngrid, lb_grid, &
                                                   periodic, radius, &
                                                   o1, o2, n1, n2, pab, &
                                                   grid) &
            BIND(C, name="grid_collocate_pgf_product_cpu")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: use_ortho
            INTEGER(KIND=C_INT), VALUE                :: func
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            REAL(KIND=C_DOUBLE), VALUE                :: zeta
            REAL(KIND=C_DOUBLE), VALUE                :: zetb
            REAL(KIND=C_DOUBLE), VALUE                :: rscale
            TYPE(C_PTR), VALUE                        :: dh
            TYPE(C_PTR), VALUE                        :: dh_inv
            TYPE(C_PTR), VALUE                        :: ra
            TYPE(C_PTR), VALUE                        :: rab
            TYPE(C_PTR), VALUE                        :: npts
            TYPE(C_PTR), VALUE                        :: ngrid
            TYPE(C_PTR), VALUE                        :: lb_grid
            TYPE(C_PTR), VALUE                        :: periodic
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: n1
            INTEGER(KIND=C_INT), VALUE                :: n2
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: grid
         END SUBROUTINE grid_collocate_pgf_product_cpu
      END INTERFACE

      !TODO remove unused dummy arguments
      MARK_USED(cell)
      MARK_USED(cube_info)

      IF (PRESENT(use_subpatch)) THEN
         my_use_subpatch = use_subpatch
      ELSE
         my_use_subpatch = .FALSE.
      ENDIF
      CPASSERT(.NOT. my_use_subpatch .OR. PRESENT(subpatch_pattern))

      my_subpatch_pattern = 0
      IF (PRESENT(subpatch_pattern)) THEN
         ! TODO: why is subpatch_pattern int_8, shouldn't 4 bytes suffice?
         my_subpatch_pattern = subpatch_pattern
      ENDIF

      IF (my_subpatch_pattern /= 0) THEN
         !TODO: implement, needed e.g. by QS/regtest-nmr-2/he2_bug_disp.inp
         CPABORT("non-zero subpatch_pattern not yet supported")
      ENDIF

      IF (2.0_dp*radius < MAXVAL(ABS(rsgrid%desc%dh))) THEN
         RETURN  ! radius too small to be mapped onto grid of given resolution
      END IF

      use_ortho = rsgrid%desc%orthorhombic .AND. .NOT. my_use_subpatch
      periodic = (rsgrid%desc%perd == 1)

      CPASSERT(LBOUND(pab, 1) == 1)
      CPASSERT(LBOUND(pab, 2) == 1)

      ! see rs_grid_create() in ./src/pw/realspace_grid_types.F
      !TODO: this might also be in the realspace_grid_type

      ! We absorb the local lower bounds into lb_grid.
      CPASSERT(LBOUND(rsgrid%r, 1) == rsgrid%lb_local(1))
      CPASSERT(LBOUND(rsgrid%r, 2) == rsgrid%lb_local(2))
      CPASSERT(LBOUND(rsgrid%r, 3) == rsgrid%lb_local(3))
      lb_grid = rsgrid%desc%lb - rsgrid%lb_local
      grid(1:, 1:, 1:) => rsgrid%r(:, :, :)
      ngrid = (/SIZE(grid, 1), SIZE(grid, 2), SIZE(grid, 3)/)

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh))
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh_inv))
      CPASSERT(IS_CONTIGUOUS(ra))
      CPASSERT(IS_CONTIGUOUS(rab))
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%npts))
      CPASSERT(IS_CONTIGUOUS(ngrid))
      CPASSERT(IS_CONTIGUOUS(lb_grid))
      CPASSERT(IS_CONTIGUOUS(periodic))
      CPASSERT(IS_CONTIGUOUS(pab))
      CPASSERT(IS_CONTIGUOUS(grid))
#endif

      CALL grid_collocate_pgf_product_cpu(use_ortho=LOGICAL(use_ortho, KIND=C_BOOL), &
                                          func=ga_gb_function, &
                                          la_max=la_max, &
                                          la_min=la_min, &
                                          lb_max=lb_max, &
                                          lb_min=lb_min, &
                                          zeta=zeta, &
                                          zetb=zetb, &
                                          rscale=scale, &
                                          dh=C_LOC(rsgrid%desc%dh(1, 1)), &
                                          dh_inv=C_LOC(rsgrid%desc%dh_inv(1, 1)), &
                                          ra=C_LOC(ra(1)), &
                                          rab=C_LOC(rab(1)), &
                                          npts=C_LOC(rsgrid%desc%npts(1)), &
                                          ngrid=C_LOC(ngrid(1)), &
                                          lb_grid=C_LOC(lb_grid(1)), &
                                          periodic=C_LOC(periodic(1)), &
                                          radius=radius, &
                                          o1=o1, &
                                          o2=o2, &
                                          n1=SIZE(pab, 1), &
                                          n2=SIZE(pab, 2), &
                                          pab=C_LOC(pab(1, 1)), &
                                          grid=C_LOC(grid(1, 1, 1)))

   END SUBROUTINE collocate_pgf_product

END MODULE grid_collocate

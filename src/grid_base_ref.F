!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE grid_base_ref
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE kinds,                           ONLY: dp
   USE spmat_types,                     ONLY: spmat_get,&
                                              spmat_p_type
   USE task_list_types,                 ONLY: grid_base_type,&
                                              task_list_type,&
                                              task_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

!&<
   INTEGER, PARAMETER :: maxfac = 30
   REAL(KIND=dp), PARAMETER, DIMENSION(0:maxfac) :: fac = (/ &
        0.10000000000000000000E+01_dp, 0.10000000000000000000E+01_dp, 0.20000000000000000000E+01_dp, &
        0.60000000000000000000E+01_dp, 0.24000000000000000000E+02_dp, 0.12000000000000000000E+03_dp, &
        0.72000000000000000000E+03_dp, 0.50400000000000000000E+04_dp, 0.40320000000000000000E+05_dp, &
        0.36288000000000000000E+06_dp, 0.36288000000000000000E+07_dp, 0.39916800000000000000E+08_dp, &
        0.47900160000000000000E+09_dp, 0.62270208000000000000E+10_dp, 0.87178291200000000000E+11_dp, &
        0.13076743680000000000E+13_dp, 0.20922789888000000000E+14_dp, 0.35568742809600000000E+15_dp, &
        0.64023737057280000000E+16_dp, 0.12164510040883200000E+18_dp, 0.24329020081766400000E+19_dp, &
        0.51090942171709440000E+20_dp, 0.11240007277776076800E+22_dp, 0.25852016738884976640E+23_dp, &
        0.62044840173323943936E+24_dp, 0.15511210043330985984E+26_dp, 0.40329146112660563558E+27_dp, &
        0.10888869450418352161E+29_dp, 0.30488834461171386050E+30_dp, 0.88417619937397019545E+31_dp, &
        0.26525285981219105864E+33_dp/)
   INTEGER, PARAMETER :: maxl = 18
   INTEGER, PARAMETER, DIMENSION(-1:maxl) :: ncoset = (/ &
        0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455, 560, 680, 816, 969, 1140, 1330/)
!&>

   PUBLIC :: grid_collocate_ref

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rs_rho Input/Output
!> \param spmats input
!> \param task_list input
!> \param compute_tau ...
!> \param kind_of input
!> \param posat input
!> \param lmax_global input
! **************************************************************************************************
   SUBROUTINE grid_collocate_ref(rs_rho, spmats, task_list, compute_tau, kind_of, posat, lmax_global)

      TYPE(grid_base_type), DIMENSION(:)                 :: rs_rho
      TYPE(spmat_p_type), DIMENSION(:), INTENT(IN)       :: spmats
      TYPE(task_list_type), POINTER                      :: task_list
      LOGICAL, INTENT(IN)                                :: compute_tau
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      INTEGER, INTENT(IN)                                :: lmax_global

      INTEGER :: bcol, brow, iatom, iatom_old, igrid_level, ikind, img, img_old, ipair, ipgf, &
         iset, iset_old, itask, jatom, jatom_old, jkind, jpgf, jset, jset_old, max_radius, maxco, &
         maxpgf, maxset, maxsgf_set, n1, n2, na1, natoms, nb1, ncoa, ncob, nimages, ntasks, sgfa, &
         sgfb
      INTEGER, DIMENSION(3)                              :: lb_grid(3), npts(3), perd(3)
      INTEGER, DIMENSION(:), POINTER                     :: first_sgfa, first_sgfb, la_max, la_min, &
                                                            lb_max, lb_min, npgfa, npgfb, nsgfa, &
                                                            nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: glb_cube, gsphere_bound, gub_cube
      LOGICAL                                            :: atom_pair_changed, use_subpatch
      REAL(KIND=dp)                                      :: adh, drmin, max_rad_ga, rab2, radius, &
                                                            rscale, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dh, dh_inv
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, pab, pin, sphi_a, sphi_b, work, &
                                                            zeta, zetb
      TYPE(task_type), DIMENSION(:), POINTER             :: tasks

      ! map all tasks on the grids
      tasks => task_list%tasks
      ntasks = task_list%ntasks

      nimages = task_list%nimages
      natoms = task_list%natoms
      maxpgf = task_list%maxpgf
      maxset = task_list%maxset
      maxco = task_list%maxco
      maxsgf_set = task_list%maxsgf_set

      ALLOCATE (pab(maxco, maxco), work(maxco, maxsgf_set))
      ALLOCATE (pin(1, 1))

      iatom_old = -1; jatom_old = -1; iset_old = -1; jset_old = -1; img_old = -1

      ! Loop over each gridlevel first, then loop and load balance over atom pairs
      loop_gridlevels: DO igrid_level = 1, task_list%ngrid_levels

         dh(1:3, 1:3) = task_list%cinfo(igrid_level)%dh(1:3, 1:3)
         adh = MAXVAL(ABS(dh))
         dh_inv(1:3, 1:3) = task_list%cinfo(igrid_level)%dh_inv(1:3, 1:3)

         perd(1:3) = task_list%cinfo(igrid_level)%perd(1:3)
         npts(1:3) = task_list%cinfo(igrid_level)%npts(1:3)
         lb_grid(1:3) = task_list%cinfo(igrid_level)%lb_grid(1:3)

         drmin = task_list%cinfo(igrid_level)%drmin
         max_radius = task_list%cinfo(igrid_level)%max_radius
         max_rad_ga = task_list%cinfo(igrid_level)%max_rad_ga

         NULLIFY (glb_cube, gub_cube, gsphere_bound)
         glb_cube => task_list%cinfo(igrid_level)%lb_cube
         gub_cube => task_list%cinfo(igrid_level)%ub_cube
         gsphere_bound => task_list%cinfo(igrid_level)%sphere_bounds

         loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
         loop_tasks: DO itask = task_list%taskstart(ipair, igrid_level), task_list%taskstop(ipair, igrid_level)
            img = tasks(itask)%image
            iatom = tasks(itask)%iatom
            jatom = tasks(itask)%jatom
            iset = tasks(itask)%iset
            jset = tasks(itask)%jset
            ipgf = tasks(itask)%ipgf
            jpgf = tasks(itask)%jpgf
            rab = tasks(itask)%rab
            radius = tasks(itask)%radius

            ikind = kind_of(iatom)
            jkind = kind_of(jatom)

            first_sgfa => task_list%basis(ikind)%first_sgf
            la_max => task_list%basis(ikind)%lmax
            la_min => task_list%basis(ikind)%lmin
            npgfa => task_list%basis(ikind)%npgf
            nsgfa => task_list%basis(ikind)%nsgf
            sphi_a => task_list%basis(ikind)%sphi
            zeta => task_list%basis(ikind)%zet

            first_sgfb => task_list%basis(jkind)%first_sgf
            lb_max => task_list%basis(jkind)%lmax
            lb_min => task_list%basis(jkind)%lmin
            npgfb => task_list%basis(jkind)%npgf
            nsgfb => task_list%basis(jkind)%nsgf
            sphi_b => task_list%basis(jkind)%sphi
            zetb => task_list%basis(jkind)%zet

            IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old .OR. img .NE. img_old) THEN
               IF (iatom .NE. iatom_old) ra(:) = posat(:, iatom)
               IF (iatom <= jatom) THEN
                  brow = iatom
                  bcol = jatom
               ELSE
                  brow = jatom
                  bcol = iatom
               END IF
               CALL spmat_get(spmats(img)%spmat, brow, bcol, p_block)
               DEALLOCATE (pin)
               n1 = SIZE(p_block, 1)
               n2 = SIZE(p_block, 2)
               IF (iatom <= jatom) THEN
                  ALLOCATE (pin(n1, n2))
                  pin(1:n1, 1:n2) = p_block(1:n1, 1:n2)
               ELSE
                  ALLOCATE (pin(n2, n1))
                  pin(1:n2, 1:n1) = TRANSPOSE(p_block(1:n1, 1:n2))
               END IF
               iatom_old = iatom
               jatom_old = jatom
               img_old = img
               atom_pair_changed = .TRUE.
            ELSE
               atom_pair_changed = .FALSE.
            ENDIF

            IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(iset)
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(jset)

               CALL dgemm("N", "N", ncoa, nsgfb(jset), nsgfa(iset), &
                          1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                          pin(sgfa, sgfb), SIZE(pin, 1), 0.0_dp, work(1, 1), maxco)
               CALL dgemm("N", "T", ncoa, ncob, nsgfb(jset), &
                          1.0_dp, work(1, 1), maxco, &
                          sphi_b(1, sgfb), SIZE(sphi_b, 1), 0.0_dp, pab(1, 1), maxco)

               iset_old = iset
               jset_old = jset
            ENDIF

            rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
            rb(:) = ra(:)+rab(:)
            zetp = zeta(ipgf, iset)+zetb(jpgf, jset)

            IF (2.0_dp*radius < adh) CYCLE

            na1 = (ipgf-1)*ncoset(la_max(iset))+1
            nb1 = (jpgf-1)*ncoset(lb_max(jset))+1

            ! takes the density matrix symmetry in account
            IF (iatom == jatom) THEN
               rscale = 1.0_dp
            ELSE
               rscale = 2.0_dp
            END IF

            ! check whether we need to use the generalised collocation scheme
            IF (tasks(itask)%dist_type == 2) THEN
               use_subpatch = .TRUE.
            ELSE
               use_subpatch = .FALSE.
            ENDIF

            CALL collocate_pgf_product_rspace( &
               la_max(iset), zeta(ipgf, iset), la_min(iset), &
               lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
               ra, rab, rab2, radius, rscale, pab, na1-1, nb1-1, &
               rs_rho(igrid_level)%grid, &
               dh, dh_inv, perd, npts, lb_grid, &
               drmin, max_radius, max_rad_ga, &
               glb_cube, gub_cube, gsphere_bound, &
               compute_tau, use_subpatch, lmax_global, maxco)

         END DO loop_tasks
         END DO loop_pairs

      END DO loop_gridlevels

      DEALLOCATE (pab, pin, work)

   END SUBROUTINE grid_collocate_ref

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
!> \param rab2 ...
!> \param radius ...
!> \param rscale ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param grid ...
!> \param dh ...
!> \param dh_inv ...
!> \param perd ...
!> \param npts ...
!> \param lb_grid ...
!> \param drmin ...
!> \param max_radius ...
!> \param max_rad_ga ...
!> \param glb_cube ...
!> \param gub_cube ...
!> \param gsphere_bound ...
!> \param compute_tau ...
!> \param use_subpatch ...
!> \param lmax_global Maximum possible value of lmax used to dimension arrays
!> \param maxco ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_rspace(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, rab2, radius, rscale, pab, o1, o2, &
                                           grid, &
                                           dh, dh_inv, perd, npts, lb_grid, &
                                           drmin, max_radius, max_rad_ga, &
                                           glb_cube, gub_cube, gsphere_bound, &
                                           compute_tau, use_subpatch, lmax_global, maxco)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: ra, rab
      REAL(KIND=dp), INTENT(IN)                          :: rab2, radius, rscale
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: dh, dh_inv
      INTEGER, DIMENSION(3), INTENT(IN)                  :: perd, npts, lb_grid
      REAL(KIND=dp), INTENT(IN)                          :: drmin
      INTEGER, INTENT(IN)                                :: max_radius
      REAL(KIND=dp), INTENT(IN)                          :: max_rad_ga
      INTEGER, DIMENSION(:, :), POINTER                  :: glb_cube, gub_cube, gsphere_bound
      LOGICAL, INTENT(IN)                                :: compute_tau, use_subpatch
      INTEGER, INTENT(IN)                                :: lmax_global, maxco

      INTEGER                                            :: gridbounds(2, 3), la_max_local, &
                                                            la_min_local, lb_max_local, &
                                                            lb_min_local, lp
      INTEGER, DIMENSION(3)                              :: ng
      REAL(KIND=dp)                                      :: f, prefactor, zetp
      REAL(kind=dp), DIMENSION(((lmax_global*2+1)*(&
         lmax_global*2+2)*(lmax_global*2+3))/6)          :: coef_xyz
      REAL(kind=dp), DIMENSION(0:lmax_global*2, 0:&
         lmax_global, 0:lmax_global, 3)                  :: alpha
      REAL(KIND=dp), DIMENSION(3)                        :: rb, rp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab_local

      zetp = zeta+zetb
      f = zetb/zetp
      rp(:) = ra(:)+f*rab(:)
      rb(:) = ra(:)+rab(:)
      prefactor = rscale*EXP(-zeta*f*rab2)

      IF (compute_tau) THEN
         la_max_local = la_max+1
         la_min_local = MAX(la_min-1, 0)
         lb_max_local = lb_max+1
         lb_min_local = MAX(lb_min-1, 0)
         NULLIFY (pab_local)
         CALL prepare_pab_tau(pab_local, pab, o1, o2, &
                              la_max, la_min, lb_max, lb_min, maxco, zeta, zetb)
      ELSE
         la_max_local = la_max
         la_min_local = la_min
         lb_max_local = lb_max
         lb_min_local = lb_min
         NULLIFY (pab_local)
         CALL prepare_pab_rho(pab_local, pab, o1, o2, la_max, la_min, lb_max, lb_min, maxco)
      END IF

      ng(:) = npts(:)
      gridbounds(1, 1) = LBOUND(GRID, 1)
      gridbounds(2, 1) = UBOUND(GRID, 1)
      gridbounds(1, 2) = LBOUND(GRID, 2)
      gridbounds(2, 2) = UBOUND(GRID, 2)
      gridbounds(1, 3) = LBOUND(GRID, 3)
      gridbounds(2, 3) = UBOUND(GRID, 3)

!   *** initialise the coefficient matrix, we transform the sum
!
!   sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} *
!           (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
!
!   into
!
!   sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
!
!   where p is center of the product gaussian, and lp = la_max + lb_max
!   (current implementation is l**7)
!
      CALL prepare_alpha(alpha, ra, rb, rp, la_max_local, lb_max_local, lmax_global)
!
!   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
!   use a three step procedure
!   we don't store zeros, so counting is done using lxyz,lxy in order to have
!   contiguous memory access in collocate_fast.F
!
      CALL prepare_coef(coef_xyz, alpha, pab_local, prefactor, &
                        la_max_local, la_min_local, lb_max_local, lb_min_local, lmax_global)

      lp = la_max_local+lb_max_local

      IF (use_subpatch) THEN
         CALL collocate_general(grid, dh, dh_inv, lp, lmax_global, ng, lb_grid, perd, &
                                rp, zetp, radius, max_rad_ga, coef_xyz)
      ELSE
         CALL collocate_ortho(grid, dh, dh_inv, ng, perd, glb_cube, gub_cube, gsphere_bound, gridbounds, lb_grid, &
                              lp, zetp, rp, radius, max_radius, drmin, coef_xyz)
      END IF

      DEALLOCATE (pab_local)

   END SUBROUTINE collocate_pgf_product_rspace

! **************************************************************************************************
!> \brief ...
!> \param pab_tau ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param la_max ...
!> \param la_min ...
!> \param lb_max ...
!> \param lb_min ...
!> \param maxco ...
!> \param zeta ...
!> \param zetb ...
! **************************************************************************************************
   SUBROUTINE prepare_pab_tau(pab_tau, pab, o1, o2, &
                              la_max, la_min, lb_max, lb_min, maxco, zeta, zetb)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab_tau
      REAL(dp), DIMENSION(:, :), INTENT(IN), TARGET      :: pab
      INTEGER, INTENT(in)                                :: o1, o2, la_max, la_min, lb_max, lb_min, &
                                                            maxco
      REAL(KIND=dp), INTENT(in)                          :: zeta, zetb

      INTEGER                                            :: istat, nla, nlb
      INTERFACE
         FUNCTION grid_prepare_pab_tau(o1, o2, la_max, la_min, lb_max, lb_min, maxco, &
                                       zeta, zetb, pab, pab_tau) &
            RESULT(istat) BIND(C, name="grid_prepare_pab_tau")
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            INTEGER(KIND=C_INT), VALUE                :: maxco
            REAL(KIND=C_DOUBLE), VALUE                :: zeta
            REAL(KIND=C_DOUBLE), VALUE                :: zetb
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: pab_tau
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_prepare_pab_tau
      END INTERFACE

      nla = ncoset(la_max+1)
      nlb = ncoset(lb_max+1)
      ALLOCATE (pab_tau(nla, nlb))

      istat = grid_prepare_pab_tau(o1, &
                                   o2, &
                                   la_max, &
                                   la_min, &
                                   lb_max, &
                                   lb_min, &
                                   maxco, &
                                   zeta, &
                                   zetb, &
                                   C_LOC(pab(1, 1)), &
                                   C_LOC(pab_tau(1, 1)))
      CPASSERT(istat == 0)

   END SUBROUTINE prepare_pab_tau

! **************************************************************************************************
!> \brief ...
!> \param pab_rho ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param la_max ...
!> \param la_min ...
!> \param lb_max ...
!> \param lb_min ...
!> \param maxco ...
! **************************************************************************************************
   SUBROUTINE prepare_pab_rho(pab_rho, pab, o1, o2, la_max, la_min, lb_max, lb_min, maxco)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab_rho
      REAL(dp), DIMENSION(:, :), INTENT(IN), TARGET      :: pab
      INTEGER, INTENT(in)                                :: o1, o2, la_max, la_min, lb_max, lb_min, &
                                                            maxco

      INTEGER                                            :: istat, nla, nlb
      INTERFACE
         FUNCTION grid_prepare_pab_rho(o1, o2, la_max, la_min, lb_max, lb_min, maxco, pab, pab_rho) &
            RESULT(istat) BIND(C, name="grid_prepare_pab_rho")
            IMPORT :: C_PTR, C_INT
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            INTEGER(KIND=C_INT), VALUE                :: maxco
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: pab_rho
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_prepare_pab_rho
      END INTERFACE

      nla = ncoset(la_max)
      nlb = ncoset(lb_max)
      ALLOCATE (pab_rho(nla, nlb))

      istat = grid_prepare_pab_rho(o1, &
                                   o2, &
                                   la_max, &
                                   la_min, &
                                   lb_max, &
                                   lb_min, &
                                   maxco, &
                                   C_LOC(pab(1, 1)), &
                                   C_LOC(pab_rho(1, 1)))
      CPASSERT(istat == 0)

   END SUBROUTINE prepare_pab_rho

! **************************************************************************************************
!> \brief ...
!> \param alpha ...
!> \param ra ...
!> \param rb ...
!> \param rp ...
!> \param la_max ...
!> \param lb_max ...
!> \param lmax ...
! **************************************************************************************************
   SUBROUTINE prepare_alpha(alpha, ra, rb, rp, la_max, lb_max, lmax)
      REAL(KIND=dp), DIMENSION(:, :, :, :), &
         INTENT(OUT), TARGET                             :: alpha
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rb, rp
      INTEGER, INTENT(IN)                                :: la_max, lb_max, lmax

      INTEGER                                            :: istat
      INTERFACE
         FUNCTION grid_prepare_alpha(ra, rb, rp, la_max, lb_max, lmax, alpha) &
            RESULT(istat) BIND(C, name="grid_prepare_alpha")
            IMPORT :: C_PTR, C_INT
            TYPE(C_PTR), VALUE                        :: ra
            TYPE(C_PTR), VALUE                        :: rb
            TYPE(C_PTR), VALUE                        :: rp
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lmax
            TYPE(C_PTR), VALUE                        :: alpha
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_prepare_alpha
      END INTERFACE

      istat = grid_prepare_alpha(C_LOC(ra(1)), &
                                 C_LOC(rb(1)), &
                                 C_LOC(rp(1)), &
                                 la_max, &
                                 lb_max, &
                                 lmax, &
                                 C_LOC(alpha(1, 1, 1, 1)))
      CPASSERT(istat == 0)

   END SUBROUTINE prepare_alpha

! **************************************************************************************************
!> \brief ...
!> \param coef_xyz ...
!> \param alpha ...
!> \param pab ...
!> \param prefactor ...
!> \param la_max ...
!> \param la_min ...
!> \param lb_max ...
!> \param lb_min ...
!> \param lmax ...
! **************************************************************************************************
   SUBROUTINE prepare_coef(coef_xyz, alpha, pab, prefactor, la_max, la_min, lb_max, lb_min, lmax)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT), TARGET :: coef_xyz
      REAL(KIND=dp), DIMENSION(:, :, :, :), INTENT(IN), &
         TARGET                                          :: alpha
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), TARGET :: pab
      REAL(KIND=dp), INTENT(IN)                          :: prefactor
      INTEGER, INTENT(IN)                                :: la_max, la_min, lb_max, lb_min, lmax

      INTEGER                                            :: istat
      INTERFACE
         FUNCTION grid_prepare_coef(la_max, la_min, lb_max, lb_min, lmax, prefactor, &
                                    alpha, pab, coef_xyz) &
            RESULT(istat) BIND(C, name="grid_prepare_coef")
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            INTEGER(KIND=C_INT), VALUE                :: lmax
            REAL(KIND=C_DOUBLE), VALUE                :: prefactor
            TYPE(C_PTR), VALUE                        :: alpha
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: coef_xyz
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_prepare_coef
      END INTERFACE

      istat = grid_prepare_coef(la_max, la_min, &
                                lb_max, lb_min, &
                                lmax, &
                                prefactor, &
                                C_LOC(alpha(1, 1, 1, 1)), &
                                C_LOC(pab(1, 1)), &
                                C_LOC(coef_xyz(1)))
      CPASSERT(istat == 0)

   END SUBROUTINE prepare_coef

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param dh ...
!> \param dh_inv ...
!> \param ng ...
!> \param perd ...
!> \param glb_cube ...
!> \param gub_cube ...
!> \param gsphere_bound ...
!> \param gridbounds ...
!> \param lb_grid ...
!> \param lp ...
!> \param zetp ...
!> \param rp ...
!> \param radius ...
!> \param max_radius ...
!> \param drmin ...
!> \param coef_xyz ...
! **************************************************************************************************
   SUBROUTINE collocate_ortho(grid, dh, dh_inv, ng, perd, glb_cube, gub_cube, gsphere_bound, gridbounds, lb_grid, &
                              lp, zetp, rp, radius, max_radius, drmin, coef_xyz)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: dh, dh_inv
      INTEGER, DIMENSION(3), INTENT(IN)                  :: ng, perd
      INTEGER, DIMENSION(:, :), POINTER                  :: glb_cube, gub_cube, gsphere_bound
      INTEGER, DIMENSION(2, 3), INTENT(IN)               :: gridbounds
      INTEGER, DIMENSION(3), INTENT(IN)                  :: lb_grid
      INTEGER, INTENT(IN)                                :: lp
      REAL(KIND=dp), INTENT(IN)                          :: zetp
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: rp
      REAL(KIND=dp), INTENT(IN)                          :: radius
      INTEGER, INTENT(IN)                                :: max_radius
      REAL(KIND=dp), INTENT(IN)                          :: drmin
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: coef_xyz

      INTEGER                                            :: cmax, i, imr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: map
      INTEGER, DIMENSION(3)                              :: cubecenter, lb_cube, ub_cube
      INTEGER, DIMENSION(:), POINTER                     :: sphere_bounds
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: pol
      REAL(KIND=dp), DIMENSION(3)                        :: dr, roffset

      dr(1) = dh(1, 1)
      dr(2) = dh(2, 2)
      dr(3) = dh(3, 3)

      imr = MAX(1, CEILING(radius/drmin))
      IF (imr .GT. max_radius) THEN
         !
         ! This is an important check. If the required radius for mapping the density is different
         ! from the actual computed one, (significant) errors can occur.
         ! This error can invariably be fixed by improving the computation of maxradius
         ! in the call to init_cube_info
         !
         STOP
      ENDIF
      lb_cube(:) = glb_cube(:, imr)
      ub_cube(:) = gub_cube(:, imr)
      sphere_bounds(1:) => gsphere_bound(:, imr)

      cmax = MAXVAL(ub_cube)

!   *** position of the gaussian product
!
!   this is the actual definition of the position on the grid
!   i.e. a point rp(:) gets here grid coordinates
!   MODULO(rp(:)/dr(:),ng(:))+1
!   hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)
!

      cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
      roffset(:) = rp(:)-REAL(cubecenter(:), dp)*dr(:)

      !   *** a mapping so that the ig corresponds to the right grid point
      ALLOCATE (map(-cmax:cmax, 3))
      DO i = 1, 3
         CALL fill_map(perd(i) == 1, lb_cube(i), ub_cube(i), cubecenter(i), &
                       lb_grid(i), LBOUND(grid, i), UBOUND(grid, i), ng(i), cmax, map(:, i))
      ENDDO

      ALLOCATE (pol(1:2, 0:lp, -cmax:0, 1:3))
      DO i = 1, 3
         CALL fill_pol(dr(i), roffset(i), lb_cube(i), lp, cmax, zetp, pol(:, :, :, i))
      END DO

      CALL collocate_core(grid, coef_xyz, pol, map, sphere_bounds, lp, cmax, gridbounds)

      ! deallocation needed to pass around a pgi bug..
      DEALLOCATE (pol)
      DEALLOCATE (map)

   END SUBROUTINE collocate_ortho

! **************************************************************************************************
!> \brief Fill mapping so that the ig corresponds to the right grid point
!> \param periodic ...
!> \param lb_cube ...
!> \param ub_cube ...
!> \param cubecenter ...
!> \param lb_grid ...
!> \param grid_lbound ...
!> \param grid_ubound ...
!> \param ng ...
!> \param cmax ...
!> \param map ...
! **************************************************************************************************
   SUBROUTINE fill_map(periodic, lb_cube, ub_cube, cubecenter, lb_grid, grid_lbound, grid_ubound, ng, cmax, map)
      LOGICAL, INTENT(IN)                                :: periodic
      INTEGER, INTENT(IN)                                :: lb_cube, ub_cube, cubecenter, lb_grid, &
                                                            grid_lbound, grid_ubound, ng, cmax
      INTEGER, DIMENSION(:), INTENT(OUT), TARGET         :: map

      INTEGER                                            :: istat
      INTERFACE
         FUNCTION grid_fill_map(periodic, lb_cube, ub_cube, cubecenter, lb_grid, &
                                grid_lbound, grid_ubound, ng, cmax, map) &
            RESULT(istat) BIND(C, name="grid_fill_map")
            IMPORT :: C_PTR, C_INT, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: periodic
            INTEGER(KIND=C_INT), VALUE                :: lb_cube
            INTEGER(KIND=C_INT), VALUE                :: ub_cube
            INTEGER(KIND=C_INT), VALUE                :: cubecenter
            INTEGER(KIND=C_INT), VALUE                :: lb_grid
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound
            INTEGER(KIND=C_INT), VALUE                :: grid_ubound
            INTEGER(KIND=C_INT), VALUE                :: ng
            INTEGER(KIND=C_INT), VALUE                :: cmax
            TYPE(C_PTR), VALUE                        :: map
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_fill_map
      END INTERFACE

      istat = grid_fill_map(LOGICAL(periodic, KIND=C_BOOL), lb_cube, ub_cube, cubecenter, lb_grid, &
                            grid_lbound, grid_ubound, ng, cmax, C_LOC(map(1)))
      CPASSERT(istat == 0)

   END SUBROUTINE fill_map

! **************************************************************************************************
!> \brief ...
!> \param dr ...
!> \param roffset ...
!> \param lb_cube ...
!> \param lp ...
!> \param cmax ...
!> \param zetp ...
!> \param pol ...
! **************************************************************************************************
   SUBROUTINE fill_pol(dr, roffset, lb_cube, lp, cmax, zetp, pol)
      REAL(KIND=dp), INTENT(IN)                          :: dr, roffset
      INTEGER, INTENT(IN)                                :: lb_cube, lp, cmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp
      REAL(kind=dp), DIMENSION(:, :, :), INTENT(OUT), &
         TARGET                                          :: pol

      INTEGER                                            :: istat
      INTERFACE
         FUNCTION grid_fill_pol(dr, roffset, lb_cube, lp, cmax, zetp, pol) &
            RESULT(istat) BIND(C, name="grid_fill_pol")
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            REAL(KIND=C_DOUBLE), VALUE                :: dr
            REAL(KIND=C_DOUBLE), VALUE                :: roffset
            INTEGER(KIND=C_INT), VALUE                :: lb_cube
            INTEGER(KIND=C_INT), VALUE                :: lp
            INTEGER(KIND=C_INT), VALUE                :: cmax
            REAL(KIND=C_DOUBLE), VALUE                :: zetp
            TYPE(C_PTR), VALUE                        :: pol
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_fill_pol
      END INTERFACE

      istat = grid_fill_pol(dr, roffset, lb_cube, lp, cmax, zetp, C_LOC(pol(1, 1, 1)))
      CPASSERT(istat == 0)

   END SUBROUTINE fill_pol

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param coef_xyz ...
!> \param pol ...
!> \param map ...
!> \param sphere_bounds ...
!> \param lp ...
!> \param cmax ...
!> \param gridbounds ...
! **************************************************************************************************
   SUBROUTINE collocate_core(grid, coef_xyz, pol, map, sphere_bounds, lp, cmax, gridbounds)
      REAL(dp), INTENT(IN), TARGET                       :: pol(:, :, :, :)
      INTEGER, INTENT(IN), TARGET                        :: map(:, :), sphere_bounds(:), lp
      REAL(dp), INTENT(IN), TARGET :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)
      INTEGER, INTENT(IN)                                :: cmax
      INTEGER, INTENT(IN), TARGET                        :: gridbounds(2, 3)
      REAL(dp), INTENT(INOUT), TARGET :: grid(gridbounds(1, 1):gridbounds(2, 1), gridbounds(1, 2): &
         gridbounds(2, 2), gridbounds(1, 3):gridbounds(2, 3))

      INTEGER                                            :: istat
      INTERFACE
         FUNCTION grid_collocate_core(grid_size_x, grid_size_y, grid_size_z, &
                                      grid_lbound_x, grid_lbound_y, grid_lbound_z, &
                                      lp, cmax, &
                                      coef_xyz, pol, map, sphere_bounds, grid) &
            RESULT(istat) BIND(C, name="grid_collocate_core")
            IMPORT :: C_PTR, C_INT
            INTEGER(KIND=C_INT), VALUE                :: grid_size_x
            INTEGER(KIND=C_INT), VALUE                :: grid_size_y
            INTEGER(KIND=C_INT), VALUE                :: grid_size_z
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound_x
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound_y
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound_z
            INTEGER(KIND=C_INT), VALUE                :: lp
            INTEGER(KIND=C_INT), VALUE                :: cmax
            TYPE(C_PTR), VALUE                        :: coef_xyz
            TYPE(C_PTR), VALUE                        :: pol
            TYPE(C_PTR), VALUE                        :: map
            TYPE(C_PTR), VALUE                        :: sphere_bounds
            TYPE(C_PTR), VALUE                        :: grid
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_collocate_core
      END INTERFACE

      ! TODO: Are LBOUND(grid, :) always 1 ???
      istat = grid_collocate_core(grid_size_x=SIZE(grid, 1), &
                                  grid_size_y=SIZE(grid, 2), &
                                  grid_size_z=SIZE(grid, 3), &
                                  grid_lbound_x=LBOUND(grid, 1), &
                                  grid_lbound_y=LBOUND(grid, 2), &
                                  grid_lbound_z=LBOUND(grid, 3), &
                                  lp=lp, &
                                  cmax=cmax, &
                                  coef_xyz=C_LOC(coef_xyz(1)), &
                                  pol=C_LOC(pol(1, 1, 1, 1)), &
                                  map=C_LOC(map(1, 1)), &
                                  sphere_bounds=C_LOC(sphere_bounds(1)), &
                                  grid=C_LOC(grid(LBOUND(grid, 1), LBOUND(grid, 2), LBOUND(grid, 3))))

      CPASSERT(istat == 0)

   END SUBROUTINE collocate_core

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param dh ...
!> \param dh_inv ...
!> \param lp ...
!> \param lmax ...
!> \param ng ...
!> \param lb_grid ...
!> \param perd ...
!> \param rp ...
!> \param zetp ...
!> \param radius ...
!> \param max_rad_ga ...
!> \param coef_xyz ...
! **************************************************************************************************
   SUBROUTINE collocate_general(grid, dh, dh_inv, lp, lmax, ng, lb_grid, perd, rp, zetp, radius, &
                                max_rad_ga, coef_xyz)
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), TARGET :: dh, dh_inv
      INTEGER, INTENT(in)                                :: lp, lmax
      INTEGER, DIMENSION(3), INTENT(in), TARGET          :: ng, lb_grid, perd
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: rp
      REAL(KIND=dp), INTENT(IN)                          :: zetp, radius, max_rad_ga
      REAL(KIND=dp), DIMENSION(:), INTENT(in), TARGET    :: coef_xyz

      INTEGER                                            :: istat
      INTERFACE
         FUNCTION grid_collocate_general(grid_size_x, grid_size_y, grid_size_z, &
                                         grid_lbound_x, grid_lbound_y, grid_lbound_z, &
                                         lp, lmax, radius, zetp, coef_xyz, dh, dh_inv, &
                                         rp, ng, lb_grid, perd, grid) &
            RESULT(istat) BIND(C, name="grid_collocate_general")
            IMPORT :: C_PTR, C_INT, C_DOUBLE
            INTEGER(KIND=C_INT), VALUE                :: grid_size_x
            INTEGER(KIND=C_INT), VALUE                :: grid_size_y
            INTEGER(KIND=C_INT), VALUE                :: grid_size_z
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound_x
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound_y
            INTEGER(KIND=C_INT), VALUE                :: grid_lbound_z
            INTEGER(KIND=C_INT), VALUE                :: lp
            INTEGER(KIND=C_INT), VALUE                :: lmax
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            REAL(KIND=C_DOUBLE), VALUE                :: zetp
            TYPE(C_PTR), VALUE                        :: coef_xyz
            TYPE(C_PTR), VALUE                        :: dh
            TYPE(C_PTR), VALUE                        :: dh_inv
            TYPE(C_PTR), VALUE                        :: rp
            TYPE(C_PTR), VALUE                        :: ng
            TYPE(C_PTR), VALUE                        :: lb_grid
            TYPE(C_PTR), VALUE                        :: perd
            TYPE(C_PTR), VALUE                        :: grid
            INTEGER(KIND=C_INT)                       :: istat
         END FUNCTION grid_collocate_general
      END INTERFACE

      IF (radius > max_rad_ga) THEN
         !
         ! This is an important check. If the required radius for mapping the density is different
         ! from the actual computed one, (significant) errors can occur.
         ! This error can invariably be fixed by improving the computation of maxradius
         ! in the call to init_cube_info
         !
         STOP
      ENDIF

      ! TODO: Are LBOUND(grid, :) always 1 ???
      istat = grid_collocate_general(grid_size_x=SIZE(grid, 1), &
                                     grid_size_y=SIZE(grid, 2), &
                                     grid_size_z=SIZE(grid, 3), &
                                     grid_lbound_x=LBOUND(grid, 1), &
                                     grid_lbound_y=LBOUND(grid, 2), &
                                     grid_lbound_z=LBOUND(grid, 3), &
                                     lp=lp, &
                                     lmax=lmax, &
                                     radius=radius, &
                                     zetp=zetp, &
                                     coef_xyz=C_LOC(coef_xyz(1)), &
                                     dh=C_LOC(dh(1, 1)), &
                                     dh_inv=C_LOC(dh_inv(1, 1)), &
                                     rp=C_LOC(rp(1)), &
                                     ng=C_LOC(ng(1)), &
                                     lb_grid=C_LOC(lb_grid(1)), &
                                     perd=C_LOC(perd(1)), &
                                     grid=C_LOC(grid(LBOUND(grid, 1), LBOUND(grid, 2), LBOUND(grid, 3))))

      CPASSERT(istat == 0)

   END SUBROUTINE collocate_general

END MODULE grid_base_ref

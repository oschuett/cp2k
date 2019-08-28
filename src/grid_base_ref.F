!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE grid_base_ref
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE task_list_methods,               ONLY: int2pair
   USE task_list_types,                 ONLY: task_list_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_base_ref'

   TYPE grid_base_type
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER       :: grid => NULL()
   END TYPE grid_base_type

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

   PUBLIC :: grid_collocate_ref
   PUBLIC :: grid_base_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param rs_rho Input/Output
!> \param deltap input
!> \param task_list input
!> \param compute_tau ...
!> \param kind_of input
!> \param posat input
!> \param lmax_global input
! **************************************************************************************************
   SUBROUTINE grid_collocate_ref(rs_rho, deltap, task_list, compute_tau, kind_of, posat, lmax_global)

      TYPE(grid_base_type), DIMENSION(:)                 :: rs_rho
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: deltap
      TYPE(task_list_type), POINTER                      :: task_list
      LOGICAL, INTENT(IN)                                :: compute_tau
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      INTEGER, INTENT(IN)                                :: lmax_global

      CHARACTER(len=*), PARAMETER :: routineN = 'grid_collocate_ref', &
         routineP = moduleN//':'//routineN

      INTEGER :: bcol, brow, iatom, iatom_old, igrid_level, igrid_level_dummy, ikind, img, &
         img_old, ipair, ipgf, iset, iset_old, itask, jatom, jatom_old, jkind, jpgf, jset, &
         jset_old, max_radius, maxco, maxpgf, maxset, maxsgf_set, na1, na2, natoms, nb1, nb2, &
         ncoa, ncob, nimages, ntasks, sgfa, sgfb
      INTEGER(kind=int_8), DIMENSION(:, :), POINTER      :: tasks
      INTEGER, DIMENSION(3)                              :: lb_grid(3), npts(3), perd(3)
      INTEGER, DIMENSION(:), POINTER                     :: first_sgfa, first_sgfb, la_max, la_min, &
                                                            lb_max, lb_min, npgfa, npgfb, nsgfa, &
                                                            nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: glb_cube, gsphere_bound, gub_cube
      LOGICAL                                            :: atom_pair_changed, found, use_subpatch
      REAL(KIND=dp)                                      :: drmin, max_rad_ga, rab2, radius, rscale, &
                                                            zetp
      REAL(KIND=dp), DIMENSION(3)                        :: ra, rab, rab_inv, rb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dh, dh_inv
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dist_ab, p_block, pab, sphi_a, sphi_b, &
                                                            work, zeta, zetb

      ! map all tasks on the grids
      tasks => task_list%tasks
      dist_ab => task_list%dist_ab
      ntasks = task_list%ntasks

      nimages = task_list%nimages
      natoms = task_list%natoms
      maxpgf = task_list%maxpgf
      maxset = task_list%maxset
      maxco = task_list%maxco
      maxsgf_set = task_list%maxsgf_set

      ALLOCATE (pab(maxco, maxco), work(maxco, maxsgf_set))

      iatom_old = -1; jatom_old = -1; iset_old = -1; jset_old = -1; img_old = -1

      ! Loop over each gridlevel first, then loop and load balance over atom pairs
      loop_gridlevels: DO igrid_level = 1, task_list%ngrid_levels

         dh(1:3, 1:3) = task_list%cinfo(igrid_level)%dh(1:3, 1:3)
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
            !decode the atom pair and basis info (igrid_level_dummy equals do loop variable by construction).
            CALL int2pair(tasks(3, itask), igrid_level_dummy, img, iatom, jatom, iset, jset, ipgf, jpgf, &
                          nimages, natoms, maxset, maxpgf)
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

               CALL dbcsr_get_block_p(matrix=deltap(img)%matrix, &
                                      row=brow, col=bcol, BLOCK=p_block, found=found)
               CPASSERT(found)

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

               IF (iatom <= jatom) THEN
                  CALL dgemm("N", "N", ncoa, nsgfb(jset), nsgfa(iset), &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             p_block(sgfa, sgfb), SIZE(p_block, 1), &
                             0.0_dp, work(1, 1), maxco)
                  CALL dgemm("N", "T", ncoa, ncob, nsgfb(jset), &
                             1.0_dp, work(1, 1), maxco, &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, pab(1, 1), maxco)
               ELSE
                  CALL dgemm("N", "N", ncob, nsgfa(iset), nsgfb(jset), &
                             1.0_dp, sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             p_block(sgfb, sgfa), SIZE(p_block, 1), &
                             0.0_dp, work(1, 1), maxco)
                  CALL dgemm("N", "T", ncob, ncoa, nsgfa(iset), &
                             1.0_dp, work(1, 1), maxco, &
                             sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             0.0_dp, pab(1, 1), maxco)
               END IF

               iset_old = iset
               jset_old = jset
            ENDIF

            rab(1:3) = dist_ab(1:3, itask)
            rab2 = rab(1)*rab(1)+rab(2)*rab(2)+rab(3)*rab(3)
            rb(:) = ra(:)+rab(:)
            zetp = zeta(ipgf, iset)+zetb(jpgf, jset)

            radius = dist_ab(4, itask)

            na1 = (ipgf-1)*ncoset(la_max(iset))+1
            na2 = ipgf*ncoset(la_max(iset))
            nb1 = (jpgf-1)*ncoset(lb_max(jset))+1
            nb2 = jpgf*ncoset(lb_max(jset))

            ! takes the density matrix symmetry in account
            IF (iatom == jatom) THEN
               rscale = 1.0_dp
            ELSE
               rscale = 2.0_dp
            END IF

            ! check whether we need to use the generalised collocation scheme
            IF (tasks(4, itask) == 2) THEN
               use_subpatch = .TRUE.
            ELSE
               use_subpatch = .FALSE.
            ENDIF

            IF (iatom <= jatom) THEN
               CALL collocate_pgf_product_rspace( &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  ra, rab, rab2, radius, rscale, pab, na1-1, nb1-1, &
                  rs_rho(igrid_level)%grid, &
                  dh, dh_inv, &
                  perd, npts, lb_grid, &
                  drmin, max_radius, max_rad_ga, &
                  glb_cube, gub_cube, gsphere_bound, &
                  compute_tau, &
                  use_subpatch, lmax_global)
            ELSE
               rab_inv = -rab
               CALL collocate_pgf_product_rspace( &
                  lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                  la_max(iset), zeta(ipgf, iset), la_min(iset), &
                  rb, rab_inv, rab2, radius, rscale, pab, nb1-1, na1-1, &
                  rs_rho(igrid_level)%grid, &
                  dh, dh_inv, &
                  perd, npts, lb_grid, &
                  drmin, max_radius, max_rad_ga, &
                  glb_cube, gub_cube, gsphere_bound, &
                  compute_tau, &
                  use_subpatch, lmax_global)
            END IF
         END DO loop_tasks
         END DO loop_pairs

      END DO loop_gridlevels

      DEALLOCATE (pab, work)

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
!> \param radius_task ...
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
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_rspace(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, rab2, radius_task, rscale, pab, o1, o2, &
                                           grid, &
                                           dh, dh_inv, &
                                           perd, npts, lb_grid, &
                                           drmin, max_radius, max_rad_ga, &
                                           glb_cube, gub_cube, gsphere_bound, &
                                           compute_tau, &
                                           use_subpatch, lmax_global)

      INTEGER, INTENT(IN)                      :: la_max
      REAL(KIND=dp), INTENT(IN)                :: zeta
      INTEGER, INTENT(IN)                      :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                :: zetb
      INTEGER, INTENT(IN)                      :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: ra, rab
      REAL(KIND=dp), INTENT(IN)                :: rab2, radius_task, rscale
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab
      INTEGER, INTENT(IN)                      :: o1, o2
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER  :: grid
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN):: dh, dh_inv
      INTEGER, DIMENSION(3), INTENT(IN)        :: perd, npts, lb_grid
      REAL(KIND=dp), INTENT(IN)                :: drmin
      INTEGER, INTENT(IN)                      :: max_radius
      REAL(KIND=dp), INTENT(IN)                :: max_rad_ga
      INTEGER, DIMENSION(:, :), POINTER         :: glb_cube
      INTEGER, DIMENSION(:, :), POINTER         :: gub_cube
      INTEGER, DIMENSION(:, :), POINTER         :: gsphere_bound
      LOGICAL, INTENT(IN)                      :: compute_tau
      LOGICAL, INTENT(IN)                      :: use_subpatch
      INTEGER, INTENT(IN)                      :: lmax_global

      CHARACTER(len=*), PARAMETER :: routineN = 'collocate_pgf_product_rspace', &
                                     routineP = moduleN//':'//routineN

      INTEGER :: cmax, gridbounds(2, 3), i, ico, icoef, ig, imr, &
                 jco, k, l, la_max_local, la_min_local, lb_max_local, lb_min_local, &
                 length, lxa, lxb, lxy, lxyz, lya, lyb, &
                 lza, lzb, offset, start
      INTEGER, DIMENSION(3)                    :: cubecenter, lb_cube, ng, &
                                                  ub_cube
      INTEGER, DIMENSION(:), POINTER           :: sphere_bounds
      REAL(KIND=dp) :: a, b, binomial_k_lxa, binomial_l_lxb, f, pg, &
                       prefactor, radius, rpg, zetp
      REAL(KIND=dp), DIMENSION(3)              :: dr, roffset, rp
      REAL(KIND=dp), DIMENSION(:, :), POINTER  :: pab_local

      INTEGER :: lxp, lyp, lzp, lp, lxpm, iaxis
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: map
      REAL(kind=dp) :: p_ele
      REAL(kind=dp), DIMENSION(0:lmax_global*2, 0:lmax_global, 0:lmax_global, 3) :: alpha
      REAL(kind=dp), DIMENSION(((lmax_global*2+1)*(lmax_global*2+2)*(lmax_global*2+3))/6) :: coef_xyz
      REAL(kind=dp), DIMENSION(((lmax_global*2+1)*(lmax_global*2+2))/2) :: coef_xyt
      REAL(kind=dp), DIMENSION(0:lmax_global*2) :: coef_xtt

      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :) :: pol_z
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :) :: pol_y
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :) :: pol_x
      REAL(KIND=dp) :: t_exp_1, t_exp_2, t_exp_min_1, t_exp_min_2, t_exp_plus_1, t_exp_plus_2

      zetp = zeta+zetb
      f = zetb/zetp
      rp(:) = ra(:)+f*rab(:)

      radius = radius_task
      prefactor = rscale*EXP(-zeta*f*rab2)

      a = MAXVAL(ABS(dh))
      IF (radius .LT. a/2.0_dp) THEN
         RETURN
      END IF

      IF (compute_tau) THEN
         la_max_local = la_max+1
         la_min_local = MAX(la_min-1, 0)
         lb_max_local = lb_max+1
         lb_min_local = MAX(lb_min-1, 0)
         NULLIFY (pab_local)
         CALL prepare_pab_tau(pab_local, pab, o1, o2, &
                              la_max, la_min, lb_max, lb_min, zeta, zetb)
      ELSE
         la_max_local = la_max
         la_min_local = la_min
         lb_max_local = lb_max
         lb_min_local = lb_min
         NULLIFY (pab_local)
         CALL prepare_pab_rho(pab_local, pab, o1, o2, la_max, la_min, lb_max, lb_min)
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
      lp = la_max_local+lb_max_local
!
!   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
!
!   *** make the alpha matrix ***
      alpha(:, :, :, :) = 0.0_dp
      DO iaxis = 1, 3
      DO lxa = 0, la_max_local
      DO lxb = 0, lb_max_local
         binomial_k_lxa = 1.0_dp
         a = 1.0_dp
         DO k = 0, lxa
            binomial_l_lxb = 1.0_dp
            b = 1.0_dp
            DO l = 0, lxb
               alpha(lxa-l+lxb-k, lxa, lxb, iaxis) = alpha(lxa-l+lxb-k, lxa, lxb, iaxis)+ &
                                                     binomial_k_lxa*binomial_l_lxb*a*b
               binomial_l_lxb = binomial_l_lxb*REAL(lxb-l, dp)/REAL(l+1, dp)
               b = b*(rp(iaxis)-(ra(iaxis)+rab(iaxis)))
            ENDDO
            binomial_k_lxa = binomial_k_lxa*REAL(lxa-k, dp)/REAL(k+1, dp)
            a = a*(-ra(iaxis)+rp(iaxis))
         ENDDO
      ENDDO
      ENDDO
      ENDDO

!
!   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
!   use a three step procedure
!   we don't store zeros, so counting is done using lxyz,lxy in order to have
!   contiguous memory access in collocate_fast.F
!
      lxyz = 0
      DO lzp = 0, lp
      DO lyp = 0, lp-lzp
      DO lxp = 0, lp-lzp-lyp
         lxyz = lxyz+1
         coef_xyz(lxyz) = 0.0_dp
      ENDDO
      ENDDO
      ENDDO
      DO lzb = 0, lb_max_local
      DO lza = 0, la_max_local
         lxy = 0
         DO lyp = 0, lp-lza-lzb
            DO lxp = 0, lp-lza-lzb-lyp
               lxy = lxy+1
               coef_xyt(lxy) = 0.0_dp
            ENDDO
            lxy = lxy+lza+lzb
         ENDDO
         DO lyb = 0, lb_max_local-lzb
         DO lya = 0, la_max_local-lza
            lxpm = (lb_max_local-lzb-lyb)+(la_max_local-lza-lya)
            coef_xtt(0:lxpm) = 0.0_dp
            DO lxb = MAX(lb_min_local-lzb-lyb, 0), lb_max_local-lzb-lyb
            DO lxa = MAX(la_min_local-lza-lya, 0), la_max_local-lza-lya
               ico = coset(lxa, lya, lza)
               jco = coset(lxb, lyb, lzb)
               p_ele = prefactor*pab_local(ico, jco)
               DO lxp = 0, lxa+lxb
                  coef_xtt(lxp) = coef_xtt(lxp)+p_ele*alpha(lxp, lxa, lxb, 1)
               ENDDO
            ENDDO
            ENDDO
            lxy = 0
            DO lyp = 0, lya+lyb
               DO lxp = 0, lp-lza-lzb-lya-lyb
                  lxy = lxy+1
                  coef_xyt(lxy) = coef_xyt(lxy)+alpha(lyp, lya, lyb, 2)*coef_xtt(lxp)
               ENDDO
               lxy = lxy+lza+lzb+lya+lyb-lyp
            ENDDO
         ENDDO
         ENDDO
         lxyz = 0
         DO lzp = 0, lza+lzb
            lxy = 0
            DO lyp = 0, lp-lza-lzb
               DO lxp = 0, lp-lza-lzb-lyp
                  lxy = lxy+1; lxyz = lxyz+1
                  coef_xyz(lxyz) = coef_xyz(lxyz)+alpha(lzp, lza, lzb, 3)*coef_xyt(lxy)
               ENDDO
               lxy = lxy+lza+lzb; lxyz = lxyz+lza+lzb-lzp
            ENDDO
            DO lyp = lp-lza-lzb+1, lp-lzp
               DO lxp = 0, lp-lyp-lzp
                  lxyz = lxyz+1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO

      IF (use_subpatch) THEN
         CALL collocate_general_opt()
      ELSE
         CALL collocate_ortho()
      END IF

!     IF (compute_tau) THEN
!        DEALLOCATE (pab_local)
!     END IF
      DEALLOCATE (pab_local)

   CONTAINS

      !
      ! this treats efficiently the orthogonal case
      !
! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_ortho()

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
            CPABORT("Called with too large radius.")
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

         ALLOCATE (map(-cmax:cmax, 3))
         cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
         roffset(:) = rp(:)-REAL(cubecenter(:), dp)*dr(:)
!   *** a mapping so that the ig corresponds to the right grid point
         DO i = 1, 3
            IF (perd(i) == 1) THEN
               start = lb_cube(i)
               DO
                  offset = MODULO(cubecenter(i)+start, ng(i))+1-start
                  length = MIN(ub_cube(i), ng(i)-offset)-start
                  DO ig = start, start+length
                     map(ig, i) = ig+offset
                  END DO
                  IF (start+length .GE. ub_cube(i)) EXIT
                  start = start+length+1
               END DO
            ELSE
               ! this takes partial grid + border regions into account
               offset = MODULO(cubecenter(i)+lb_cube(i)+lb_grid(i), ng(i))+1-lb_cube(i)
               ! check for out of bounds
               IF (ub_cube(i)+offset > UBOUND(grid, i) .OR. lb_cube(i)+offset < LBOUND(grid, i)) THEN
                  CPASSERT(.FALSE.)
               ENDIF
               DO ig = lb_cube(i), ub_cube(i)
                  map(ig, i) = ig+offset
               END DO
            END IF
         ENDDO
         ALLOCATE (pol_z(1:2, 0:lp, -cmax:0))
         ALLOCATE (pol_y(1:2, 0:lp, -cmax:0))
         ALLOCATE (pol_x(0:lp, -cmax:cmax))

!        prep.f90
!
!   compute the values of all (x-xp)**lp*exp(..)
!
! still requires the old trick:
!  new trick to avoid to many exps (reuse the result from the previous gridpoint):
!  exp( -a*(x+d)**2)=exp(-a*x**2)*(-2*a*x*d)*exp(-a*d**2)
!  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)

         iaxis = 3
         t_exp_1 = EXP(-zetp*dr(iaxis)**2)
         t_exp_2 = t_exp_1**2
         t_exp_min_1 = EXP(-zetp*(+dr(iaxis)-roffset(iaxis))**2)
         t_exp_min_2 = EXP(-2*zetp*(+dr(iaxis)-roffset(iaxis))*(-dr(iaxis)))
         t_exp_plus_1 = EXP(-zetp*(-roffset(iaxis))**2)
         t_exp_plus_2 = EXP(-2*zetp*(-roffset(iaxis))*(+dr(iaxis)))
         DO ig = 0, lb_cube(iaxis), -1
            rpg = REAL(ig, dp)*dr(iaxis)-roffset(iaxis)
            t_exp_min_1 = t_exp_min_1*t_exp_min_2*t_exp_1
            t_exp_min_2 = t_exp_min_2*t_exp_2
            pg = t_exp_min_1
            ! pg  = EXP(-zetp*rpg**2)
            DO icoef = 0, lp
               pol_z(1, icoef, ig) = pg
               pg = pg*(rpg)
            ENDDO

            rpg = REAL(1-ig, dp)*dr(iaxis)-roffset(iaxis)
            t_exp_plus_1 = t_exp_plus_1*t_exp_plus_2*t_exp_1
            t_exp_plus_2 = t_exp_plus_2*t_exp_2
            pg = t_exp_plus_1
            ! pg  = EXP(-zetp*rpg**2)
            DO icoef = 0, lp
               pol_z(2, icoef, ig) = pg
               pg = pg*(rpg)
            ENDDO
         ENDDO

         iaxis = 2
         t_exp_1 = EXP(-zetp*dr(iaxis)**2)
         t_exp_2 = t_exp_1**2
         t_exp_min_1 = EXP(-zetp*(+dr(iaxis)-roffset(iaxis))**2)
         t_exp_min_2 = EXP(-2*zetp*(+dr(iaxis)-roffset(iaxis))*(-dr(iaxis)))
         t_exp_plus_1 = EXP(-zetp*(-roffset(iaxis))**2)
         t_exp_plus_2 = EXP(-2*zetp*(-roffset(iaxis))*(+dr(iaxis)))
         DO ig = 0, lb_cube(iaxis), -1
            rpg = REAL(ig, dp)*dr(iaxis)-roffset(iaxis)
            t_exp_min_1 = t_exp_min_1*t_exp_min_2*t_exp_1
            t_exp_min_2 = t_exp_min_2*t_exp_2
            pg = t_exp_min_1
            ! pg  = EXP(-zetp*rpg**2)
            DO icoef = 0, lp
               pol_y(1, icoef, ig) = pg
               pg = pg*(rpg)
            ENDDO

            rpg = REAL(1-ig, dp)*dr(iaxis)-roffset(iaxis)
            t_exp_plus_1 = t_exp_plus_1*t_exp_plus_2*t_exp_1
            t_exp_plus_2 = t_exp_plus_2*t_exp_2
            pg = t_exp_plus_1
            ! pg  = EXP(-zetp*rpg**2)
            DO icoef = 0, lp
               pol_y(2, icoef, ig) = pg
               pg = pg*(rpg)
            ENDDO
         ENDDO

         iaxis = 1
         t_exp_1 = EXP(-zetp*dr(iaxis)**2)
         t_exp_2 = t_exp_1**2
         t_exp_min_1 = EXP(-zetp*(+dr(iaxis)-roffset(iaxis))**2)
         t_exp_min_2 = EXP(-2*zetp*(+dr(iaxis)-roffset(iaxis))*(-dr(iaxis)))
         t_exp_plus_1 = EXP(-zetp*(-roffset(iaxis))**2)
         t_exp_plus_2 = EXP(-2*zetp*(-roffset(iaxis))*(+dr(iaxis)))
         DO ig = 0, lb_cube(1), -1

            rpg = REAL(ig, dp)*dr(1)-roffset(1)
            t_exp_min_1 = t_exp_min_1*t_exp_min_2*t_exp_1
            t_exp_min_2 = t_exp_min_2*t_exp_2
            pg = t_exp_min_1
            !pg  = EXP(-zetp*rpg**2)
            DO icoef = 0, lp
               pol_x(icoef, ig) = pg
               pg = pg*(rpg)
            ENDDO

            rpg = REAL(1-ig, dp)*dr(1)-roffset(1)
            t_exp_plus_1 = t_exp_plus_1*t_exp_plus_2*t_exp_1
            t_exp_plus_2 = t_exp_plus_2*t_exp_2
            pg = t_exp_plus_1
            ! pg  = EXP(-zetp*rpg**2)
            DO icoef = 0, lp
               pol_x(icoef, 1-ig) = pg
               pg = pg*(rpg)
            ENDDO
         ENDDO

         CALL collocate_core_default(grid, coef_xyz, pol_x(0:, -cmax:), pol_y(1:, 0:, -cmax:), pol_z(1:, 0:, -cmax:), &
                                     map(-cmax:, 1:), sphere_bounds, lp, cmax, gridbounds)

         ! deallocation needed to pass around a pgi bug..
         DEALLOCATE (pol_z)
         DEALLOCATE (pol_y)
         DEALLOCATE (pol_x)
         DEALLOCATE (map)

      END SUBROUTINE collocate_ortho

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
      SUBROUTINE collocate_general_opt()

      INTEGER :: i, i_index, il, ilx, ily, ilz, index_max(3), index_min(3), ismax, ismin, j, &
         j_index, jl, jlx, jly, jlz, k, k_index, kl, klx, kly, klz, lpx, lpy, lpz, lx, ly, lz, &
         offset(3)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: grid_map
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: coef_map
      REAL(KIND=dp)                                      :: a, b, c, d, di, dip, dj, djp, dk, dkp, &
                                                            exp0i, exp1i, exp2i, gp(3), &
                                                            hmatgrid(3, 3), point(3), pointj(3), &
                                                            pointk(3), res, resc(3), v(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coef_ijk
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: hmatgridp

!
! transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
! sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
! sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
!

         ALLOCATE (coef_ijk(((lp+1)*(lp+2)*(lp+3))/6))

         ! aux mapping array to simplify life
         ALLOCATE (coef_map(0:lp, 0:lp, 0:lp))
         coef_map = HUGE(coef_map)
         lxyz = 0
         DO lzp = 0, lp
         DO lyp = 0, lp-lzp
         DO lxp = 0, lp-lzp-lyp
            lxyz = lxyz+1
            coef_ijk(lxyz) = 0.0_dp
            coef_map(lxp, lyp, lzp) = lxyz
         ENDDO
         ENDDO
         ENDDO

         ! cell hmat in grid points
         hmatgrid = dh

         ! center in grid coords
         gp = MATMUL(dh_inv, rp)
         cubecenter(:) = FLOOR(gp)

         ! transform using multinomials
         ALLOCATE (hmatgridp(3, 3, 0:lp))
         hmatgridp(:, :, 0) = 1.0_dp
         DO k = 1, lp
            hmatgridp(:, :, k) = hmatgridp(:, :, k-1)*hmatgrid(:, :)
         ENDDO

         lpx = lp
         DO klx = 0, lpx
         DO jlx = 0, lpx-klx
         DO ilx = 0, lpx-klx-jlx
            lx = ilx+jlx+klx
            lpy = lp-lx
            DO kly = 0, lpy
            DO jly = 0, lpy-kly
            DO ily = 0, lpy-kly-jly
               ly = ily+jly+kly
               lpz = lp-lx-ly
               DO klz = 0, lpz
               DO jlz = 0, lpz-klz
               DO ilz = 0, lpz-klz-jlz
                  lz = ilz+jlz+klz

                  il = ilx+ily+ilz
                  jl = jlx+jly+jlz
                  kl = klx+kly+klz
                  coef_ijk(coef_map(il, jl, kl)) = &
                     coef_ijk(coef_map(il, jl, kl))+coef_xyz(coef_map(lx, ly, lz))* &
                     hmatgridp(1, 1, ilx)*hmatgridp(1, 2, jlx)*hmatgridp(1, 3, klx)* &
                     hmatgridp(2, 1, ily)*hmatgridp(2, 2, jly)*hmatgridp(2, 3, kly)* &
                     hmatgridp(3, 1, ilz)*hmatgridp(3, 2, jlz)*hmatgridp(3, 3, klz)* &
                     fac(lx)*fac(ly)*fac(lz)/ &
                     (fac(ilx)*fac(ily)*fac(ilz)*fac(jlx)*fac(jly)*fac(jlz)*fac(klx)*fac(kly)*fac(klz))
               ENDDO
               ENDDO
               ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDDO

         IF (radius > max_rad_ga) THEN
            !
            ! This is an important check. If the required radius for mapping the density is different
            ! from the actual computed one, (significant) errors can occur.
            ! This error can invariably be fixed by improving the computation of maxradius
            ! in the call to init_cube_info
            !
            WRITE (*, *) max_rad_ga, radius
            CPABORT("Called with too large radius.")
         ENDIF

         ! get the min max indices that contain at least the cube that contains a sphere around rp of radius radius
         ! if the cell is very non-orthogonal this implies that many useless points are included
         ! this estimate can be improved (i.e. not box but sphere should be used)
         index_min = HUGE(index_min)
         index_max = -HUGE(index_max)
         DO i = -1, 1
         DO j = -1, 1
         DO k = -1, 1
            point(1) = rp(1)+i*radius
            point(2) = rp(2)+j*radius
            point(3) = rp(3)+k*radius
            resc(1) = dh_inv(1, 1)*point(1)+dh_inv(1, 2)*point(2)+dh_inv(1, 3)*point(3)
            resc(2) = dh_inv(2, 1)*point(1)+dh_inv(2, 2)*point(2)+dh_inv(2, 3)*point(3)
            resc(3) = dh_inv(3, 1)*point(1)+dh_inv(3, 2)*point(2)+dh_inv(3, 3)*point(3)
            index_min = MIN(index_min, FLOOR(resc))
            index_max = MAX(index_max, CEILING(resc))
         ENDDO
         ENDDO
         ENDDO

         offset(:) = MODULO(index_min(:)+lb_grid(:), ng(:))+1

         ALLOCATE (grid_map(index_min(1):index_max(1)))
         DO i = index_min(1), index_max(1)
            grid_map(i) = MODULO(i, ng(1))+1
            IF (perd(1) == 1) THEN
               grid_map(i) = MODULO(i, ng(1))+1
            ELSE
               grid_map(i) = i-index_min(1)+offset(1)
            ENDIF
         ENDDO

         ! go over the grid, but cycle if the point is not within the radius
         DO k = index_min(3), index_max(3)
            dk = k-gp(3)
            pointk = hmatgrid(:, 3)*dk

            IF (perd(3) == 1) THEN
               k_index = MODULO(k, ng(3))+1
            ELSE
               k_index = k-index_min(3)+offset(3)
            ENDIF

            coef_xyt = 0.0_dp
            lxyz = 0
            dkp = 1.0_dp
            DO kl = 0, lp
               lxy = 0
               DO jl = 0, lp-kl
                  DO il = 0, lp-kl-jl
                     lxyz = lxyz+1; lxy = lxy+1
                     coef_xyt(lxy) = coef_xyt(lxy)+coef_ijk(lxyz)*dkp
                  ENDDO
                  lxy = lxy+kl
               ENDDO
               dkp = dkp*dk
            ENDDO

            DO j = index_min(2), index_max(2)
               dj = j-gp(2)
               pointj = pointk+hmatgrid(:, 2)*dj
               IF (perd(2) == 1) THEN
                  j_index = MODULO(j, ng(2))+1
               ELSE
                  j_index = j-index_min(2)+offset(2)
               ENDIF

               coef_xtt = 0.0_dp
               lxy = 0
               djp = 1.0_dp
               DO jl = 0, lp
                  DO il = 0, lp-jl
                     lxy = lxy+1
                     coef_xtt(il) = coef_xtt(il)+coef_xyt(lxy)*djp
                  ENDDO
                  djp = djp*dj
               ENDDO

               ! find bounds for the inner loop
               ! based on a quadratic equation in i
               ! a*i**2+b*i+c=radius**2
               v = pointj-gp(1)*hmatgrid(:, 1)
               a = DOT_PRODUCT(hmatgrid(:, 1), hmatgrid(:, 1))
               b = 2*DOT_PRODUCT(v, hmatgrid(:, 1))
               c = DOT_PRODUCT(v, v)
               d = b*b-4*a*(c-radius**2)

               IF (d < 0) THEN
                  CYCLE
               ELSE
                  d = SQRT(d)
                  ismin = CEILING((-b-d)/(2*a))
                  ismax = FLOOR((-b+d)/(2*a))
               ENDIF
               ! prepare for computing -zetp*rsq
               a = -zetp*a
               b = -zetp*b
               c = -zetp*c
               i = ismin-1

               ! the recursion relation might have to be done
               ! from the center of the gaussian (in both directions)
               ! instead as the current implementation from an edge
               exp2i = EXP((a*i+b)*i+c)
               exp1i = EXP(2*a*i+a+b)
               exp0i = EXP(2*a)

               DO i = ismin, ismax
                  di = i-gp(1)

                  ! polynomial terms
                  res = 0.0_dp
                  dip = 1.0_dp
                  DO il = 0, lp
                     res = res+coef_xtt(il)*dip
                     dip = dip*di
                  ENDDO

                  ! the exponential recursion
                  exp2i = exp2i*exp1i
                  exp1i = exp1i*exp0i
                  res = res*exp2i

                  i_index = grid_map(i)
                  grid(i_index, j_index, k_index) = grid(i_index, j_index, k_index)+res
               ENDDO
            ENDDO
         ENDDO
         !t2=nanotime_ia32()
         !write(*,*) t2-t1
         ! deallocation needed to pass around a pgi bug..
         DEALLOCATE (coef_ijk)
         DEALLOCATE (coef_map)
         DEALLOCATE (hmatgridp)
         DEALLOCATE (grid_map)

      END SUBROUTINE collocate_general_opt

   END SUBROUTINE collocate_pgf_product_rspace

! **************************************************************************************************
!> \brief ...
!> \param grid ...
!> \param coef_xyz ...
!> \param pol_x ...
!> \param pol_y ...
!> \param pol_z ...
!> \param map ...
!> \param sphere_bounds ...
!> \param lp ...
!> \param cmax ...
!> \param gridbounds ...
! **************************************************************************************************
   SUBROUTINE collocate_core_default(grid, coef_xyz, pol_x, pol_y, pol_z, map, sphere_bounds, lp, cmax, gridbounds)
      INTEGER, INTENT(IN)                      :: sphere_bounds(:), lp
      REAL(dp), INTENT(IN) :: coef_xyz(((lp+1)*(lp+2)*(lp+3))/6)
      INTEGER, INTENT(IN)                      :: cmax
      REAL(dp), INTENT(IN)                     :: pol_x(0:lp, -cmax:cmax), &
                                                  pol_y(1:2, 0:lp, -cmax:0), &
                                                  pol_z(1:2, 0:lp, -cmax:0)
      INTEGER, INTENT(IN)                      :: map(-cmax:cmax, 1:3), &
                                                  gridbounds(2, 3)
      REAL(dp), INTENT(INOUT) :: grid(gridbounds(1, 1):gridbounds(2, 1), &
                                      gridbounds(1, 2):gridbounds(2, 2), gridbounds(1, 3):gridbounds(2, 3))

      INTEGER                                  :: i, ig, igmax, igmin, j, j2, &
                                                  jg, jg2, jgmin, k, k2, kg, &
                                                  kg2, kgmin, lxp, lxy, lxyz, &
                                                  lyp, lzp, sci
      REAL(dp)                                 :: coef_x(4, 0:lp), &
                                                  coef_xy(2, (lp+1)*(lp+2)/2), &
                                                  s01, s02, s03, s04

      sci = 1

      kgmin = sphere_bounds(sci)
      sci = sci+1
      DO kg = kgmin, 0
         kg2 = 1-kg
         k = map(kg, 3)
         k2 = map(kg2, 3)

         coef_xy = 0.0_dp
         lxyz = 0
         DO lzp = 0, lp
            lxy = 0
            DO lyp = 0, lp-lzp
               DO lxp = 0, lp-lzp-lyp
                  lxyz = lxyz+1; lxy = lxy+1
                  coef_xy(1, lxy) = coef_xy(1, lxy)+coef_xyz(lxyz)*pol_z(1, lzp, kg)
                  coef_xy(2, lxy) = coef_xy(2, lxy)+coef_xyz(lxyz)*pol_z(2, lzp, kg)
               ENDDO
               lxy = lxy+lzp
            ENDDO
         ENDDO

         jgmin = sphere_bounds(sci)
         sci = sci+1
         DO jg = jgmin, 0
            jg2 = 1-jg
            j = map(jg, 2)
            j2 = map(jg2, 2)
            igmin = sphere_bounds(sci)
            sci = sci+1
            igmax = 1-igmin

            coef_x = 0.0_dp
            lxy = 0
            DO lyp = 0, lp
            DO lxp = 0, lp-lyp
               lxy = lxy+1
               coef_x(1, lxp) = coef_x(1, lxp)+coef_xy(1, lxy)*pol_y(1, lyp, jg)
               coef_x(2, lxp) = coef_x(2, lxp)+coef_xy(2, lxy)*pol_y(1, lyp, jg)
               coef_x(3, lxp) = coef_x(3, lxp)+coef_xy(1, lxy)*pol_y(2, lyp, jg)
               coef_x(4, lxp) = coef_x(4, lxp)+coef_xy(2, lxy)*pol_y(2, lyp, jg)
            ENDDO
            ENDDO

            DO ig = igmin, igmax
               i = map(ig, 1)
               s01 = 0.0_dp
               s02 = 0.0_dp
               s03 = 0.0_dp
               s04 = 0.0_dp
               DO lxp = 0, lp
                  s01 = s01+coef_x(1, lxp)*pol_x(lxp, ig)
                  s02 = s02+coef_x(2, lxp)*pol_x(lxp, ig)
                  s03 = s03+coef_x(3, lxp)*pol_x(lxp, ig)
                  s04 = s04+coef_x(4, lxp)*pol_x(lxp, ig)
               ENDDO
               grid(i, j, k) = grid(i, j, k)+s01
               grid(i, j2, k) = grid(i, j2, k)+s03
               grid(i, j, k2) = grid(i, j, k2)+s02
               grid(i, j2, k2) = grid(i, j2, k2)+s04
            END DO

         END DO
      END DO

   END SUBROUTINE collocate_core_default

! **************************************************************************************************

! **************************************************************************************************
!> \brief ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
!> \return ...
! **************************************************************************************************
   FUNCTION coset(lx, ly, lz) RESULT(cval)
      INTEGER, INTENT(IN)                                :: lx, ly, lz
      INTEGER                                            :: cval

      INTEGER                                            :: l

      l = lx+ly+lz
      cval = ncoset(l-1)+1+((l-lx)*(l-lx+1))/2+lz

   END FUNCTION coset

! **************************************************************************************************

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
!> \param zeta ...
!> \param zetb ...
! **************************************************************************************************
   SUBROUTINE prepare_pab_tau(pab_tau, pab, o1, o2, &
                              la_max, la_min, lb_max, lb_min, zeta, zetb)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab_tau
      REAL(dp), DIMENSION(:, :), INTENT(in)              :: pab
      INTEGER, INTENT(in)                                :: o1, o2, la_max, la_min, lb_max, lb_min
      REAL(KIND=dp), INTENT(in)                          :: zeta, zetb

      INTEGER                                            :: ico, ico_l, jco, jco_l, lxa, lxb, lya, &
                                                            lyb, lza, lzb, nla, nlb

      ! create a new pab_tau so that mapping pab_tau with pgf_a pgf_b
      ! is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
      ! (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})
      ! cleaner would possibly be to touch pzyx directly (avoiding the following allocate)
      nla = ncoset(la_max+1)
      nlb = ncoset(lb_max+1)
      ALLOCATE (pab_tau(nla, nlb))
      pab_tau = 0.0_dp
      DO lxa = 0, la_max
      DO lxb = 0, lb_max
         DO lya = 0, la_max-lxa
         DO lyb = 0, lb_max-lxb
            DO lza = MAX(la_min-lxa-lya, 0), la_max-lxa-lya
            DO lzb = MAX(lb_min-lxb-lyb, 0), lb_max-lxb-lyb

               ico = coset(lxa, lya, lza)
               jco = coset(lxb, lyb, lzb)

               ! x  (all safe if lxa = 0, as the spurious added terms have zero prefactor)

               ico_l = coset(MAX(lxa-1, 0), lya, lza)
               jco_l = coset(MAX(lxb-1, 0), lyb, lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)+lxa*lxb*pab(o1+ico, o2+jco)
               ico_l = coset(MAX(lxa-1, 0), lya, lza)
               jco_l = coset((lxb+1), lyb, lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)-2.0_dp*lxa*zetb*pab(o1+ico, o2+jco)
               ico_l = coset((lxa+1), lya, lza)
               jco_l = coset(MAX(lxb-1, 0), lyb, lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)-2.0_dp*zeta*lxb*pab(o1+ico, o2+jco)
               ico_l = coset((lxa+1), lya, lza)
               jco_l = coset((lxb+1), lyb, lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)+4.0_dp*zeta*zetb*pab(o1+ico, o2+jco)

               ! y

               ico_l = coset(lxa, MAX(lya-1, 0), lza)
               jco_l = coset(lxb, MAX(lyb-1, 0), lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)+lya*lyb*pab(o1+ico, o2+jco)
               ico_l = coset(lxa, MAX(lya-1, 0), lza)
               jco_l = coset(lxb, (lyb+1), lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)-2.0_dp*lya*zetb*pab(o1+ico, o2+jco)
               ico_l = coset(lxa, (lya+1), lza)
               jco_l = coset(lxb, MAX(lyb-1, 0), lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)-2.0_dp*zeta*lyb*pab(o1+ico, o2+jco)
               ico_l = coset(lxa, (lya+1), lza)
               jco_l = coset(lxb, (lyb+1), lzb)
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)+4.0_dp*zeta*zetb*pab(o1+ico, o2+jco)

               ! z

               ico_l = coset(lxa, lya, MAX(lza-1, 0))
               jco_l = coset(lxb, lyb, MAX(lzb-1, 0))
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)+lza*lzb*pab(o1+ico, o2+jco)
               ico_l = coset(lxa, lya, MAX(lza-1, 0))
               jco_l = coset(lxb, lyb, (lzb+1))
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)-2.0_dp*lza*zetb*pab(o1+ico, o2+jco)
               ico_l = coset(lxa, lya, (lza+1))
               jco_l = coset(lxb, lyb, MAX(lzb-1, 0))
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)-2.0_dp*zeta*lzb*pab(o1+ico, o2+jco)
               ico_l = coset(lxa, lya, (lza+1))
               jco_l = coset(lxb, lyb, (lzb+1))
               pab_tau(ico_l, jco_l) = pab_tau(ico_l, jco_l)+4.0_dp*zeta*zetb*pab(o1+ico, o2+jco)

            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      pab_tau = 0.5_dp*pab_tau

   END SUBROUTINE prepare_pab_tau

! **************************************************************************************************

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
! **************************************************************************************************
   SUBROUTINE prepare_pab_rho(pab_rho, pab, o1, o2, la_max, la_min, lb_max, lb_min)
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab_rho
      REAL(dp), DIMENSION(:, :), INTENT(in)              :: pab
      INTEGER, INTENT(in)                                :: o1, o2, la_max, la_min, lb_max, lb_min

      INTEGER                                            :: ico, jco, lxa, lxb, lya, lyb, lza, lzb, &
                                                            nla, nlb

      nla = ncoset(la_max)
      nlb = ncoset(lb_max)
      ALLOCATE (pab_rho(nla, nlb))
      pab_rho = 0.0_dp
      DO lxa = 0, la_max
      DO lxb = 0, lb_max
         DO lya = 0, la_max-lxa
         DO lyb = 0, lb_max-lxb
            DO lza = MAX(la_min-lxa-lya, 0), la_max-lxa-lya
            DO lzb = MAX(lb_min-lxb-lyb, 0), lb_max-lxb-lyb
               ico = coset(lxa, lya, lza)
               jco = coset(lxb, lyb, lzb)
               pab_rho(ico, jco) = pab(o1+ico, o2+jco)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
      ENDDO

   END SUBROUTINE prepare_pab_rho

END MODULE grid_base_ref

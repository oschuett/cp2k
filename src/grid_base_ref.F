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
                                              C_NULL_PTR,&
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

   INTEGER, PARAMETER :: maxl = 18
   INTEGER, PARAMETER, DIMENSION(-1:maxl) :: ncoset = (/ &
                                         0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455, 560, 680, 816, 969, 1140, 1330/)

   PUBLIC :: grid_collocate_ref

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

      INTEGER :: bcol, brow, iatom, iatom_old, igrid_level, ikind, img, img_old, imr, ipair, ipgf, &
         iset, iset_old, itask, jatom, jatom_old, jkind, jpgf, jset, jset_old, max_radius, maxco, &
         maxpgf, maxset, maxsgf_set, n1, n2, na1, natoms, nb1, ncoa, ncob, nimages, ntasks, sgfa, &
         sgfb
      INTEGER, DIMENSION(3)                              :: lb_grid(3), npts(3), perd(3)
      INTEGER, DIMENSION(:), POINTER                     :: first_sgfa, first_sgfb, la_max, la_min, &
                                                            lb_max, lb_min, npgfa, npgfb, nsgfa, &
                                                            nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: glb_cube, gsphere_bounds, gub_cube
      LOGICAL                                            :: atom_pair_changed, orthorhombic, &
                                                            radius_inconsistent, use_ortho
      REAL(KIND=dp)                                      :: adh, drmin, max_rad_ga, radius, rscale, &
                                                            zetp
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

         glb_cube => task_list%cinfo(igrid_level)%lb_cube
         gub_cube => task_list%cinfo(igrid_level)%ub_cube
         gsphere_bounds => task_list%cinfo(igrid_level)%sphere_bounds
         orthorhombic = task_list%cinfo(igrid_level)%orthorhombic

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

               !TODO: Only needed for tests/QS/regtest-kp-1/c_gapw*.inp, see below.
               IF (iatom > jatom) THEN
                  pab = TRANSPOSE(pab)
               ENDIF

               iset_old = iset
               jset_old = jset
            ENDIF

            rb(:) = ra(:) + rab(:)
            zetp = zeta(ipgf, iset) + zetb(jpgf, jset)
            imr = MAX(1, CEILING(radius/drmin))

            IF (2.0_dp*radius < adh) CYCLE

            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1

            ! takes the density matrix symmetry in account
            IF (iatom == jatom) THEN
               rscale = 1.0_dp
            ELSE
               rscale = 2.0_dp
            END IF

            ! check whether we need to use the generalised collocation scheme
            IF (tasks(itask)%dist_type == 2) THEN
               use_ortho = .FALSE.
               radius_inconsistent = radius > max_rad_ga
            ELSE
               use_ortho = orthorhombic
               radius_inconsistent = imr > max_radius
            ENDIF

            IF (radius_inconsistent) THEN
               !
               ! This is an important check. If the required radius for mapping the density is different
               ! from the actual computed one, (significant) errors can occur.
               ! This error can invariably be fixed by improving the computation of maxradius
               ! in the call to init_cube_info
               !
               CPABORT("Radius consistency check failed.")
            ENDIF

            !TODO: IF/ELSE only needed for tests/QS/regtest-kp-1/c_gapw*.inp, see above.
            IF (iatom <= jatom) THEN
               CALL collocate_pgf_product_rspace(la_max=la_max(iset), &
                                                 zeta=zeta(ipgf, iset), &
                                                 la_min=la_min(iset), &
                                                 lb_max=lb_max(jset), &
                                                 zetb=zetb(jpgf, jset), &
                                                 lb_min=lb_min(jset), &
                                                 ra=ra, &
                                                 rab=rab, &
                                                 radius=radius, &
                                                 rscale=rscale, &
                                                 pab=pab, &
                                                 o1=na1 - 1, &
                                                 o2=nb1 - 1, &
                                                 grid=rs_rho(igrid_level)%grid, &
                                                 dh=dh, &
                                                 dh_inv=dh_inv, &
                                                 perd=perd, &
                                                 npts=npts, &
                                                 lb_grid=lb_grid, &
                                                 glb_cube=glb_cube, &
                                                 gub_cube=gub_cube, &
                                                 gsphere_bounds=gsphere_bounds, &
                                                 imr=imr, &
                                                 compute_tau=compute_tau, &
                                                 use_ortho=use_ortho, &
                                                 lmax_global=lmax_global, &
                                                 maxco=maxco)
            ELSE
               CALL collocate_pgf_product_rspace(la_max=lb_max(jset), &
                                                 zeta=zetb(jpgf, jset), &
                                                 la_min=lb_min(jset), &
                                                 lb_max=la_max(iset), &
                                                 zetb=zeta(ipgf, iset), &
                                                 lb_min=la_min(iset), &
                                                 ra=rb, &
                                                 rab=-rab, &
                                                 radius=radius, &
                                                 rscale=rscale, &
                                                 pab=pab, &
                                                 o1=nb1 - 1, &
                                                 o2=na1 - 1, &
                                                 grid=rs_rho(igrid_level)%grid, &
                                                 dh=dh, &
                                                 dh_inv=dh_inv, &
                                                 perd=perd, &
                                                 npts=npts, &
                                                 lb_grid=lb_grid, &
                                                 glb_cube=glb_cube, &
                                                 gub_cube=gub_cube, &
                                                 gsphere_bounds=gsphere_bounds, &
                                                 imr=imr, &
                                                 compute_tau=compute_tau, &
                                                 use_ortho=use_ortho, &
                                                 lmax_global=lmax_global, &
                                                 maxco=maxco)
            ENDIF

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
!> \param glb_cube ...
!> \param gub_cube ...
!> \param gsphere_bounds ...
!> \param imr ...
!> \param compute_tau ...
!> \param use_ortho ...
!> \param lmax_global Maximum possible value of lmax used to dimension arrays
!> \param maxco ...
! **************************************************************************************************
   SUBROUTINE collocate_pgf_product_rspace(la_max, zeta, la_min, &
                                           lb_max, zetb, lb_min, &
                                           ra, rab, radius, rscale, pab, o1, o2, &
                                           grid, &
                                           dh, dh_inv, perd, npts, lb_grid, &
                                           glb_cube, gub_cube, gsphere_bounds, imr, &
                                           compute_tau, use_ortho, lmax_global, maxco)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rab
      REAL(KIND=dp), INTENT(IN)                          :: radius, rscale
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN), TARGET :: dh, dh_inv
      INTEGER, DIMENSION(3), INTENT(IN), TARGET          :: perd, npts, lb_grid
      INTEGER, DIMENSION(:, :), POINTER                  :: glb_cube, gub_cube, gsphere_bounds
      INTEGER, INTENT(IN)                                :: imr
      LOGICAL, INTENT(IN)                                :: compute_tau, use_ortho
      INTEGER, INTENT(IN)                                :: lmax_global, maxco

      INTEGER                                            :: nspheres
      INTEGER, DIMENSION(3), TARGET                      :: ngrid
      LOGICAL(KIND=C_BOOL), DIMENSION(3), TARGET         :: periodic
      TYPE(C_PTR)                                        :: lb_cube, sphere_bounds, ub_cube
      INTERFACE
         SUBROUTINE grid_collocate_pgf_product_rspace(compute_tau, use_ortho, &
                                                      la_max, la_min, lb_max, lb_min, &
                                                      zeta, zetb, rscale, &
                                                      dh, dh_inv, ra, rab, npts, ngrid, lb_grid, periodic, &
                                                      lmax, radius, &
                                                      lb_cube, ub_cube, nspheres, sphere_bounds, &
                                                      maxco, o1, o2, pab, &
                                                      grid) &
            BIND(C, name="grid_collocate_pgf_product_rspace")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: compute_tau
            LOGICAL(KIND=C_BOOL), VALUE               :: use_ortho
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
            INTEGER(KIND=C_INT), VALUE                :: lmax
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            TYPE(C_PTR), VALUE                        :: lb_cube
            TYPE(C_PTR), VALUE                        :: ub_cube
            INTEGER(KIND=C_INT), VALUE                :: nspheres
            TYPE(C_PTR), VALUE                        :: sphere_bounds
            INTEGER(KIND=C_INT), VALUE                :: maxco
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: grid
         END SUBROUTINE grid_collocate_pgf_product_rspace
      END INTERFACE

      IF (use_ortho) THEN
         lb_cube = C_LOC(glb_cube(1, imr))
         ub_cube = C_LOC(gub_cube(1, imr))
         CPASSERT(SIZE(glb_cube, 1) == 3 .AND. SIZE(gub_cube, 1) == 3)
         sphere_bounds = C_LOC(gsphere_bounds(1, imr))
         nspheres = SIZE(gsphere_bounds, 1)
      ELSE
         lb_cube = C_NULL_PTR
         ub_cube = C_NULL_PTR
         sphere_bounds = C_NULL_PTR
         nspheres = 0
      ENDIF

      periodic = (perd == 1)

      CPASSERT(SIZE(pab, 1) == maxco .AND. LBOUND(pab, 1) == 1)
      CPASSERT(SIZE(pab, 2) == maxco .AND. LBOUND(pab, 2) == 1)

      ! see rs_grid_create() in ./src/pw/realspace_grid_types.F
      CPASSERT(LBOUND(grid, 1) == 1)
      CPASSERT(LBOUND(grid, 2) == 1)
      CPASSERT(LBOUND(grid, 3) == 1)
      ngrid = (/SIZE(grid, 1), SIZE(grid, 2), SIZE(grid, 3)/)

      !TODO: enable for gfortran 9
      !CPASSERT(IS_CONTIGUOUS(dh))
      !CPASSERT(IS_CONTIGUOUS(dh_inv))
      !CPASSERT(IS_CONTIGUOUS(ra))
      !CPASSERT(IS_CONTIGUOUS(rab))
      !CPASSERT(IS_CONTIGUOUS(npts))
      !CPASSERT(IS_CONTIGUOUS(lb_grid))
      !CPASSERT(IS_CONTIGUOUS(periodic))
      !CPASSERT(IS_CONTIGUOUS(pab))
      !CPASSERT(IS_CONTIGUOUS(grid))

      CALL grid_collocate_pgf_product_rspace(compute_tau=LOGICAL(compute_tau, KIND=C_BOOL), &
                                             use_ortho=LOGICAL(use_ortho, KIND=C_BOOL), &
                                             la_max=la_max, &
                                             la_min=la_min, &
                                             lb_max=lb_max, &
                                             lb_min=lb_min, &
                                             zeta=zeta, &
                                             zetb=zetb, &
                                             rscale=rscale, &
                                             dh=C_LOC(dh(1, 1)), &
                                             dh_inv=C_LOC(dh_inv(1, 1)), &
                                             ra=C_LOC(ra(1)), &
                                             rab=C_LOC(rab(1)), &
                                             npts=C_LOC(npts(1)), &
                                             ngrid=C_LOC(ngrid(1)), &
                                             lb_grid=C_LOC(lb_grid(1)), &
                                             periodic=C_LOC(periodic(1)), &
                                             lmax=lmax_global, &
                                             radius=radius, &
                                             lb_cube=lb_cube, &
                                             ub_cube=ub_cube, &
                                             nspheres=nspheres, &
                                             sphere_bounds=sphere_bounds, &
                                             maxco=maxco, &
                                             o1=o1, &
                                             o2=o2, &
                                             pab=C_LOC(pab(1, 1)), &
                                             grid=C_LOC(grid(1, 1, 1)))

   END SUBROUTINE collocate_pgf_product_rspace

END MODULE grid_base_ref

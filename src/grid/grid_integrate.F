!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE grid_integrate
   USE ISO_C_BINDING,                   ONLY: C_BOOL,&
                                              C_DOUBLE,&
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_PTR,&
                                              C_PTR
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE kinds,                           ONLY: dp
   USE realspace_grid_types,            ONLY: realspace_grid_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_integrate'

! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

   PUBLIC :: integrate_pgf_product

CONTAINS

! **************************************************************************************************
!> \brief Helper routines for getting rsgrid properties and asserting underlying assumptions.
!> \param rsgrid ...
!> \param npts_global ...
!> \param npts_local ...
!> \param shift_local ...
!> \param border_width ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE get_rsgrid_properties(rsgrid, npts_global, npts_local, shift_local, border_width)
      TYPE(realspace_grid_type), INTENT(IN)              :: rsgrid
      INTEGER, DIMENSION(:)                              :: npts_global, npts_local, shift_local, &
                                                            border_width

      INTEGER                                            :: i

      ! See rs_grid_create() in ./src/pw/realspace_grid_types.F.
      CPASSERT(LBOUND(rsgrid%r, 1) == rsgrid%lb_local(1))
      CPASSERT(UBOUND(rsgrid%r, 1) == rsgrid%ub_local(1))
      CPASSERT(LBOUND(rsgrid%r, 2) == rsgrid%lb_local(2))
      CPASSERT(UBOUND(rsgrid%r, 2) == rsgrid%ub_local(2))
      CPASSERT(LBOUND(rsgrid%r, 3) == rsgrid%lb_local(3))
      CPASSERT(UBOUND(rsgrid%r, 3) == rsgrid%ub_local(3))

      ! While the rsgrid code assumes that the grid starts at rsgrid%lb,
      ! the collocate code assumes that the grid starts at (1,1,1) in Fortran, or (0,0,0) in C.
      ! So, a point rp(:) gets the following grid coordinates MODULO(rp(:)/dr(:),npts_global(:))

      ! Number of global grid points in each direction.
      npts_global = rsgrid%desc%ub - rsgrid%desc%lb + 1

      ! Number of local grid points in each direction.
      npts_local = rsgrid%ub_local - rsgrid%lb_local + 1

      ! Number of points the local grid is shifted wrt global grid.
      shift_local = rsgrid%lb_local - rsgrid%desc%lb

      ! Convert rsgrid%desc%border and rsgrid%desc%perd into the more convenient border_width array.
      DO i = 1, 3
         IF (rsgrid%desc%perd(i) == 1) THEN
            ! Periodic meaning the grid in this direction is entriely present on every processor.
            CPASSERT(npts_local(i) == npts_global(i))
            CPASSERT(shift_local(i) == 0)
            ! No need for halo regions.
            border_width(i) = 0
         ELSE
            ! Not periodic meaning the grid in this direction is distributed among processors.
            CPASSERT(npts_local(i) <= npts_global(i))
            ! Check bounds of grid section that is owned by this processor.
            CPASSERT(rsgrid%lb_real(i) == rsgrid%lb_local(i) + rsgrid%desc%border)
            CPASSERT(rsgrid%ub_real(i) == rsgrid%ub_local(i) - rsgrid%desc%border)
            ! We have halo regions.
            border_width(i) = rsgrid%desc%border
         ENDIF
      ENDDO
   END SUBROUTINE get_rsgrid_properties

! **************************************************************************************************
!> \brief low level function to compute matrix elements of primitive gaussian functions
!> \param la_max ...
!> \param zeta ...
!> \param la_min ...
!> \param lb_max ...
!> \param zetb ...
!> \param lb_min ...
!> \param ra ...
!> \param rab ...
!> \param rsgrid ...
!> \param cell ...
!> \param cube_info ...
!> \param hab ...
!> \param pab ...
!> \param o1 ...
!> \param o2 ...
!> \param radius ...
!> \param calculate_forces ...
!> \param force_a ...
!> \param force_b ...
!> \param compute_tau ...
!> \param use_virial ...
!> \param my_virial_a ...
!> \param my_virial_b ...
!> \param hdab ...
!> \param a_hdab ...
!> \param use_subpatch ...
!> \param subpatch_pattern ...
! **************************************************************************************************
   SUBROUTINE integrate_pgf_product(la_max, zeta, la_min, &
                                    lb_max, zetb, lb_min, &
                                    ra, rab, rsgrid, cell, &
                                    cube_info, hab, pab, o1, o2, &
                                    radius, &
                                    calculate_forces, force_a, force_b, &
                                    compute_tau, &
                                    use_virial, my_virial_a, &
                                    my_virial_b, hdab, a_hdab, use_subpatch, subpatch_pattern)

      INTEGER, INTENT(IN)                                :: la_max
      REAL(KIND=dp), INTENT(IN)                          :: zeta
      INTEGER, INTENT(IN)                                :: la_min, lb_max
      REAL(KIND=dp), INTENT(IN)                          :: zetb
      INTEGER, INTENT(IN)                                :: lb_min
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), TARGET    :: ra, rab
      TYPE(realspace_grid_type), POINTER                 :: rsgrid
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cube_info_type), INTENT(IN)                   :: cube_info
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: pab
      INTEGER, INTENT(IN)                                :: o1, o2
      REAL(KIND=dp), INTENT(IN)                          :: radius
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(3), INTENT(INOUT), &
         OPTIONAL                                        :: force_a, force_b
      LOGICAL, INTENT(IN), OPTIONAL                      :: compute_tau, use_virial
      REAL(KIND=dp), DIMENSION(3, 3), OPTIONAL           :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: hdab
      REAL(KIND=dp), DIMENSION(:, :, :, :), OPTIONAL, &
         POINTER                                         :: a_hdab
      LOGICAL, OPTIONAL                                  :: use_subpatch
      INTEGER, INTENT(IN), OPTIONAL                      :: subpatch_pattern

      INTEGER                                            :: border_mask
      INTEGER, DIMENSION(3), TARGET                      :: border_width, npts_global, npts_local, &
                                                            shift_local
      LOGICAL                                            :: my_use_virial
      LOGICAL(KIND=C_BOOL)                               :: my_compute_tau, orthorhombic
      REAL(KIND=dp), DIMENSION(3, 2), TARGET             :: forces
      REAL(KIND=dp), DIMENSION(3, 3, 2), TARGET          :: virials
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: grid
      TYPE(C_PTR)                                        :: a_hdab_cptr, forces_cptr, hdab_cptr, &
                                                            pab_cptr, virials_cptr
      INTERFACE
         SUBROUTINE grid_ref_integrate_pgf_product_c(orthorhombic, compute_tau, &
                                                     border_mask, &
                                                     la_max, la_min, lb_max, lb_min, &
                                                     zeta, zetb, dh, dh_inv, ra, rab, &
                                                     npts_global, npts_local, shift_local, border_width, &
                                                     radius, o1, o2, n1, n2, grid, hab, pab, &
                                                     forces, virials, hdab, a_hdab) &
            BIND(C, name="grid_ref_integrate_pgf_product")
            IMPORT :: C_PTR, C_INT, C_DOUBLE, C_BOOL
            LOGICAL(KIND=C_BOOL), VALUE               :: orthorhombic
            LOGICAL(KIND=C_BOOL), VALUE               :: compute_tau
            INTEGER(KIND=C_INT), VALUE                :: border_mask
            INTEGER(KIND=C_INT), VALUE                :: la_max
            INTEGER(KIND=C_INT), VALUE                :: la_min
            INTEGER(KIND=C_INT), VALUE                :: lb_max
            INTEGER(KIND=C_INT), VALUE                :: lb_min
            REAL(KIND=C_DOUBLE), VALUE                :: zeta
            REAL(KIND=C_DOUBLE), VALUE                :: zetb
            TYPE(C_PTR), VALUE                        :: dh
            TYPE(C_PTR), VALUE                        :: dh_inv
            TYPE(C_PTR), VALUE                        :: ra
            TYPE(C_PTR), VALUE                        :: rab
            TYPE(C_PTR), VALUE                        :: npts_global
            TYPE(C_PTR), VALUE                        :: npts_local
            TYPE(C_PTR), VALUE                        :: shift_local
            TYPE(C_PTR), VALUE                        :: border_width
            REAL(KIND=C_DOUBLE), VALUE                :: radius
            INTEGER(KIND=C_INT), VALUE                :: o1
            INTEGER(KIND=C_INT), VALUE                :: o2
            INTEGER(KIND=C_INT), VALUE                :: n1
            INTEGER(KIND=C_INT), VALUE                :: n2
            TYPE(C_PTR), VALUE                        :: grid
            TYPE(C_PTR), VALUE                        :: hab
            TYPE(C_PTR), VALUE                        :: pab
            TYPE(C_PTR), VALUE                        :: forces
            TYPE(C_PTR), VALUE                        :: virials
            TYPE(C_PTR), VALUE                        :: hdab
            TYPE(C_PTR), VALUE                        :: a_hdab
         END SUBROUTINE grid_ref_integrate_pgf_product_c
      END INTERFACE

      MARK_USED(cube_info)
      MARK_USED(cell)

      IF (radius == 0.0_dp) THEN
         RETURN
      ENDIF

      border_mask = 0
      IF (PRESENT(use_subpatch)) THEN
         IF (use_subpatch) THEN
            CPASSERT(PRESENT(subpatch_pattern))
            border_mask = IAND(63, NOT(subpatch_pattern))  ! invert last 6 bits
         END IF
      END IF

      ! my_compute_tau defaults to .FALSE.
      ! IF (.true.) it will compute 0.5 * (nabla x_a).(v(r) nabla x_b)
      IF (PRESENT(compute_tau)) THEN
         my_compute_tau = LOGICAL(compute_tau, C_BOOL)
      ELSE
         my_compute_tau = .FALSE.
      ENDIF

      IF (PRESENT(use_virial)) THEN
         my_use_virial = use_virial
      ELSE
         my_use_virial = .FALSE.
      ENDIF

      IF (calculate_forces) THEN
         CPASSERT(PRESENT(pab))
         forces(:, :) = 0.0_dp
         forces_cptr = C_LOC(forces(1, 1))
         pab_cptr = C_LOC(pab(1, 1))
      ELSE
         forces_cptr = C_NULL_PTR
         pab_cptr = C_NULL_PTR
      ENDIF

      IF (my_use_virial) THEN
         CPASSERT(PRESENT(pab))
         virials(:, :, :) = 0.0_dp
         virials_cptr = C_LOC(virials(1, 1, 1))
      ELSE
         virials_cptr = C_NULL_PTR
      ENDIF

      IF (PRESENT(hdab) .AND. calculate_forces) THEN
         CPASSERT(PRESENT(pab))
         hdab(:, :, :) = 0.0_dp
         hdab_cptr = C_LOC(hdab(1, 1, 1))
      ELSE
         hdab_cptr = C_NULL_PTR
      ENDIF

      IF (PRESENT(a_hdab) .AND. my_use_virial) THEN
         CPASSERT(PRESENT(pab))
         a_hdab_cptr = C_LOC(a_hdab(1, 1, 1, 1))
      ELSE
         a_hdab_cptr = C_NULL_PTR
      ENDIF

      orthorhombic = LOGICAL(rsgrid%desc%orthorhombic, C_BOOL)

      CALL get_rsgrid_properties(rsgrid, npts_global=npts_global, &
                                 npts_local=npts_local, &
                                 shift_local=shift_local, &
                                 border_width=border_width)

      grid(1:, 1:, 1:) => rsgrid%r(:, :, :)  ! pointer assignment

#if __GNUC__ >= 9
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh))
      CPASSERT(IS_CONTIGUOUS(rsgrid%desc%dh_inv))
      CPASSERT(IS_CONTIGUOUS(ra))
      CPASSERT(IS_CONTIGUOUS(rab))
      CPASSERT(IS_CONTIGUOUS(npts_global))
      CPASSERT(IS_CONTIGUOUS(npts_local))
      CPASSERT(IS_CONTIGUOUS(shift_local))
      CPASSERT(IS_CONTIGUOUS(border_width))
      CPASSERT(IS_CONTIGUOUS(grid))
      CPASSERT(IS_CONTIGUOUS(hab))
      CPASSERT(IS_CONTIGUOUS(forces))
      CPASSERT(IS_CONTIGUOUS(virials))
      IF (PRESENT(pab)) THEN
         CPASSERT(IS_CONTIGUOUS(pab))
      ENDIF
      IF (PRESENT(hdab)) THEN
         CPASSERT(IS_CONTIGUOUS(hdab))
      ENDIF
      IF (PRESENT(a_hdab)) THEN
         CPASSERT(IS_CONTIGUOUS(a_hdab))
      ENDIF
#endif

      CALL grid_ref_integrate_pgf_product_c(orthorhombic=orthorhombic, &
                                            compute_tau=my_compute_tau, &
                                            border_mask=border_mask, &
                                            la_max=la_max, &
                                            la_min=la_min, &
                                            lb_max=lb_max, &
                                            lb_min=lb_min, &
                                            zeta=zeta, &
                                            zetb=zetb, &
                                            dh=C_LOC(rsgrid%desc%dh(1, 1)), &
                                            dh_inv=C_LOC(rsgrid%desc%dh_inv(1, 1)), &
                                            ra=C_LOC(ra(1)), &
                                            rab=C_LOC(rab(1)), &
                                            npts_global=C_LOC(npts_global(1)), &
                                            npts_local=C_LOC(npts_local(1)), &
                                            shift_local=C_LOC(shift_local(1)), &
                                            border_width=C_LOC(border_width(1)), &
                                            radius=radius, &
                                            o1=o1, &
                                            o2=o2, &
                                            n1=SIZE(hab, 1), &
                                            n2=SIZE(hab, 2), &
                                            grid=C_LOC(grid(1, 1, 1)), &
                                            hab=C_LOC(hab(1, 1)), &
                                            pab=pab_cptr, &
                                            forces=forces_cptr, &
                                            virials=virials_cptr, &
                                            hdab=hdab_cptr, &
                                            a_hdab=a_hdab_cptr)

      IF (calculate_forces .AND. PRESENT(force_a)) &
         force_a = force_a + forces(:, 1)
      IF (calculate_forces .AND. PRESENT(force_b)) &
         force_b = force_b + forces(:, 2)
      IF (my_use_virial .AND. PRESENT(my_virial_a)) &
         my_virial_a = my_virial_a + virials(:, :, 1)
      IF (my_use_virial .AND. PRESENT(my_virial_b)) &
         my_virial_b = my_virial_b + virials(:, :, 2)

   END SUBROUTINE integrate_pgf_product

END MODULE grid_integrate

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE grid_ibase_ref
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_add_block_node,&
                                              dbcsr_distribution_type,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_work_create
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE orbital_pointers,                ONLY: ncoset
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
   USE realspace_grid_types,            ONLY: realspace_grid_p_type
   USE task_list_methods,               ONLY: int2pair
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_ibase_ref'

   PUBLIC :: grid_integrate_ref

! **************************************************************************************************

CONTAINS

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
!> \param my_compute_tau ...
!> \param pab_required ...
!> \param map_consistent ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param posat ...
!> \param calculate_forces ...
!> \param force ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE grid_integrate_ref(rs_v, dhmat, deltap, task_list, cube_info, cell, &
                                 eps_gvg_rspace, admm_scal_fac, &
                                 my_compute_tau, pab_required, map_consistent, atom_of_kind, kind_of, posat, &
                                 calculate_forces, force, use_virial, virial)

      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dhmat, deltap
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: eps_gvg_rspace, admm_scal_fac
      LOGICAL, INTENT(IN)                                :: my_compute_tau, pab_required, &
                                                            map_consistent
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      LOGICAL, INTENT(IN)                                :: calculate_forces
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(len=*), PARAMETER :: routineN = 'grid_integrate_ref', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, atom_b, bcol, brow, iatom, igrid_level, ikind, ilevel, img, ipair, ipgf, &
         ipgf_new, iset, iset_new, iset_old, itask, jatom, jkind, jpgf, jpgf_new, jset, jset_new, &
         jset_old, maxco, maxpgf, maxset, maxsgf_set, na1, na2, natom, nb1, nb2, ncoa, ncob, &
         nimages, sgfa, sgfb
      INTEGER(kind=int_8), DIMENSION(:, :), POINTER      :: tasks
      INTEGER, DIMENSION(:), POINTER                     :: first_sgfa, first_sgfb, la_max, la_min, &
                                                            lb_max, lb_min, npgfa, npgfb, nsgfa, &
                                                            nsgfb
      LOGICAL                                            :: atom_pair_changed, atom_pair_done, &
                                                            found, new_set_pair_coming, &
                                                            use_subpatch
      REAL(KIND=dp)                                      :: dab, rab2, scalef, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra, rab, rab_inv, rb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dist_ab, h_block, hab, p_block, pab, &
                                                            sphi_a, sphi_b, work, zeta, zetb

!     TYPE(dbcsr_distribution_type)                      :: dist

      ! short cuts to task list variables
      tasks => task_list%tasks
      dist_ab => task_list%dist_ab
      maxco = task_list%maxco
      maxsgf_set = task_list%maxsgf_set
      maxset = task_list%maxset
      maxpgf = task_list%maxpgf
      nimages = task_list%nimages
      natom = task_list%natoms

      ! Allocate work storage
      ALLOCATE (work(maxco, maxsgf_set), hab(maxco, maxco))
      IF (pab_required) THEN
         ALLOCATE (pab(maxco, maxco))
      END IF

      iset_old = -1; jset_old = -1

      ! Here we loop over gridlevels first, finalising the matrix after each grid level is
      ! completed.  On each grid level, we loop over atom pairs, which will only access
      loop_gridlevels: DO igrid_level = 1, task_list%ngrid_levels

         DO img = 1, nimages
            CALL dbcsr_work_create(dhmat(img)%matrix, work_mutable=.TRUE.)
         END DO

         loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
         loop_tasks: DO itask = task_list%taskstart(ipair, igrid_level), task_list%taskstop(ipair, igrid_level)

            CALL int2pair(tasks(3, itask), ilevel, img, iatom, jatom, iset, jset, ipgf, jpgf, &
                          nimages, natom, maxset, maxpgf)

            ! At the start of a block of tasks, get atom data (and kind data, if needed)
            IF (itask .EQ. task_list%taskstart(ipair, igrid_level)) THEN

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

               ra(:) = posat(:, iatom)

               IF (iatom <= jatom) THEN
                  brow = iatom
                  bcol = jatom
               ELSE
                  brow = jatom
                  bcol = iatom
               END IF

               NULLIFY (h_block)
               CALL dbcsr_get_block_p(dhmat(img)%matrix, brow, bcol, h_block, found)
               IF (.NOT. ASSOCIATED(h_block)) THEN
                  CALL dbcsr_add_block_node(dhmat(img)%matrix, brow, bcol, h_block)
               END IF

               IF (pab_required) THEN
                  CALL dbcsr_get_block_p(matrix=deltap(img)%matrix, &
                                         row=brow, col=bcol, BLOCK=p_block, found=found)
                  CPASSERT(found)
               END IF

               IF (calculate_forces) THEN
                  atom_a = atom_of_kind(iatom)
                  atom_b = atom_of_kind(jatom)
                  force_a(:) = 0.0_dp
                  force_b(:) = 0.0_dp
               ENDIF
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  my_virial_b = 0.0_dp
               ENDIF

               atom_pair_changed = .TRUE.

            ELSE

               atom_pair_changed = .FALSE.

            ENDIF

            IF (atom_pair_changed .OR. iset_old .NE. iset .OR. jset_old .NE. jset) THEN

               ncoa = npgfa(iset)*ncoset(la_max(iset))
               sgfa = first_sgfa(iset)
               ncob = npgfb(jset)*ncoset(lb_max(jset))
               sgfb = first_sgfb(jset)
               IF (pab_required) THEN
                  IF (iatom <= jatom) THEN
                     CALL dgemm("N", "N", ncoa, nsgfb(jset), nsgfa(iset), &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                p_block(sgfa, sgfb), SIZE(p_block, 1), &
                                0.0_dp, work(1, 1), SIZE(work, 1))
                     CALL dgemm("N", "T", ncoa, ncob, nsgfb(jset), &
                                1.0_dp, work(1, 1), SIZE(work, 1), &
                                sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                                0.0_dp, pab(1, 1), SIZE(pab, 1))
                  ELSE
                     CALL dgemm("N", "N", ncob, nsgfa(iset), nsgfb(jset), &
                                1.0_dp, sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                                p_block(sgfb, sgfa), SIZE(p_block, 1), &
                                0.0_dp, work(1, 1), SIZE(work, 1))
                     CALL dgemm("N", "T", ncob, ncoa, nsgfa(iset), &
                                1.0_dp, work(1, 1), SIZE(work, 1), &
                                sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                0.0_dp, pab(1, 1), SIZE(pab, 1))
                  END IF
               END IF

               IF (iatom <= jatom) THEN
                  hab(1:ncoa, 1:ncob) = 0._dp
               ELSE
                  hab(1:ncob, 1:ncoa) = 0._dp
               ENDIF

               iset_old = iset
               jset_old = jset

            ENDIF

            rab(1) = dist_ab(1, itask)
            rab(2) = dist_ab(2, itask)
            rab(3) = dist_ab(3, itask)
            rab2 = DOT_PRODUCT(rab, rab)
            rb(1) = ra(1)+rab(1)
            rb(2) = ra(2)+rab(2)
            rb(3) = ra(3)+rab(3)
            zetp = zeta(ipgf, iset)+zetb(jpgf, jset)
            dab = SQRT(rab2)

            na1 = (ipgf-1)*ncoset(la_max(iset))+1
            na2 = ipgf*ncoset(la_max(iset))
            nb1 = (jpgf-1)*ncoset(lb_max(jset))+1
            nb2 = jpgf*ncoset(lb_max(jset))

            IF (tasks(4, itask) .EQ. 2) THEN
               use_subpatch = .TRUE.
            ELSE
               use_subpatch = .FALSE.
            ENDIF

            IF (pab_required) THEN
               IF (iatom <= jatom) THEN
                  CALL integrate_pgf_product_rspace( &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     ra, rab, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, pab=pab, o1=na1-1, o2=nb1-1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_a, force_b=force_b, &
                     compute_tau=my_compute_tau, map_consistent=map_consistent, &
                     use_virial=use_virial, my_virial_a=my_virial_a, &
                     my_virial_b=my_virial_b, use_subpatch=use_subpatch, subpatch_pattern=tasks(6, itask))
               ELSE
                  rab_inv = -rab
                  CALL integrate_pgf_product_rspace( &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     rb, rab_inv, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, pab=pab, o1=nb1-1, o2=na1-1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_b, force_b=force_a, &
                     compute_tau=my_compute_tau, map_consistent=map_consistent, &
                     use_virial=use_virial, my_virial_a=my_virial_b, &
                     my_virial_b=my_virial_a, use_subpatch=use_subpatch, subpatch_pattern=tasks(6, itask))
               END IF
            ELSE
               IF (iatom <= jatom) THEN
                  CALL integrate_pgf_product_rspace( &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     ra, rab, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, o1=na1-1, o2=nb1-1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_a, force_b=force_b, &
                     compute_tau=my_compute_tau, &
                     map_consistent=map_consistent, use_subpatch=use_subpatch, subpatch_pattern=tasks(6, itask))
               ELSE
                  rab_inv = -rab
                  CALL integrate_pgf_product_rspace( &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     rb, rab_inv, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, o1=nb1-1, o2=na1-1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_b, force_b=force_a, &
                     compute_tau=my_compute_tau, &
                     map_consistent=map_consistent, use_subpatch=use_subpatch, subpatch_pattern=tasks(6, itask))
               END IF
            END IF

            new_set_pair_coming = .FALSE.
            atom_pair_done = .FALSE.
            IF (itask < task_list%taskstop(ipair, igrid_level)) THEN
               CALL int2pair(tasks(3, itask+1), ilevel, img, iatom, jatom, iset_new, jset_new, ipgf_new, jpgf_new, &
                             nimages, natom, maxset, maxpgf)
               IF (iset_new .NE. iset .OR. jset_new .NE. jset) THEN
                  new_set_pair_coming = .TRUE.
               ENDIF
            ELSE
               ! do not forget the last block
               new_set_pair_coming = .TRUE.
               atom_pair_done = .TRUE.
            ENDIF

            ! contract the block into h if we're done with the current set pair
            IF (new_set_pair_coming) THEN
               IF (iatom <= jatom) THEN
                  CALL dgemm("N", "N", ncoa, nsgfb(jset), ncob, &
                             1.0_dp, hab(1, 1), SIZE(hab, 1), &
                             sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))
                  CALL dgemm("T", "N", nsgfa(iset), nsgfb(jset), ncoa, &
                             1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             work(1, 1), SIZE(work, 1), &
                             1.0_dp, h_block(sgfa, sgfb), SIZE(h_block, 1))
               ELSE
                  CALL dgemm("N", "N", ncob, nsgfa(iset), ncoa, &
                             1.0_dp, hab(1, 1), SIZE(hab, 1), &
                             sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                             0.0_dp, work(1, 1), SIZE(work, 1))
                  CALL dgemm("T", "N", nsgfb(jset), nsgfa(iset), ncob, &
                             1.0_dp, sphi_b(1, sgfb), SIZE(sphi_b, 1), &
                             work(1, 1), SIZE(work, 1), &
                             1.0_dp, h_block(sgfb, sgfa), SIZE(h_block, 1))
               END IF
            END IF

            IF (atom_pair_done) THEN
               IF (iatom == jatom) THEN
                  scalef = 1.0_dp
               ELSE
                  scalef = 2.0_dp
               END IF
               IF (calculate_forces) THEN
                  force(ikind)%rho_elec(:, atom_a) = &
                     force(ikind)%rho_elec(:, atom_a)+scalef*admm_scal_fac*force_a(:)
                  force(jkind)%rho_elec(:, atom_b) = &
                     force(jkind)%rho_elec(:, atom_b)+scalef*admm_scal_fac*force_b(:)
               ENDIF
               IF (use_virial) THEN
                  IF (use_virial .AND. calculate_forces) THEN
                     virial%pv_virial = virial%pv_virial+scalef*admm_scal_fac*my_virial_a
                     virial%pv_virial = virial%pv_virial+scalef*admm_scal_fac*my_virial_b
                  END IF
               END IF
            ENDIF
         END DO loop_tasks
         END DO loop_pairs

         DO img = 1, nimages
            CALL dbcsr_finalize(dhmat(img)%matrix)
         END DO

      END DO loop_gridlevels

      ! Deallocate work storage
      DEALLOCATE (work, hab)
      IF (pab_required) THEN
         DEALLOCATE (pab)
      END IF

   END SUBROUTINE grid_integrate_ref

END MODULE grid_ibase_ref

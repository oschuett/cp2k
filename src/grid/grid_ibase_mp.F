!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE grid_ibase_mp
   USE cell_types,                      ONLY: cell_type
   USE cube_utils,                      ONLY: cube_info_type
   USE dbcsr_api,                       ONLY: dbcsr_add_block_node,&
                                              dbcsr_distribution_get,&
                                              dbcsr_distribution_type,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_work_create
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: ncoset
   USE qs_integrate_potential_low,      ONLY: integrate_pgf_product_rspace
   USE realspace_grid_types,            ONLY: realspace_grid_p_type
   USE task_list_types,                 ONLY: task_list_type,&
                                              task_type
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_ibase_mp'

   PUBLIC :: grid_integrate_mp

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
!> \param kind_of ...
!> \param posat ...
!> \param calculate_forces ...
!> \param force ...
!> \param use_virial ...
!> \param virial ...
! **************************************************************************************************
   SUBROUTINE grid_integrate_mp(rs_v, dhmat, deltap, task_list, cube_info, cell, &
                                eps_gvg_rspace, admm_scal_fac, &
                                my_compute_tau, pab_required, map_consistent, kind_of, posat, &
                                calculate_forces, force, use_virial, virial)

      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dhmat, deltap
      TYPE(task_list_type), POINTER                      :: task_list
      TYPE(cube_info_type), DIMENSION(:), POINTER        :: cube_info
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(IN)                          :: eps_gvg_rspace, admm_scal_fac
      LOGICAL, INTENT(IN)                                :: my_compute_tau, pab_required, &
                                                            map_consistent
      INTEGER, DIMENSION(:), INTENT(IN)                  :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: posat
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: force
      LOGICAL, INTENT(IN)                                :: use_virial
      TYPE(virial_type), POINTER                         :: virial

      CHARACTER(len=*), PARAMETER :: routineN = 'grid_integrate_mp', &
         routineP = moduleN//':'//routineN

      INTEGER :: bcol, brow, iatom, igrid_level, ikind, ilevel, img, ipair, ipgf, ipgf_new, iset, &
         iset_new, iset_old, itask, ithread, jatom, jkind, jpgf, jpgf_new, jset, jset_new, &
         jset_old, maxco, maxpgf, maxset, maxsgf_set, na1, na2, natom, nb1, nb2, ncoa, ncob, &
         nimages, nthread, sgfa, sgfb
      INTEGER, DIMENSION(:), POINTER                     :: first_sgfa, first_sgfb, la_max, la_min, &
                                                            lb_max, lb_min, npgfa, npgfb, nsgfa, &
                                                            nsgfb
      LOGICAL                                            :: atom_pair_changed, atom_pair_done, &
                                                            found, has_threads, &
                                                            new_set_pair_coming, use_subpatch
      REAL(KIND=dp)                                      :: dab, rab2, scalef, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra, rab, rab_inv, rb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, hab, p_block, pab, sphi_a, &
                                                            sphi_b, work, zeta, zetb
      TYPE(dbcsr_distribution_type)                      :: dist
      TYPE(task_type), DIMENSION(:), POINTER             :: tasks

      ! short cuts to task list variables
      tasks => task_list%tasks
      maxco = task_list%maxco
      maxsgf_set = task_list%maxsgf_set
      maxset = task_list%maxset
      maxpgf = task_list%maxpgf
      nimages = task_list%nimages
      natom = task_list%natoms

      nthread = 1
!$    nthread = omp_get_max_threads()

!$OMP PARALLEL DEFAULT(NONE), &
!$OMP SHARED(tasks,natom,maxset), &
!$OMP SHARED(maxco,maxsgf_set,maxpgf,dhmat,deltap,use_virial,admm_scal_fac), &
!$OMP SHARED(pab_required,calculate_forces,ncoset,rs_v,cube_info,my_compute_tau), &
!$OMP SHARED(map_consistent,eps_gvg_rspace,force,virial,cell,posat,kind_of), &
!$OMP SHARED(task_list,nthread,nimages), &
!$OMP PRIVATE(ithread,work,hab,pab,iset_old,jset_old), &
!$OMP PRIVATE(iatom,jatom,iset,jset,ikind,jkind,ilevel,ipgf,jpgf), &
!$OMP PRIVATE(img,brow,bcol,first_sgfa,la_max,la_min,npgfa,nsgfa), &
!$OMP PRIVATE(sphi_a,zeta,first_sgfb,lb_max,lb_min,npgfb), &
!$OMP PRIVATE(nsgfb,sphi_b,zetb,found), &
!$OMP PRIVATE(force_a,force_b,my_virial_a,my_virial_b,atom_pair_changed,h_block), &
!$OMP PRIVATE(p_block,ncoa,sgfa,ncob,sgfb,rab,rab2,ra,rb,zetp,dab,igrid_level), &
!$OMP PRIVATE(na1,na2,nb1,nb2,use_subpatch,rab_inv,new_set_pair_coming,atom_pair_done), &
!$OMP PRIVATE(iset_new,jset_new,ipgf_new,jpgf_new,scalef), &
!$OMP PRIVATE(itask,dist,has_threads)

      ithread = 0
!$    ithread = omp_get_thread_num()

      ! Allocate work storage
      ALLOCATE (work(maxco, maxsgf_set), hab(maxco, maxco))
      IF (pab_required) THEN
         ALLOCATE (pab(maxco, maxco))
      END IF

      iset_old = -1; jset_old = -1

      ! Here we loop over gridlevels first, finalising the matrix after each grid level is
      ! completed.  On each grid level, we loop over atom pairs, which will only access
      ! a single block of each matrix, so with OpenMP, each matrix block is only touched
      ! by a single thread for each grid level
      loop_gridlevels: DO igrid_level = 1, task_list%ngrid_levels

         DO img = 1, nimages
            CALL dbcsr_work_create(dhmat(img)%matrix, work_mutable=.TRUE., n=nthread)
            CALL dbcsr_get_info(dhmat(img)%matrix, distribution=dist)
            CALL dbcsr_distribution_get(dist, has_threads=has_threads)
!$          IF (.NOT. has_threads) &
!$             CPABORT("No thread distribution defined.")
         END DO
!$OMP BARRIER

!$OMP DO schedule (dynamic, MAX(1,task_list%npairs(igrid_level)/(nthread*50)))
         loop_pairs: DO ipair = 1, task_list%npairs(igrid_level)
         loop_tasks: DO itask = task_list%taskstart(ipair, igrid_level), task_list%taskstop(ipair, igrid_level)

            ilevel = tasks(itask)%grid_level
            img = tasks(itask)%image
            iatom = tasks(itask)%iatom
            jatom = tasks(itask)%jatom
            iset = tasks(itask)%iset
            jset = tasks(itask)%jset
            ipgf = tasks(itask)%ipgf
            jpgf = tasks(itask)%jpgf
            rab = tasks(itask)%rab

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

               force_a(:) = 0.0_dp
               force_b(:) = 0.0_dp
               my_virial_a = 0.0_dp
               my_virial_b = 0.0_dp
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

            rab2 = DOT_PRODUCT(rab, rab)
            rb(1) = ra(1) + rab(1)
            rb(2) = ra(2) + rab(2)
            rb(3) = ra(3) + rab(3)
            zetp = zeta(ipgf, iset) + zetb(jpgf, jset)
            dab = SQRT(rab2)

            na1 = (ipgf - 1)*ncoset(la_max(iset)) + 1
            na2 = ipgf*ncoset(la_max(iset))
            nb1 = (jpgf - 1)*ncoset(lb_max(jset)) + 1
            nb2 = jpgf*ncoset(lb_max(jset))

            ! check whether we need to use fawzi's generalised collocation scheme
            IF (rs_v(igrid_level)%rs_grid%desc%distributed) THEN
               !tasks(4,:) is 0 for replicated, 1 for distributed 2 for exceptional distributed tasks
               IF (tasks(itask)%dist_type .EQ. 2) THEN
                  use_subpatch = .TRUE.
               ELSE
                  use_subpatch = .FALSE.
               ENDIF
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
                     hab, pab=pab, o1=na1 - 1, o2=nb1 - 1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_a, force_b=force_b, &
                     compute_tau=my_compute_tau, map_consistent=map_consistent, &
                     use_virial=use_virial, my_virial_a=my_virial_a, &
                     my_virial_b=my_virial_b, use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern)
               ELSE
                  rab_inv = -rab
                  CALL integrate_pgf_product_rspace( &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     rb, rab_inv, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, pab=pab, o1=nb1 - 1, o2=na1 - 1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_b, force_b=force_a, &
                     compute_tau=my_compute_tau, map_consistent=map_consistent, &
                     use_virial=use_virial, my_virial_a=my_virial_b, &
                     my_virial_b=my_virial_a, use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern)
               END IF
            ELSE
               IF (iatom <= jatom) THEN
                  CALL integrate_pgf_product_rspace( &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     ra, rab, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, o1=na1 - 1, o2=nb1 - 1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_a, force_b=force_b, &
                     compute_tau=my_compute_tau, &
                     map_consistent=map_consistent, use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern)
               ELSE
                  rab_inv = -rab
                  CALL integrate_pgf_product_rspace( &
                     lb_max(jset), zetb(jpgf, jset), lb_min(jset), &
                     la_max(iset), zeta(ipgf, iset), la_min(iset), &
                     rb, rab_inv, rab2, rs_v(igrid_level)%rs_grid, cell, &
                     cube_info(igrid_level), &
                     hab, o1=nb1 - 1, o2=na1 - 1, &
                     eps_gvg_rspace=eps_gvg_rspace, &
                     calculate_forces=calculate_forces, &
                     force_a=force_b, force_b=force_a, &
                     compute_tau=my_compute_tau, &
                     map_consistent=map_consistent, use_subpatch=use_subpatch, subpatch_pattern=tasks(itask)%subpatch_pattern)
               END IF
            END IF

            new_set_pair_coming = .FALSE.
            atom_pair_done = .FALSE.
            IF (itask < task_list%taskstop(ipair, igrid_level)) THEN
               ilevel = tasks(itask + 1)%grid_level
               img = tasks(itask + 1)%image
               iatom = tasks(itask + 1)%iatom
               jatom = tasks(itask + 1)%jatom
               iset_new = tasks(itask + 1)%iset
               jset_new = tasks(itask + 1)%jset
               ipgf_new = tasks(itask + 1)%ipgf
               jpgf_new = tasks(itask + 1)%jpgf
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
!$OMP CRITICAL(force_critical)
               IF (iatom == jatom) THEN
                  scalef = 1.0_dp
               ELSE
                  scalef = 2.0_dp
               END IF
               IF (calculate_forces) THEN
                  force(:, iatom) = force(:, iatom) + scalef*admm_scal_fac*force_a(:)
                  force(:, jatom) = force(:, jatom) + scalef*admm_scal_fac*force_b(:)
                  IF (use_virial) THEN
                     virial%pv_virial = virial%pv_virial + scalef*admm_scal_fac*my_virial_a
                     virial%pv_virial = virial%pv_virial + scalef*admm_scal_fac*my_virial_b
                  END IF
               END IF
!$OMP END CRITICAL (force_critical)
            ENDIF
         END DO loop_tasks
         END DO loop_pairs
!$OMP END DO

         DO img = 1, nimages
            CALL dbcsr_finalize(dhmat(img)%matrix)
         END DO

      END DO loop_gridlevels

      ! Deallocate work storage
      DEALLOCATE (work, hab)
      IF (pab_required) THEN
         DEALLOCATE (pab)
      END IF

!$OMP END PARALLEL

   END SUBROUTINE grid_integrate_mp

END MODULE grid_ibase_mp

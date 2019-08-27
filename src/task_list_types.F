!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types for task lists
!> \par History
!>      01.2008 [Joost VandeVondele] refactered out of qs_collocate / qs_integrate
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE task_list_types

   USE kinds,                           ONLY: dp,&
                                              int_8
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'task_list_types'

! **************************************************************************************************
   TYPE basis_type
      INTEGER, DIMENSION(:), POINTER                :: first_sgf => Null()
      INTEGER, DIMENSION(:), POINTER                :: lmax => Null()
      INTEGER, DIMENSION(:), POINTER                :: lmin => Null()
      INTEGER, DIMENSION(:), POINTER                :: npgf => Null()
      INTEGER, DIMENSION(:), POINTER                :: nsgf => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER       :: sphi => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER       :: zet => Null()
   END TYPE basis_type

   TYPE rgrid_type
      REAL(KIND=dp), DIMENSION(3, 3)                :: dh = 0.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                :: dh_inv = 0.0_dp
      INTEGER, DIMENSION(3)                         :: perd = 0
      INTEGER, DIMENSION(3)                         :: npts = 0
      INTEGER, DIMENSION(3)                         :: lb_grid = 0
      REAL(KIND=dp)                                 :: drmin = 0.0_dp
      INTEGER                                       :: max_radius = 0
      REAL(KIND=dp)                                 :: max_rad_ga = 0.0_dp
      INTEGER, DIMENSION(:, :), POINTER             :: lb_cube => Null()
      INTEGER, DIMENSION(:, :), POINTER             :: ub_cube => Null()
      INTEGER, DIMENSION(:, :), POINTER             :: sphere_bounds => Null()
   END TYPE rgrid_type

   TYPE task_list_type
      INTEGER(kind=int_8), DIMENSION(:, :), POINTER :: tasks => Null()
      REAL(KIND=dp), DIMENSION(:, :), POINTER       :: dist_ab => Null()
      INTEGER(kind=int_8), DIMENSION(:), POINTER    :: atom_pair_send => Null()
      INTEGER(kind=int_8), DIMENSION(:), POINTER    :: atom_pair_recv => Null()
      INTEGER                                       :: ntasks = 0
      INTEGER, DIMENSION(:, :), POINTER             :: taskstart => Null()
      INTEGER, DIMENSION(:, :), POINTER             :: taskstop => Null()
      INTEGER, DIMENSION(:), POINTER                :: npairs => Null()
      INTEGER                                       :: ngrid_levels = 0
      INTEGER                                       :: nimages = 0
      INTEGER                                       :: natoms = 0
      INTEGER                                       :: maxset = 0
      INTEGER                                       :: maxpgf = 0
      INTEGER                                       :: lmax_all = 0
      INTEGER                                       :: maxco = 0
      INTEGER                                       :: maxsgf_set = 0
      REAL(KIND=dp)                                 :: eps_rho_rspace = 0.0_dp
      TYPE(basis_type), DIMENSION(:), POINTER       :: basis => Null()
      TYPE(rgrid_type), DIMENSION(:), POINTER       :: cinfo => Null()
   END TYPE task_list_type

! **************************************************************************************************

   PUBLIC :: task_list_type, allocate_task_list, &
             deallocate_task_list, clear_task_list

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief allocates and initialised the components of the task_list_type
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE allocate_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      ALLOCATE (task_list)

   END SUBROUTINE allocate_task_list

! **************************************************************************************************
!> \brief deallocates the components and the object itself
!> \param task_list ...
!> \par History
!>      01.2008 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE deallocate_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      CALL clear_task_list(task_list)

      DEALLOCATE (task_list)

   END SUBROUTINE deallocate_task_list

! **************************************************************************************************
!> \brief deallocates the components, but not the object itself
!> \param task_list ...
!> \par History
!>      08.2019 factored out of deallocate_task_list [Ole Schuett]
! **************************************************************************************************
   SUBROUTINE clear_task_list(task_list)
      TYPE(task_list_type), POINTER                      :: task_list

      INTEGER                                            :: ib

      IF (ASSOCIATED(task_list%tasks)) &
         DEALLOCATE (task_list%tasks)
      IF (ASSOCIATED(task_list%dist_ab)) &
         DEALLOCATE (task_list%dist_ab)
      IF (ASSOCIATED(task_list%atom_pair_send)) &
         DEALLOCATE (task_list%atom_pair_send)
      IF (ASSOCIATED(task_list%atom_pair_recv)) &
         DEALLOCATE (task_list%atom_pair_recv)
      IF (ASSOCIATED(task_list%taskstart)) &
         DEALLOCATE (task_list%taskstart)
      IF (ASSOCIATED(task_list%taskstop)) &
         DEALLOCATE (task_list%taskstop)
      IF (ASSOCIATED(task_list%npairs)) &
         DEALLOCATE (task_list%npairs)

      IF (ASSOCIATED(task_list%basis)) THEN
         DO ib = 1, SIZE(task_list%basis)
            IF (ASSOCIATED(task_list%basis(ib)%first_sgf)) &
               DEALLOCATE (task_list%basis(ib)%first_sgf)
            IF (ASSOCIATED(task_list%basis(ib)%lmax)) &
               DEALLOCATE (task_list%basis(ib)%lmax)
            IF (ASSOCIATED(task_list%basis(ib)%lmin)) &
               DEALLOCATE (task_list%basis(ib)%lmin)
            IF (ASSOCIATED(task_list%basis(ib)%npgf)) &
               DEALLOCATE (task_list%basis(ib)%npgf)
            IF (ASSOCIATED(task_list%basis(ib)%nsgf)) &
               DEALLOCATE (task_list%basis(ib)%nsgf)
            IF (ASSOCIATED(task_list%basis(ib)%sphi)) &
               DEALLOCATE (task_list%basis(ib)%sphi)
            IF (ASSOCIATED(task_list%basis(ib)%zet)) &
               DEALLOCATE (task_list%basis(ib)%zet)
         END DO
         DEALLOCATE (task_list%basis)
      ENDIF

      IF (ASSOCIATED(task_list%cinfo)) THEN
         DO ib = 1, SIZE(task_list%cinfo)
            IF (ASSOCIATED(task_list%cinfo(ib)%lb_cube)) &
               DEALLOCATE (task_list%cinfo(ib)%lb_cube)
            IF (ASSOCIATED(task_list%cinfo(ib)%ub_cube)) &
               DEALLOCATE (task_list%cinfo(ib)%ub_cube)
            IF (ASSOCIATED(task_list%cinfo(ib)%sphere_bounds)) &
               DEALLOCATE (task_list%cinfo(ib)%sphere_bounds)
         END DO
         DEALLOCATE (task_list%cinfo)
      ENDIF

   END SUBROUTINE clear_task_list

END MODULE task_list_types

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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

!>        tasks(1,i) := destination
!>        tasks(2,i) := source
!>        tasks(3,i) := compressed type (iatom, jatom, ....)
!>        tasks(4,i) := type (0: replicated, 1: distributed local, 2: distributed generalised)
!>        tasks(5,i) := cost
!>        tasks(6,i) := alternate destination code (0 if none available)

   TYPE task_type
      INTEGER(kind=int_8)                           :: destination = 0 ! index 1
      INTEGER(kind=int_8)                           :: source = 0 ! index 2
      INTEGER(kind=int_8)                           :: dist_type = 0 ! index 4
      INTEGER(kind=int_8)                           :: cost = 0 ! index 5
      INTEGER(kind=int_8)                           :: subpatch_pattern = 0 ! index 6

      INTEGER                                       :: grid_level = 0
      INTEGER                                       :: image = 0
      INTEGER                                       :: iatom = 0
      INTEGER                                       :: jatom = 0
      INTEGER                                       :: iset = 0
      INTEGER                                       :: jset = 0
      INTEGER                                       :: ipgf = 0
      INTEGER                                       :: jpgf = 0

      REAL(KIND=dp), DIMENSION(3)                   :: rab = 0.0_dp
      REAL(KIND=dp)                                 :: radius = 0.0_dp
   END TYPE task_type

   TYPE task_list_type
      TYPE(task_type), DIMENSION(:), POINTER        :: tasks => Null()
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
      REAL(KIND=dp)                                 :: eps_rho_rspace
      TYPE(basis_type), DIMENSION(:), POINTER       :: basis => Null()
      TYPE(rgrid_type), DIMENSION(:), POINTER       :: cinfo => Null()
   END TYPE task_list_type

   TYPE grid_base_type
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER       :: grid => NULL()
   END TYPE grid_base_type

   INTEGER, PARAMETER                               :: task_size_in_int8 = 17

! **************************************************************************************************

   PUBLIC :: task_list_type, task_type, &
             allocate_task_list, deallocate_task_list, clear_task_list, reallocate_tasks, &
             grid_base_type, &
             serialize_task, deserialize_task, task_size_in_int8

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

! **************************************************************************************************
!> \brief Grow an array of tasks while preserving the existing entries.
!> \param tasks ...
!> \param new_size ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE reallocate_tasks(tasks, new_size)
      TYPE(task_type), DIMENSION(:), INTENT(inout), &
         POINTER                                         :: tasks
      INTEGER, INTENT(in)                                :: new_size

      TYPE(task_type), DIMENSION(:), POINTER             :: old_tasks

      IF (.NOT. ASSOCIATED(tasks)) THEN
         ALLOCATE (tasks(new_size))
         RETURN
      ENDIF

      IF (new_size <= SIZE(tasks)) RETURN

      old_tasks => tasks
      NULLIFY (tasks)

      ALLOCATE (tasks(new_size))
      tasks(1:SIZE(old_tasks)) = old_tasks(:)
      DEALLOCATE (old_tasks)
   END SUBROUTINE reallocate_tasks

! **************************************************************************************************
!> \brief ...
!> \param task ...
!> \param serialized_task ...
! **************************************************************************************************
   SUBROUTINE serialize_task(task, serialized_task)
      TYPE(task_type), INTENT(IN)                        :: task
      INTEGER(KIND=int_8), &
         DIMENSION(task_size_in_int8), INTENT(OUT)       :: serialized_task

      serialized_task(1) = task%destination
      serialized_task(2) = task%source
      serialized_task(3) = task%dist_type
      serialized_task(4) = task%cost
      serialized_task(5) = task%subpatch_pattern

      serialized_task(6) = task%grid_level
      serialized_task(7) = task%image
      serialized_task(8) = task%iatom
      serialized_task(9) = task%jatom
      serialized_task(10) = task%iset
      serialized_task(11) = task%jset
      serialized_task(12) = task%ipgf
      serialized_task(13) = task%jpgf

      serialized_task(14) = TRANSFER(task%rab(1), mold=1_int_8)
      serialized_task(15) = TRANSFER(task%rab(2), mold=1_int_8)
      serialized_task(16) = TRANSFER(task%rab(3), mold=1_int_8)
      serialized_task(17) = TRANSFER(task%radius, mold=1_int_8)
   END SUBROUTINE serialize_task

! **************************************************************************************************
!> \brief ...
!> \param task ...
!> \param serialized_task ...
! **************************************************************************************************
   SUBROUTINE deserialize_task(task, serialized_task)
      TYPE(task_type), INTENT(OUT)                       :: task
      INTEGER(KIND=int_8), &
         DIMENSION(task_size_in_int8), INTENT(IN)        :: serialized_task

      task%destination = serialized_task(1)
      task%source = serialized_task(2)
      task%dist_type = serialized_task(3)
      task%cost = serialized_task(4)
      task%subpatch_pattern = serialized_task(5)

      task%grid_level = INT(serialized_task(6))
      task%image = INT(serialized_task(7))
      task%iatom = INT(serialized_task(8))
      task%jatom = INT(serialized_task(9))
      task%iset = INT(serialized_task(10))
      task%jset = INT(serialized_task(11))
      task%ipgf = INT(serialized_task(12))
      task%jpgf = INT(serialized_task(13))

      task%rab(1) = TRANSFER(serialized_task(14), mold=1.0_dp)
      task%rab(2) = TRANSFER(serialized_task(15), mold=1.0_dp)
      task%rab(3) = TRANSFER(serialized_task(16), mold=1.0_dp)
      task%radius = TRANSFER(serialized_task(17), mold=1.0_dp)
   END SUBROUTINE deserialize_task

END MODULE task_list_types

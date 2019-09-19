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

!>        tasks(1,i) := destination
!>        tasks(2,i) := source
!>        tasks(3,i) := compressed type (iatom, jatom, ....)
!>        tasks(4,i) := type (0: replicated, 1: distributed local, 2: distributed generalised)
!>        tasks(5,i) := cost
!>        tasks(6,i) := alternate destination code (0 if none available)

   TYPE task_type
      INTEGER(kind=int_8)                           :: destination = 0 ! index 1
      INTEGER(kind=int_8)                           :: source = 0 ! index 2
      INTEGER(kind=int_8)                           :: pair = 0 ! index 3
      INTEGER(kind=int_8)                           :: dist_type = 0 ! index 4
      INTEGER(kind=int_8)                           :: cost = 0 ! index 5
      INTEGER(kind=int_8)                           :: subpatch_pattern = 0 ! index 6
   END TYPE task_type

   TYPE task_list_type
      TYPE(task_type), DIMENSION(:), POINTER        :: tasks => Null()
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
      REAL(KIND=dp)                                 :: eps_rho_rspace
      TYPE(basis_type), DIMENSION(:), POINTER       :: basis => Null()
      TYPE(rgrid_type), DIMENSION(:), POINTER       :: cinfo => Null()
   END TYPE task_list_type

   TYPE grid_base_type
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER       :: grid => NULL()
   END TYPE grid_base_type

! **************************************************************************************************

   PUBLIC :: task_list_type, task_type, &
             allocate_task_list, deallocate_task_list, clear_task_list, reallocate_tasks, &
             grid_base_type, pair2int, int2pair

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
!> \brief converts a pgf index pair (ipgf,iset,iatom),(jpgf,jset,jatom)
!>      to a unique integer.
!>      a list of integers can be sorted, and will result in a list of pgf pairs
!>      for which all atom pairs are grouped, and for each atom pair all set pairs are grouped
!>      and for each set pair, all pgfs are grouped
!> \param res ...
!> \param ilevel ...
!> \param image ...
!> \param iatom ...
!> \param jatom ...
!> \param iset ...
!> \param jset ...
!> \param ipgf ...
!> \param jpgf ...
!> \param nimages ...
!> \param natom ...
!> \param maxset ...
!> \param maxpgf ...
!> \par History
!>      11.2007 created [Joost]
!> \note
!>      will hopefully not overflow any more
! **************************************************************************************************
   SUBROUTINE pair2int(res, ilevel, image, iatom, jatom, iset, jset, ipgf, jpgf, &
                       nimages, natom, maxset, maxpgf)
      INTEGER(KIND=int_8), INTENT(OUT)                   :: res
      INTEGER, INTENT(IN)                                :: ilevel, image, iatom, jatom, iset, jset, &
                                                            ipgf, jpgf, nimages, natom, maxset, &
                                                            maxpgf

      INTEGER(KIND=int_8)                                :: maxpgf8, maxset8, natom8, nimages8, &
                                                            nlev1, nlev2, nlev3, nlev4

      natom8 = natom; maxset8 = maxset; maxpgf8 = maxpgf; nimages8 = nimages
      !
      ! this encoding yields the right order of the tasks for collocation after the sort
      ! in distribute_tasks. E.g. for a atom pair, all sets and pgfs are computed in one go.
      ! The exception is the gridlevel. Tasks are first ordered wrt to grid_level. This implies
      ! that a given density matrix block will be decontracted several times, but cache effects on the
      ! grid make up for this.
      !
      nlev1 = maxpgf8**2
      nlev2 = maxset8**2*nlev1
      nlev3 = natom8**2*nlev2
      nlev4 = nimages8*nlev3
      !
      res = ilevel*nlev4+ &
            (image-1)*nlev3+ &
            ((iatom-1)*natom8+jatom-1)*nlev2+ &
            ((iset-1)*maxset8+jset-1)*nlev1+ &
            (ipgf-1)*maxpgf8+jpgf-1

   END SUBROUTINE pair2int

! **************************************************************************************************
!> \brief ...
!> \param res ...
!> \param ilevel ...
!> \param image ...
!> \param iatom ...
!> \param jatom ...
!> \param iset ...
!> \param jset ...
!> \param ipgf ...
!> \param jpgf ...
!> \param nimages ...
!> \param natom ...
!> \param maxset ...
!> \param maxpgf ...
! **************************************************************************************************
   SUBROUTINE int2pair(res, ilevel, image, iatom, jatom, iset, jset, ipgf, jpgf, &
                       nimages, natom, maxset, maxpgf)
      INTEGER(KIND=int_8), INTENT(IN)                    :: res
      INTEGER, INTENT(OUT)                               :: ilevel, image, iatom, jatom, iset, jset, &
                                                            ipgf, jpgf
      INTEGER, INTENT(IN)                                :: nimages, natom, maxset, maxpgf

      INTEGER(KIND=int_8) :: iatom8, ijatom, ijset, img, ipgf8, iset8, jatom8, jpgf8, jset8, &
         maxpgf8, maxset8, natom8, nimages8, nlev1, nlev2, nlev3, nlev4, tmp

      natom8 = natom; maxset8 = maxset; maxpgf8 = maxpgf; nimages8 = nimages
      !
      nlev1 = maxpgf8**2
      nlev2 = maxset8**2*nlev1
      nlev3 = natom8**2*nlev2
      nlev4 = nimages8*nlev3
      !
      ilevel = INT(res/nlev4)
      tmp = MOD(res, nlev4)
      img = tmp/nlev3+1
      tmp = MOD(tmp, nlev3)
      ijatom = tmp/nlev2
      iatom8 = ijatom/natom8+1
      jatom8 = MOD(ijatom, natom8)+1
      tmp = MOD(tmp, nlev2)
      ijset = tmp/nlev1
      iset8 = ijset/maxset8+1
      jset8 = MOD(ijset, maxset8)+1
      tmp = MOD(tmp, nlev1)
      ipgf8 = tmp/maxpgf8+1
      jpgf8 = MOD(tmp, maxpgf8)+1
      !
      image = INT(img)
      iatom = INT(iatom8); jatom = INT(jatom8); iset = INT(iset8); jset = INT(jset8)
      ipgf = INT(ipgf8); jpgf = INT(jpgf8)

   END SUBROUTINE int2pair

END MODULE task_list_types

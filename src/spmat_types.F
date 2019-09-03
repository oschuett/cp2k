!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief A minimal sparse block matrix type (non-parallel)
!> \par
!> \author JGH
! **************************************************************************************************
MODULE spmat_types

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'spmat_types'

! **************************************************************************************************
   TYPE block_type
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE         :: mat
   END TYPE block_type
   TYPE index_type
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: col
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: blk
   END TYPE index_type

   TYPE spmat_type
      ! total number of rows and columns
      INTEGER                                            :: nrow = 0
      INTEGER                                            :: ncol = 0
      ! total number of non-zero blocks stored
      INTEGER                                            :: nnz = 0
      ! list of column and block indices for a given row
      TYPE(index_type), DIMENSION(:), ALLOCATABLE        :: rind
      ! list of block data
      TYPE(block_type), DIMENSION(:), ALLOCATABLE        :: spdata
   END TYPE spmat_type

   TYPE spmat_p_type
      TYPE(spmat_type)                                   :: spmat
   END TYPE spmat_p_type
! **************************************************************************************************
   ! chunk of blocks allocated
   INTEGER, PARAMETER                                    :: tblk = 100
! **************************************************************************************************

   PUBLIC :: spmat_p_type, spmat_type, spmat_allocate, spmat_release
   PUBLIC :: spmat_get, spmat_put, spmat_add

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param spmat ...
!> \param nrow ...
!> \param ncol ...
! **************************************************************************************************
   SUBROUTINE spmat_allocate(spmat, nrow, ncol)
      TYPE(spmat_type), INTENT(INOUT)                    :: spmat
      INTEGER, INTENT(IN)                                :: nrow, ncol

      spmat%nrow = nrow
      spmat%ncol = ncol

      spmat%nnz = 0
      ALLOCATE (spmat%rind(nrow))
      ALLOCATE (spmat%spdata(tblk))

   END SUBROUTINE spmat_allocate

! **************************************************************************************************
!> \brief ...
!> \param spmat ...
! **************************************************************************************************
   SUBROUTINE spmat_release(spmat)
      TYPE(spmat_type), INTENT(INOUT)                    :: spmat

      INTEGER                                            :: i

      spmat%nrow = 0
      spmat%ncol = 0

      spmat%nnz = 0
      IF (ALLOCATED(spmat%rind)) THEN
         DO i = 1, SIZE(spmat%rind)
            IF (ALLOCATED(spmat%rind(i)%col)) DEALLOCATE (spmat%rind(i)%col)
            IF (ALLOCATED(spmat%rind(i)%blk)) DEALLOCATE (spmat%rind(i)%blk)
         END DO
         DEALLOCATE (spmat%rind)
      END IF
      IF (ALLOCATED(spmat%spdata)) THEN
         DO i = 1, SIZE(spmat%spdata)
            IF (ALLOCATED(spmat%spdata(i)%mat)) DEALLOCATE (spmat%spdata(i)%mat)
         END DO
         DEALLOCATE (spmat%spdata)
      END IF

   END SUBROUTINE spmat_release

! **************************************************************************************************
!> \brief ...
!> \param spmat ...
!> \param irow ...
!> \param icol ...
!> \param blk ...
! **************************************************************************************************
   SUBROUTINE spmat_get(spmat, irow, icol, blk)
      TYPE(spmat_type), INTENT(IN), TARGET               :: spmat
      INTEGER, INTENT(IN)                                :: irow, icol
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: blk

      INTEGER                                            :: iblk, iloc

      CPASSERT(irow <= spmat%nrow .AND. icol <= spmat%ncol)
      CPASSERT(irow > 0 .AND. icol > 0)

      iloc = findpos(spmat%rind(irow)%col, icol)
      IF (iloc /= 0) THEN
         iblk = spmat%rind(irow)%blk(iloc)
         blk => spmat%spdata(iblk)%mat
      ELSE
         NULLIFY (blk)
      END IF

   END SUBROUTINE spmat_get

! **************************************************************************************************
!> \brief ...
!> \param spmat ...
!> \param irow ...
!> \param icol ...
!> \param blk ...
! **************************************************************************************************
   SUBROUTINE spmat_put(spmat, irow, icol, blk)
      TYPE(spmat_type), INTENT(INOUT), TARGET            :: spmat
      INTEGER, INTENT(IN)                                :: irow, icol
      REAL(KIND=dp), DIMENSION(:, :), INTENT(in)         :: blk

      INTEGER                                            :: i, iblk, iloc, n1, n2, nnz, nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blkn, coln
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: mat
      TYPE(block_type), ALLOCATABLE, DIMENSION(:)        :: spnew

      CPASSERT(irow <= spmat%nrow .AND. icol <= spmat%ncol)
      CPASSERT(irow > 0 .AND. icol > 0)
      iloc = findpos(spmat%rind(irow)%col, icol)
      IF (iloc /= 0) THEN
         iblk = spmat%rind(irow)%blk(iloc)
         mat => spmat%spdata(iblk)%mat
         CPASSERT(SIZE(mat, 1) == SIZE(blk, 1))
         CPASSERT(SIZE(mat, 2) == SIZE(blk, 2))
         mat = blk
      ELSE
         spmat%nnz = spmat%nnz+1
         IF (ALLOCATED(spmat%rind(irow)%col)) THEN
            nr = SIZE(spmat%rind(irow)%col)+1
         ELSE
            nr = 1
         END IF
         ALLOCATE (coln(nr), blkn(nr))
         IF (nr > 1) THEN
            coln(1:nr-1) = spmat%rind(irow)%col(:)
            blkn(1:nr-1) = spmat%rind(irow)%blk(:)
         END IF
         nnz = spmat%nnz
         coln(nr) = icol
         blkn(nr) = nnz
         CALL move_alloc(coln, spmat%rind(irow)%col)
         CALL move_alloc(blkn, spmat%rind(irow)%blk)
         n1 = SIZE(blk, 1)
         n2 = SIZE(blk, 2)
         IF (nnz > SIZE(spmat%spdata)) THEN
            ALLOCATE (spnew(nnz+tblk))
            DO i = 1, nnz-1
               CALL move_alloc(spmat%spdata(i)%mat, spnew(i)%mat)
            END DO
            ALLOCATE (spnew(nnz)%mat(n1, n2))
            spnew(nnz)%mat(1:n1, 1:n2) = blk(1:n1, 1:n2)
            CALL move_alloc(spnew, spmat%spdata)
         ELSE
            ALLOCATE (spmat%spdata(nnz)%mat(n1, n2))
            spmat%spdata(nnz)%mat(1:n1, 1:n2) = blk(1:n1, 1:n2)
         END IF
      END IF

   END SUBROUTINE spmat_put

   ! update matrix with: sp = alpha*sp + beta*blk
! **************************************************************************************************
!> \brief ...
!> \param spmat ...
!> \param irow ...
!> \param icol ...
!> \param blk ...
!> \param alpha ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE spmat_add(spmat, irow, icol, blk, alpha, beta)
      TYPE(spmat_type), INTENT(INOUT), TARGET            :: spmat
      INTEGER, INTENT(IN)                                :: irow, icol
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: blk
      REAL(KIND=dp), INTENT(IN)                          :: alpha, beta

      INTEGER                                            :: iblk, iloc
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: mat

      CPASSERT(irow <= spmat%nrow .AND. icol <= spmat%ncol)
      CPASSERT(irow > 0 .AND. icol > 0)
      iloc = findpos(spmat%rind(irow)%col, icol)
      CPASSERT(iloc /= 0)
      iblk = spmat%rind(irow)%blk(iloc)
      mat => spmat%spdata(iblk)%mat
      CPASSERT(SIZE(mat, 1) == SIZE(blk, 1))
      CPASSERT(SIZE(mat, 2) == SIZE(blk, 2))
      IF (alpha == 0.0_dp) THEN
         mat = beta*blk
      ELSE
         mat = alpha*mat+beta*blk
      END IF

   END SUBROUTINE spmat_add

! **************************************************************************************************
!> \brief ...
!> \param array ...
!> \param ival ...
!> \return ...
! **************************************************************************************************
   FUNCTION findpos(array, ival) RESULT(ipos)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(IN)     :: array
      INTEGER, INTENT(IN)                                :: ival
      INTEGER                                            :: ipos

      INTEGER                                            :: i

      ipos = 0
      IF (ALLOCATED(array)) THEN
         DO i = 1, SIZE(array)
            IF (ival == array(i)) THEN
               ipos = i
               EXIT
            END IF
         END DO
      END IF

   END FUNCTION findpos

END MODULE spmat_types

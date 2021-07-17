!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE block_tas_util

   !! often used utilities for tall-and-skinny matrices

   USE kinds, ONLY: int_4, int_8
   USE util, ONLY: sort

#include "../../base/base_uses.f90"

! #if TO_VERSION(1, 11) <= TO_VERSION(LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
!    USE libxsmm, ONLY: libxsmm_diff
! #  define PURE_ARRAY_EQ
! #else
#  define PURE_ARRAY_EQ PURE
! #endif

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_tas_util'

   PUBLIC :: &
      array_eq, &
      index_unique, &
      swap

   INTERFACE swap
      MODULE PROCEDURE swap_i8
      MODULE PROCEDURE swap_i
   END INTERFACE

   INTERFACE array_eq
      MODULE PROCEDURE array_eq_i8
      MODULE PROCEDURE array_eq_i
   END INTERFACE

CONTAINS

   SUBROUTINE swap_i8(arr)
      INTEGER(KIND=int_8), DIMENSION(2), INTENT(INOUT)   :: arr

      INTEGER(KIND=int_8)                                :: tmp

      tmp = arr(1)
      arr(1) = arr(2)
      arr(2) = tmp
   END SUBROUTINE

   SUBROUTINE swap_i(arr)
      INTEGER, DIMENSION(2), INTENT(INOUT)               :: arr

      INTEGER                                            :: tmp

      tmp = arr(1)
      arr(1) = arr(2)
      arr(2) = tmp
   END SUBROUTINE

   SUBROUTINE index_unique(index_in, index_out)
      !! Get all unique elements in index_in
      INTEGER, DIMENSION(:, :), INTENT(IN)   :: index_in
      INTEGER, ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: index_out

      INTEGER                                            :: blk, count, orig_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)  :: index_tmp
      INTEGER, DIMENSION(2)                  :: prev_index
      INTEGER, DIMENSION(1:SIZE(index_in, 1) &
                         , 1:SIZE(index_in, 2))                          :: index_sorted

      orig_size = SIZE(index_in, 1)
      ALLOCATE (index_tmp(orig_size, 2))
      index_sorted(:, :) = index_in(:, :)
      CALL sort_indices(orig_size, index_sorted(:, 1), index_sorted(:, 2))
      count = 0
      prev_index(:) = [0, 0]
      DO blk = 1, orig_size
         IF (ANY(index_sorted(blk, :) .NE. prev_index(:))) THEN
            count = count + 1
            index_tmp(count, :) = index_sorted(blk, :)
            prev_index(:) = index_sorted(blk, :)
         ENDIF
      ENDDO

      ALLOCATE (index_out(count, 2))
      index_out(:, :) = index_tmp(1:count, :)
   END SUBROUTINE

   !================================================================================================
   !TODO: Simplify using array_sort.fypp and merge with index_unique above.
   !================================================================================================
   SUBROUTINE sort_indices(n, row_i, col_i, blk_p, blk_d)
      !! Sorts the rows & columns of a work matrix
      !!
      !! Description
      !! Sorts the row and column indices so that the rows monotonically
      !! increase and the columns monotonically increase within each row.
      !! Passing the blk_p array rearranges the block pointers accordingly.
      !! This must be done if they are pointing to valid data, otherwise
      !! they become invalid.

      INTEGER, INTENT(IN)                                :: n
         !! number of blocks (elements) to sort
      INTEGER, DIMENSION(1:), INTENT(INOUT)              :: row_i, col_i
         !! row indices
         !! column indices
      INTEGER, DIMENSION(1:), INTENT(INOUT), OPTIONAL    :: blk_p, blk_d
         !! block pointers
         !! data storage

      CHARACTER(len=*), PARAMETER :: routineN = 'sort_indices'
      INTEGER(KIND=int_8), PARAMETER                     :: lmask8 = 4294967295_int_8

      INTEGER                                            :: error_handle, i
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: sort_keys
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: buf, buf_d

!   ---------------------------------------------------------------------------

      IF (SIZE(row_i) .EQ. 0) RETURN

      CALL timeset(routineN, error_handle)

      IF (SIZE(row_i) < n) CPABORT('row_i too small')
      IF (SIZE(col_i) < n) CPABORT('col_i too small')
      IF (PRESENT(blk_p)) THEN
         IF (SIZE(blk_p) < n) CPABORT('blk_p too small')
         ALLOCATE (buf(n))
         buf(1:n) = blk_p(1:n)
      ENDIF
      IF (PRESENT(blk_d)) THEN
         ALLOCATE (buf_d(n))
         buf_d(1:n) = blk_d(1:n)
      ENDIF
      ! Create an ordering for both rows and columns. If the blk_p must
      ! be rearranged, then the col_i array will be used as a
      ! permutation vector.
      ALLOCATE (sort_keys(n))
      sort_keys(:) = IOR(ISHFT(INT(row_i(1:n), int_8), 32), INT(col_i(1:n), int_8))
      IF (PRESENT(blk_p)) col_i(1:n) = (/(i, i=1, n)/)
      ! Now do a nice quicksort.
      CALL sort(sort_keys, n, col_i)
      ! Since blk_d is usually not present we can have two loops that
      ! are essentially the same.
      IF (PRESENT(blk_p)) THEN
         DO i = 1, n
            blk_p(i) = buf(col_i(i))
         ENDDO
         DEALLOCATE (buf)
      END IF
      IF (PRESENT(blk_d)) THEN
         DO i = 1, n
            blk_d(i) = buf_d(col_i(i))
         ENDDO
         DEALLOCATE (buf_d)
      ENDIF
      DO i = 1, n
         col_i(i) = INT(IAND(sort_keys(i), lmask8), int_4)
         row_i(i) = INT(ISHFT(sort_keys(i), -32), int_4)
      ENDDO
      DEALLOCATE (sort_keys)
      CALL timestop(error_handle)

   END SUBROUTINE sort_indices


   PURE_ARRAY_EQ FUNCTION array_eq_i(arr1, arr2)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: arr1, arr2
      LOGICAL                                            :: array_eq_i
! #if TO_VERSION(1, 11) <= TO_VERSION(LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
!       array_eq_i = .NOT. libxsmm_diff(arr1, arr2)
! #else
      array_eq_i = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i = ALL(arr1 == arr2)
! #endif
   END FUNCTION

   PURE_ARRAY_EQ FUNCTION array_eq_i8(arr1, arr2)
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: arr1, arr2
      LOGICAL                                            :: array_eq_i8
! #if TO_VERSION(1, 11) <= TO_VERSION(LIBXSMM_CONFIG_VERSION_MAJOR, LIBXSMM_CONFIG_VERSION_MINOR)
!       array_eq_i8 = .NOT. libxsmm_diff(arr1, arr2)
! #else
      array_eq_i8 = .FALSE.
      IF (SIZE(arr1) .EQ. SIZE(arr2)) array_eq_i8 = ALL(arr1 == arr2)
! #endif
   END FUNCTION

END MODULE

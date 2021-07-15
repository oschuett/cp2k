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

   USE block_matrix_api, ONLY: block_matrix_transpose, &
                               block_matrix_no_transpose, &
                               block_matrix_sort_indices

   USE kinds, ONLY: int_8

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
      invert_transpose_flag, &
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
      CALL block_matrix_sort_indices(orig_size, index_sorted(:, 1), index_sorted(:, 2))
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

   SUBROUTINE invert_transpose_flag(trans_flag)
      CHARACTER(LEN=1), INTENT(INOUT)                    :: trans_flag

      IF (trans_flag == block_matrix_transpose) THEN
         trans_flag = block_matrix_no_transpose
      ELSEIF (trans_flag == block_matrix_no_transpose) THEN
         trans_flag = block_matrix_transpose
      ENDIF
   END SUBROUTINE

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

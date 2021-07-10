!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE block_tas_types
   !! DBCSR tall-and-skinny base types.
   !! Mostly wrappers around existing DBCSR routines.

   USE block_tas_global, ONLY: &
      block_tas_distribution, block_tas_rowcol_data
   USE block_matrix_api, ONLY: &
      block_matrix_distribution_obj, block_matrix_iterator, block_matrix_type, block_scalar_type
   USE kinds, ONLY: int_8, dp

#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_tas_types'

   PUBLIC :: &
      block_tas_distribution_type, &
      block_tas_iterator, &
      block_tas_split_info, &
      block_tas_type, &
      block_tas_mm_storage

   ! info on MPI Cartesian grid that is split on MPI subgroups.
   ! info on distribution of matrix rows / columns to different subgroups.
   TYPE block_tas_split_info
      INTEGER :: mp_comm ! global communicator
      INTEGER, DIMENSION(2) :: pdims ! dimensions of process grid
      INTEGER :: igroup ! which subgroup do I belong to
      INTEGER :: ngroup ! how many groups in total
      INTEGER :: split_rowcol ! split row or column?
      INTEGER :: pgrid_split_size ! how many process rows/cols in subgroups
      INTEGER :: group_size ! group size (how many cores) of subgroups
      INTEGER :: mp_comm_group ! sub communicator
      INTEGER, ALLOCATABLE :: ngroup_opt ! optimal number of groups (split factor)
      LOGICAL, DIMENSION(2) :: strict_split = [.FALSE., .FALSE.]
      ! if .true., split factor should not be modified (2 parameters for current and general settings)
      INTEGER, POINTER :: refcount => NULL() ! lightweight reference counting for communicators
   END TYPE

   TYPE block_tas_distribution_type
      TYPE(block_tas_split_info) :: info
      TYPE(block_matrix_distribution_obj) :: block_matrix_dist
      CLASS(block_tas_distribution), ALLOCATABLE :: row_dist
      CLASS(block_tas_distribution), ALLOCATABLE :: col_dist
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:) :: local_rowcols
   END TYPE

   ! storage for batched matrix multiplication
   TYPE block_tas_mm_storage
      TYPE(block_tas_type), POINTER :: store_batched => NULL() ! intermediate replicated matrix
      TYPE(block_tas_type), POINTER :: store_batched_repl => NULL() ! intermediate replicated matrix
      LOGICAL :: batched_out ! whether replicated matrix has been changed in mm and should be copied to actual matrix
      LOGICAL :: batched_trans
      REAL(dp) :: batched_beta
   END TYPE

   ! type for tall-and-skinny matrices
   TYPE block_tas_type
      TYPE(block_tas_distribution_type)  :: dist
      CLASS(block_tas_rowcol_data), ALLOCATABLE :: row_blk_size
      CLASS(block_tas_rowcol_data), ALLOCATABLE :: col_blk_size

      TYPE(block_matrix_type) :: matrix ! matrix on subgroup
      INTEGER(KIND=int_8) :: nblkrows ! total number of rows
      INTEGER(KIND=int_8) :: nblkcols ! total number of columns
      INTEGER(KIND=int_8) :: nblkrowscols_split ! nblkrows or nblkcols depending on which is splitted
      INTEGER(KIND=int_8) :: nfullrows ! total number of full (not blocked) rows
      INTEGER(KIND=int_8) :: nfullcols ! total number of full (not blocked) columns
      LOGICAL :: valid = .FALSE. ! has been created?

      ! storage and flags for batched matrix multiplication
      INTEGER :: do_batched = 0 ! state flag for batched multiplication
      TYPE(block_tas_mm_storage), ALLOCATABLE :: mm_storage ! storage for batched processing of matrix matrix multiplication.
      LOGICAL :: has_opt_pgrid = .FALSE. ! whether pgrid was automatically optimized
   END TYPE

   TYPE block_tas_iterator
      TYPE(block_tas_split_info) :: info
      TYPE(block_tas_distribution_type) :: dist
      TYPE(block_matrix_iterator) :: iter
   END TYPE block_tas_iterator

END MODULE

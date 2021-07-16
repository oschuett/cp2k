!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE block_tas_base
   !! Tall-and-skinny matrices: base routines similar to DBCSR API, mostly wrappers around existing
   !! DBCSR routines.

   USE block_matrix_api, ONLY: &
      block_matrix_row_block_sizes, block_matrix_col_block_sizes, &
      block_matrix_get_block_p, block_matrix_put_block, block_matrix_reserve_blocks, &
      block_matrix_distribution_col_dist, block_matrix_distribution_new, &
      block_matrix_distribution_row_dist, block_matrix_distribution_hold, &
      block_matrix_iterator_blocks_left, block_matrix_iterator_next_block, &
      block_matrix_iterator_start, block_matrix_iterator_stop, &
      block_matrix_distribution_release, block_matrix_release, &
      block_matrix_get_num_blocks, block_matrix_get_nze, &
      block_matrix_get_info, block_matrix_set, block_matrix_filter, block_matrix_clear, &
      block_matrix_distribution_obj, block_matrix_iterator, block_matrix_type, &
      block_matrix_create, block_matrix_finalize, block_matrix_get_name
   USE block_tas_types, ONLY: &
      block_tas_distribution_type, block_tas_iterator, block_tas_split_info, block_tas_type
   USE block_tas_global, ONLY: &
      block_tas_blk_size_arb, block_tas_dist_arb, block_tas_distribution, block_tas_rowcol_data
   USE block_tas_split, ONLY: &
      block_index_global_to_local, block_index_local_to_global, colsplit, &
      block_tas_info_hold, block_tas_release_info, block_tas_create_split, &
      group_to_mrowcol, rowsplit, block_tas_get_split_info
   USE block_tas_util, ONLY: index_unique
   USE kinds, ONLY: &
      default_string_length, int_8, real_8, real_4, dp
   USE message_passing, ONLY: &
      mp_cart_rank, mp_environ, mp_sum, mp_max
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_tas_base'

   PUBLIC :: &
      ! DBCSR wrappers / interface routines
      block_tas_blk_sizes, &
      block_tas_clear, &
      block_tas_copy, &
      block_tas_create, &
      block_tas_destroy, &
      block_tas_distribution_destroy, &
      block_tas_distribution_new, &
      block_tas_filter, &
      block_tas_finalize, &
      block_tas_get_block_p, &
      block_tas_get_info, &
      block_tas_get_num_blocks, &
      block_tas_get_nze, &
      block_tas_get_nze_total, &
      block_tas_get_num_blocks_total, &
      block_tas_get_stored_coordinates, &
      block_tas_info, &
      block_tas_iterator_blocks_left, &
      block_tas_iterator_next_block, &
      block_tas_iterator_start, &
      block_tas_iterator_stop, &
      block_tas_nblkcols_local, &
      block_tas_nblkcols_total, &
      block_tas_nblkrows_local, &
      block_tas_nblkrows_total, &
      block_tas_nfullrows_total, &
      block_tas_nfullcols_total, &
      block_tas_put_block, &
      block_tas_reserve_blocks, &
      block_tas_set, &
      block_repl_get_stored_coordinates
   PUBLIC :: &
      ! conversion routines
      block_tas_convert_to_dbcsr, &
      block_tas_convert_to_tas

   INTERFACE block_tas_create
      MODULE PROCEDURE block_tas_create_new
      MODULE PROCEDURE block_tas_create_template
   END INTERFACE

   INTERFACE block_tas_reserve_blocks
      MODULE PROCEDURE block_tas_reserve_blocks_template
      MODULE PROCEDURE block_tas_reserve_blocks_index
   END INTERFACE

   INTERFACE block_tas_iterator_next_block
      MODULE PROCEDURE block_tas_iterator_next_block_d
      MODULE PROCEDURE block_tas_iterator_next_block_index
   END INTERFACE

CONTAINS

   SUBROUTINE block_tas_create_new(matrix, name, dist, row_blk_size, col_blk_size, own_dist)
      !! Create new tall-and-skinny matrix.
      !! Exactly like block_create_new but with custom types for row_blk_size and col_blk_size instead of
      !! arrays.

      TYPE(block_tas_type), INTENT(OUT)              :: matrix
      CHARACTER(len=*), INTENT(IN)                   :: name
      TYPE(block_tas_distribution_type), INTENT(INOUT) :: dist
      CLASS(block_tas_rowcol_data), INTENT(IN)       :: row_blk_size, col_blk_size
      LOGICAL, INTENT(IN), OPTIONAL                  :: own_dist
         !! whether matrix should own distribution

      TYPE(block_tas_split_info)                     :: info

      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS     :: row_blk_size_vec, col_blk_size_vec
      INTEGER                                        :: nrows, ncols, irow, col, icol, row
      CHARACTER(LEN=*), PARAMETER                    :: routineN = 'block_tas_create_new'
      INTEGER                                        :: handle

      CALL timeset(routineN, handle)

      CALL block_tas_copy_distribution(dist, matrix%dist, own_dist)
      matrix%nblkrows = row_blk_size%nmrowcol
      matrix%nblkcols = col_blk_size%nmrowcol

      CPASSERT(matrix%nblkrows == dist%row_dist%nmrowcol)
      CPASSERT(matrix%nblkcols == dist%col_dist%nmrowcol)

      matrix%nfullrows = row_blk_size%nfullrowcol
      matrix%nfullcols = col_blk_size%nfullrowcol

      ALLOCATE (matrix%row_blk_size, source=row_blk_size)
      ALLOCATE (matrix%col_blk_size, source=col_blk_size)

      info = block_tas_info(matrix)

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         matrix%nblkrowscols_split = matrix%nblkrows

         ASSOCIATE (rows => dist%local_rowcols)
            nrows = SIZE(rows)
            ncols = INT(dist%col_dist%nmrowcol)
            ALLOCATE (row_blk_size_vec(nrows))
            ALLOCATE (col_blk_size_vec(ncols))
            DO irow = 1, nrows
               row_blk_size_vec(irow) = row_blk_size%data(rows(irow))
            ENDDO
            DO col = 1, ncols
               col_blk_size_vec(col) = col_blk_size%data(INT(col, KIND=int_8))
            ENDDO
         END ASSOCIATE
      CASE (colsplit)
         matrix%nblkrowscols_split = matrix%nblkcols

         ASSOCIATE (cols => dist%local_rowcols)
            ncols = SIZE(cols)
            nrows = INT(dist%row_dist%nmrowcol)
            ALLOCATE (row_blk_size_vec(nrows))
            ALLOCATE (col_blk_size_vec(ncols))
            DO icol = 1, ncols
               col_blk_size_vec(icol) = col_blk_size%data(cols(icol))
            ENDDO
            DO row = 1, nrows
               row_blk_size_vec(row) = row_blk_size%data(INT(row, KIND=int_8))
            ENDDO
         END ASSOCIATE
      END SELECT

      CALL block_matrix_create(matrix=matrix%matrix, &
                        name=name, &
                        dist=dist%block_matrix_dist, &
                        row_blk_size=row_blk_size_vec, &
                        col_blk_size=col_blk_size_vec)

      DEALLOCATE(row_blk_size_vec, col_blk_size_vec)
      matrix%valid = .TRUE.
      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE block_tas_create_template(matrix_in, matrix, name)
      !! Create matrix from template
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix_in
      TYPE(block_tas_type), INTENT(OUT)                    :: matrix
      CHARACTER(len=*), INTENT(IN), OPTIONAL               :: name

      CALL block_matrix_create(matrix%matrix, template=matrix_in%matrix, name=name)
      CALL block_matrix_finalize(matrix%matrix)

      CALL block_tas_copy_distribution(matrix_in%dist, matrix%dist)
      ALLOCATE (matrix%row_blk_size, source=matrix_in%row_blk_size)
      ALLOCATE (matrix%col_blk_size, source=matrix_in%col_blk_size)
      matrix%nblkrows = matrix_in%nblkrows
      matrix%nblkcols = matrix_in%nblkcols
      matrix%nblkrowscols_split = matrix_in%nblkrowscols_split
      matrix%nfullrows = matrix_in%nfullrows
      matrix%nfullcols = matrix_in%nfullcols
      matrix%valid = .TRUE.

   END SUBROUTINE

   SUBROUTINE block_tas_destroy(matrix)
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix

      CALL block_matrix_release(matrix%matrix)
      CALL block_tas_distribution_destroy(matrix%dist)
      DEALLOCATE (matrix%row_blk_size)
      DEALLOCATE (matrix%col_blk_size)
      matrix%valid = .FALSE.
   END SUBROUTINE

   SUBROUTINE block_tas_copy(matrix_b, matrix_a, summation)
      !! Copy matrix_a to matrix_b

      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix_b
      TYPE(block_tas_type), INTENT(IN)                     :: matrix_a
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation
         !! Whether to sum matrices b = a + b

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_copy'
      INTEGER                                            :: handle
      TYPE(block_tas_iterator)                             :: iter
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      INTEGER(KIND=int_8)                                :: row, column

      CALL timeset(routineN, handle)
      CPASSERT(matrix_b%valid)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL block_tas_clear(matrix_b)
      ELSE
         CALL block_tas_clear(matrix_b)
      ENDIF

      CALL block_tas_reserve_blocks(matrix_a, matrix_b)

      CALL block_tas_iterator_start(iter, matrix_a)
      DO WHILE (block_tas_iterator_blocks_left(iter))
         CALL block_tas_iterator_next_block(iter, row, column, block)
         CALL block_tas_put_block(matrix_b, row, column, block, summation=summation)
      ENDDO
      CALL block_tas_iterator_stop(iter)
      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_tas_reserve_blocks_template(matrix_in, matrix_out)
      !! Make sure that matrix_out has same blocks reserved as matrix_in. This assumes that both
      !! matrices have same number of block rows and block columns.

      TYPE(block_tas_type), INTENT(IN) :: matrix_in
      TYPE(block_tas_type), INTENT(INOUT) :: matrix_out
      TYPE(block_tas_iterator)                        :: iter
      INTEGER :: iblk, nblk
      INTEGER(KIND=int_8)                  :: row, column
      INTEGER(KIND=int_8), DIMENSION(:), ALLOCATABLE :: row_res, col_res

      nblk = block_tas_get_num_blocks(matrix_in)
      ALLOCATE (row_res(nblk), col_res(nblk))
      CALL block_tas_iterator_start(iter, matrix_in)
      DO iblk = 1, nblk
         CALL block_tas_iterator_next_block(iter, row, column)
         row_res(iblk) = row
         col_res(iblk) = column
      ENDDO
      CPASSERT(.NOT. block_tas_iterator_blocks_left(iter))
      CALL block_tas_iterator_stop(iter)

      CALL block_tas_reserve_blocks(matrix_out, row_res, col_res)
   END SUBROUTINE

   SUBROUTINE block_tas_finalize(matrix)
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix

      CALL block_matrix_finalize(matrix%matrix)
   END SUBROUTINE

   SUBROUTINE block_tas_distribution_new(dist, mp_comm, &
                                         row_dist, col_dist, split_info, nosplit)
      !! create new distribution.
      !! Exactly like block_matrix_distribution_new but with custom types for row_dist and col_dist instead of
      !! arrays.

      TYPE(block_tas_distribution_type), INTENT(OUT)   :: dist
      INTEGER, INTENT(IN)                              :: mp_comm
      CLASS(block_tas_distribution), INTENT(IN)        :: row_dist, col_dist
      TYPE(block_tas_split_info), INTENT(IN), OPTIONAL :: split_info
         !! Strategy of how to split process grid (optional). If not present a default split heuristic is applied.
      LOGICAL, INTENT(IN), OPTIONAL                    :: nosplit
         !! if .TRUE. don't split process grid (optional)
      !LOGICAL, INTENT(IN), OPTIONAL                    :: strict_split

      TYPE(block_tas_split_info)                       :: split_info_prv

      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS       :: row_dist_vec
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS       :: col_dist_vec
      INTEGER                                          :: nrows, ncols, irow, col, icol, row, &
                                                          split_rowcol, nsplit, handle
      LOGICAL                                          :: opt_nsplit
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_distribution_new'

      CALL timeset(routineN, handle)
      IF (PRESENT(split_info)) THEN
         CALL block_tas_info_hold(split_info)
         split_info_prv = split_info
      ELSE
         ! default split heuristic: split into submatrices that have roughly same block dimensions
         IF (row_dist%nmrowcol >= col_dist%nmrowcol) THEN
            split_rowcol = rowsplit
            nsplit = INT((row_dist%nmrowcol - 1)/col_dist%nmrowcol + 1)
         ELSE
            split_rowcol = colsplit
            nsplit = INT((col_dist%nmrowcol - 1)/row_dist%nmrowcol + 1)
         ENDIF
         opt_nsplit = .TRUE.
         IF (PRESENT(nosplit)) THEN
            IF (nosplit) THEN
               nsplit = 1
               opt_nsplit = .FALSE.
            ENDIF
         ENDIF
         CALL block_tas_create_split(split_info_prv, mp_comm, split_rowcol, nsplit=nsplit, opt_nsplit=opt_nsplit)
      ENDIF

      SELECT CASE (split_info_prv%split_rowcol)
      CASE (rowsplit)
         CALL group_to_mrowcol(split_info_prv, row_dist, split_info_prv%igroup, dist%local_rowcols)
         nrows = SIZE(dist%local_rowcols)
         ncols = INT(col_dist%nmrowcol)
         ALLOCATE (row_dist_vec(nrows))
         ALLOCATE (col_dist_vec(ncols))
         DO irow = 1, nrows
            row_dist_vec(irow) = row_dist%dist(dist%local_rowcols(irow)) - split_info_prv%pgrid_split_size*split_info_prv%igroup
         ENDDO
         DO col = 1, ncols
            col_dist_vec(col) = col_dist%dist(INT(col, KIND=int_8))
         ENDDO
      CASE (colsplit)
         CALL group_to_mrowcol(split_info_prv, col_dist, split_info_prv%igroup, dist%local_rowcols)
         ncols = SIZE(dist%local_rowcols)
         nrows = INT(row_dist%nmrowcol)
         ALLOCATE (col_dist_vec(ncols))
         ALLOCATE (row_dist_vec(nrows))
         DO icol = 1, ncols
            col_dist_vec(icol) = col_dist%dist(dist%local_rowcols(icol)) - split_info_prv%pgrid_split_size*split_info_prv%igroup
         ENDDO
         DO row = 1, nrows
            row_dist_vec(row) = row_dist%dist(INT(row, KIND=int_8))
         ENDDO
      END SELECT

      dist%info = split_info_prv

      CALL block_matrix_distribution_new(dist%block_matrix_dist, split_info_prv%mp_comm_group, &
                                         row_dist_vec, col_dist_vec)
      DEALLOCATE(row_dist_vec, col_dist_vec)
      ALLOCATE (dist%row_dist, source=row_dist)
      ALLOCATE (dist%col_dist, source=col_dist)

      !IF(PRESENT(strict_split)) dist%strict_split = strict_split

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_tas_distribution_destroy(dist)
      TYPE(block_tas_distribution_type), INTENT(INOUT)     :: dist

      ! Note: Issue with Cray CCE compiler
      ! commented out the following deallocate statements on polymorphic variables,
      ! these cause segfaults with CCE compiler at a later point

      !IF (ALLOCATED(dist%row_dist)) THEN
      !   DEALLOCATE (dist%row_dist)
      !ENDIF
      !IF (ALLOCATED(dist%col_dist)) THEN
      !   DEALLOCATE (dist%col_dist)
      !ENDIF

      IF (ALLOCATED(dist%local_rowcols)) THEN
         DEALLOCATE (dist%local_rowcols)
      ENDIF
      CALL block_tas_release_info(dist%info)
      CALL block_matrix_distribution_release(dist%block_matrix_dist)
   END SUBROUTINE

   SUBROUTINE block_tas_get_stored_coordinates(matrix, row, column, processor)
      !! As block_get_stored_coordinates

      TYPE(block_tas_type), INTENT(IN)                   :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                    :: row, column
         !! global matrix blocked row
         !! global matrix blocked column
      INTEGER, INTENT(OUT)                               :: processor
         !! process ID

      INTEGER, DIMENSION(2)                              :: pcoord
      TYPE(block_tas_split_info)                         :: info

      pcoord(1) = matrix%dist%row_dist%dist(row)
      pcoord(2) = matrix%dist%col_dist%dist(column)
      info = block_tas_info(matrix)

      ! workaround for inefficient mpi_cart_rank
      processor = pcoord(1)*info%pdims(2) + pcoord(2)

   END SUBROUTINE

   SUBROUTINE block_repl_get_stored_coordinates(matrix, row, column, processors)
      !! Get all processors for a given row/col combination if matrix is replicated on each process
      !! subgroup.

      TYPE(block_tas_type), INTENT(IN)                     :: matrix
         !! tall-and-skinny matrix whose DBCSR submatrices are replicated matrices
      INTEGER, INTENT(IN)                                :: row, column
         !! row of a submatrix
         !! column of a submatrix
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: processors

      INTEGER                                            :: igroup
      INTEGER(KIND=int_8)                                :: col_s, row_s
      INTEGER, DIMENSION(2)                              :: pcoord
      TYPE(block_tas_split_info)                         :: info

      row_s = INT(row, KIND=int_8); col_s = INT(column, KIND=int_8)

      info = block_tas_info(matrix)
      pcoord(1) = matrix%dist%row_dist%dist(row_s)
      pcoord(2) = matrix%dist%col_dist%dist(col_s)

      DO igroup = 0, info%ngroup - 1
         CALL mp_cart_rank(info%mp_comm, pcoord, processors(igroup + 1))
         SELECT CASE (info%split_rowcol)
         CASE (rowsplit)
            row_s = row_s + block_tas_nblkrows_local(matrix)
            pcoord(1) = matrix%dist%row_dist%dist(row_s)
         CASE (colsplit)
            col_s = col_s + block_tas_nblkcols_local(matrix)
            pcoord(2) = matrix%dist%col_dist%dist(col_s)
         END SELECT
      ENDDO
   END SUBROUTINE

   SUBROUTINE block_tas_convert_to_dbcsr(matrix_rect, matrix_dbcsr)
      !! Convert a tall-and-skinny matrix into a normal DBCSR matrix.
      !! This is not recommended for matrices with a very large dimension.

      TYPE(block_tas_type), INTENT(IN)                     :: matrix_rect
      TYPE(block_matrix_type), INTENT(OUT)                      :: matrix_dbcsr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_convert_to_dbcsr'

      INTEGER                                            :: handle
      INTEGER(KIND=int_8)                                :: col, row
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: col_dist_vec, col_size_vec, &
                                                            row_dist_vec, row_size_vec
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(block_matrix_distribution_obj)                       :: dist
      TYPE(block_tas_iterator)                           :: iter
      TYPE(block_tas_split_info)                         :: info
      INTEGER                                            :: rb_count, nblks_local
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nz_rows, nz_cols

      CALL timeset(routineN, handle)

      info = block_tas_info(matrix_rect)

      ALLOCATE (row_dist_vec(matrix_rect%nblkrows))
      ALLOCATE (row_size_vec(matrix_rect%nblkrows))
      ALLOCATE (col_dist_vec(matrix_rect%nblkcols))
      ALLOCATE (col_size_vec(matrix_rect%nblkcols))

      DO row = 1, matrix_rect%nblkrows
         row_dist_vec(row) = matrix_rect%dist%row_dist%dist(row)
         row_size_vec(row) = matrix_rect%row_blk_size%data(row)
      ENDDO

      DO col = 1, matrix_rect%nblkcols
         col_dist_vec(col) = matrix_rect%dist%col_dist%dist(col)
         col_size_vec(col) = matrix_rect%col_blk_size%data(col)
      ENDDO

      CALL block_matrix_distribution_new(dist, info%mp_comm, row_dist_vec, col_dist_vec)
      DEALLOCATE(row_dist_vec, col_dist_vec)

      CALL block_matrix_create(matrix=matrix_dbcsr, &
                        name=TRIM(block_matrix_get_name(matrix_rect%matrix)), &
                        dist=dist, &
                        row_blk_size=row_size_vec, &
                        col_blk_size=col_size_vec)

      CALL block_matrix_distribution_release(dist)

      DEALLOCATE(row_size_vec, col_size_vec)
      nblks_local = block_tas_get_num_blocks(matrix_rect)
      CALL block_tas_iterator_start(iter, matrix_rect)
      ALLOCATE (nz_rows(nblks_local), nz_cols(nblks_local))
      rb_count = 0
      DO WHILE (block_tas_iterator_blocks_left(iter))
         CALL block_tas_iterator_next_block(iter, row, col)
         rb_count = rb_count + 1
         nz_rows(rb_count) = INT(row)
         nz_cols(rb_count) = INT(col)
      ENDDO

      CALL block_matrix_reserve_blocks(matrix_dbcsr, nz_rows, nz_cols)

      CALL block_tas_iterator_stop(iter)
      CALL block_tas_iterator_start(iter, matrix_rect)
      DO WHILE (block_tas_iterator_blocks_left(iter))
         CALL block_tas_iterator_next_block(iter, row, col, block)
         CALL block_matrix_put_block(matrix_dbcsr, INT(row), INT(col), block)
      ENDDO

      CALL block_tas_iterator_stop(iter)

      CALL block_matrix_finalize(matrix_dbcsr)

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_tas_convert_to_tas(info, matrix_rect, matrix_dbcsr)
      !! Converts a DBCSR matrix into the tall-and-skinny matrix type.

      TYPE(block_tas_split_info), INTENT(IN)               :: info
         !! Strategy of how to split process grid
      TYPE(block_tas_type), INTENT(OUT)                    :: matrix_rect
      TYPE(block_matrix_type), INTENT(IN)                         :: matrix_dbcsr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_convert_to_tas'

      CHARACTER(len=default_string_length)               :: name
      INTEGER                                            :: col, handle, numnodes, row
      INTEGER(KIND=int_8)                                :: nbcols, nbrows
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: col_blk_size, row_blk_size
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(block_matrix_distribution_obj)                       :: block_matrix_dist
      TYPE(block_matrix_iterator)                               :: iter
      TYPE(block_tas_blk_size_arb)                         :: col_blk_size_obj, row_blk_size_obj
      TYPE(block_tas_dist_arb)                             :: col_dist_obj, row_dist_obj
      TYPE(block_tas_distribution_type)                    :: dist

      NULLIFY (col_blk_size, row_blk_size)
      CALL timeset(routineN, handle)
      CALL mp_environ(numnodes, pdims, pcoord, info%mp_comm)
      CALL block_matrix_get_info(matrix_dbcsr, distribution=block_matrix_dist, name=name, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size)

      nbrows = SIZE(block_matrix_row_block_sizes(matrix_dbcsr))
      nbcols = SIZE(block_matrix_col_block_sizes(matrix_dbcsr))

      row_dist_obj = block_tas_dist_arb(block_matrix_distribution_row_dist(block_matrix_dist), pdims(1), nbrows)
      col_dist_obj = block_tas_dist_arb(block_matrix_distribution_col_dist(block_matrix_dist), pdims(2), nbcols)

      row_blk_size_obj = block_tas_blk_size_arb(row_blk_size)
      col_blk_size_obj = block_tas_blk_size_arb(col_blk_size)

      CALL block_tas_distribution_new(dist, info%mp_comm, row_dist_obj, col_dist_obj)

      CALL block_tas_create(matrix_rect, TRIM(name)//"_compressed", &
                            dist, row_blk_size_obj, col_blk_size_obj)

      CALL block_matrix_iterator_start(iter, matrix_dbcsr)
      DO WHILE (block_matrix_iterator_blocks_left(iter))
         CALL block_matrix_iterator_next_block(iter, row, col, block)
         CALL block_tas_put_block(matrix_rect, INT(row, KIND=int_8), INT(col, KIND=int_8), block)
      ENDDO
      CALL block_matrix_iterator_stop(iter)

      CALL block_tas_finalize(matrix_rect)

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_tas_iterator_start(iter, matrix_in)
      !! As block_matrix_iterator_start
      TYPE(block_tas_iterator), INTENT(INOUT)              :: iter
      TYPE(block_tas_type), INTENT(IN)                     :: matrix_in

      CALL block_matrix_iterator_start(iter%iter, matrix_in%matrix)

      iter%dist = matrix_in%dist
   END SUBROUTINE

   FUNCTION block_tas_iterator_blocks_left(iter)
      !! As block_matrix_iterator_blocks_left
      TYPE(block_tas_iterator), INTENT(IN)                 :: iter
      LOGICAL                                            :: block_tas_iterator_blocks_left

      block_tas_iterator_blocks_left = block_matrix_iterator_blocks_left(iter%iter)
   END FUNCTION

   SUBROUTINE block_tas_iterator_stop(iter)
      !! As block_matrix_iterator_stop
      TYPE(block_tas_iterator), INTENT(INOUT)              :: iter

      CALL block_matrix_iterator_stop(iter%iter)
   END SUBROUTINE

   SUBROUTINE block_tas_iterator_next_block_index(iterator, row, column, &
                                                  blk_p, row_size, col_size)
      !! As block_matrix_iterator_next_block

      TYPE(block_tas_iterator), INTENT(INOUT)              :: iterator
      INTEGER(KIND=int_8), INTENT(OUT)                   :: row, column
         !! global block row
         !! global block column
      INTEGER, INTENT(OUT), OPTIONAL                     :: blk_p, row_size, col_size

      INTEGER                                            :: column_group, row_group

      CALL block_matrix_iterator_next_block(iterator%iter, row_group, column_group, blk_p, row_size, col_size)

      CALL block_index_local_to_global(iterator%dist%info, iterator%dist, row_group=row_group, column_group=column_group, &
                                       row=row, column=column)

   END SUBROUTINE

   SUBROUTINE block_tas_reserve_blocks_index(matrix, rows, columns)
      !! As block_matrix_reserve_blocks
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix
      INTEGER(KIND=int_8), DIMENSION(:), INTENT(IN)      :: rows, columns

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_reserve_blocks_index'

      INTEGER                                            :: handle, icol, irow
      INTEGER, DIMENSION(SIZE(rows), 2)                :: rowcols_group
      INTEGER, DIMENSION(:, :), ALLOCATABLE             :: rowcols_group_u

      CALL timeset(routineN, handle)

      DO irow = 1, SIZE(rows)
         CALL block_index_global_to_local(block_tas_info(matrix), matrix%dist, row=rows(irow), &
                                          row_group=rowcols_group(irow, 1))
      ENDDO

      DO icol = 1, SIZE(columns)
         CALL block_index_global_to_local(block_tas_info(matrix), matrix%dist, column=columns(icol), &
                                          column_group=rowcols_group(icol, 2))
      ENDDO

      CALL index_unique(rowcols_group, rowcols_group_u) ! make sure that index is unique, not sure
      ! if this is really needed or whether DBCSR
      ! takes care of duplicate indices

      CALL block_matrix_reserve_blocks(matrix%matrix, rowcols_group_u(:, 1), rowcols_group_u(:, 2))

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_tas_copy_distribution(dist_in, dist_out, own_dist)
      !! Copy a distribution

      TYPE(block_tas_distribution_type), INTENT(INOUT)     :: dist_in
      TYPE(block_tas_distribution_type), INTENT(OUT)       :: dist_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_dist
         !! Whether distribution should be owned by dist_out

      LOGICAL                                            :: own_dist_prv

      IF (PRESENT(own_dist)) THEN
         own_dist_prv = own_dist
      ELSE
         own_dist_prv = .FALSE.
      ENDIF

      IF (.NOT. own_dist_prv) THEN
         CALL block_matrix_distribution_hold(dist_in%block_matrix_dist)
         CALL block_tas_info_hold(dist_in%info)
      ENDIF

      dist_out = dist_in
   END SUBROUTINE

   SUBROUTINE block_tas_blk_sizes(matrix, row, col, row_size, col_size)
      !! Get block size for a given row & column
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8), INTENT(IN)                    :: row, col
      INTEGER, INTENT(OUT)                               :: row_size, col_size
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_blk_sizes'
      INTEGER :: handle

      CALL timeset(routineN, handle)

      row_size = matrix%row_blk_size%data(row)
      col_size = matrix%col_blk_size%data(col)

      CALL timestop(handle)
   END SUBROUTINE

   FUNCTION block_tas_info(matrix)
      !! get info on mpi grid splitting
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      TYPE(block_tas_split_info)                           :: block_tas_info

      block_tas_info = matrix%dist%info
   END FUNCTION

   FUNCTION block_tas_nblkrows_total(matrix) RESULT(nblkrows_total)
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nblkrows_total

      nblkrows_total = matrix%nblkrows
   END FUNCTION

   FUNCTION block_tas_nfullrows_total(matrix) RESULT(nfullrows_total)
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nfullrows_total

      nfullrows_total = matrix%nfullrows
   END FUNCTION

   FUNCTION block_tas_nblkcols_total(matrix) RESULT(nblkcols_total)
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nblkcols_total

      nblkcols_total = matrix%nblkcols
   END FUNCTION

   FUNCTION block_tas_nfullcols_total(matrix) RESULT(nfullcols_total)
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER(KIND=int_8)                                :: nfullcols_total

      nfullcols_total = matrix%nfullcols
   END FUNCTION

   FUNCTION block_tas_nblkcols_local(matrix) RESULT(nblkcols_local)
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER                                            :: nblkcols_local

      nblkcols_local = SIZE(block_matrix_col_block_sizes(matrix%matrix))
   END FUNCTION

   FUNCTION block_tas_nblkrows_local(matrix) RESULT(nblkrows_local)
      TYPE(block_tas_type), INTENT(IN)                     :: matrix
      INTEGER                                            :: nblkrows_local

      nblkrows_local = SIZE(block_matrix_row_block_sizes(matrix%matrix))
   END FUNCTION

   PURE FUNCTION block_tas_get_num_blocks(matrix) RESULT(num_blocks)
      !! As block_get_num_blocks: get number of local blocks
      TYPE(block_tas_type), INTENT(IN) :: matrix
      INTEGER                        :: num_blocks

      num_blocks = block_matrix_get_num_blocks(matrix%matrix)
   END FUNCTION

   FUNCTION block_tas_get_num_blocks_total(matrix) RESULT(num_blocks)
      !! get total number of blocks
      TYPE(block_tas_type), INTENT(IN) :: matrix
      INTEGER(KIND=int_8)              :: num_blocks
      TYPE(block_tas_split_info)       :: info

      info = block_tas_info(matrix)
      num_blocks = block_tas_get_num_blocks(matrix)
      CALL mp_sum(num_blocks, info%mp_comm)

   END FUNCTION

   PURE FUNCTION block_tas_get_nze(matrix)
      !! As block_get_nze: get number of local non-zero elements
      TYPE(block_tas_type), INTENT(IN) :: matrix
      INTEGER :: block_tas_get_nze

      block_tas_get_nze = block_matrix_get_nze(matrix%matrix)

   END FUNCTION

   FUNCTION block_tas_get_nze_total(matrix)
      !! Get total number of non-zero elements
      TYPE(block_tas_type), INTENT(IN) :: matrix
      INTEGER(KIND=int_8)              :: block_tas_get_nze_total
      TYPE(block_tas_split_info)       :: info

      block_tas_get_nze_total = block_tas_get_nze(matrix)
      info = block_tas_info(matrix)
      CALL mp_sum(block_tas_get_nze_total, info%mp_comm)
   END FUNCTION

   SUBROUTINE block_tas_clear(matrix)
      !! Clear matrix (erase all data)
      TYPE(block_tas_type), INTENT(INOUT) :: matrix

      CALL block_matrix_clear(matrix%matrix)
   END SUBROUTINE

   SUBROUTINE block_tas_get_info(matrix, nblkrows_total, nblkcols_total, &
                                 nfullrows_total, nfullcols_total, &
                                 nblkrows_local, nblkcols_local, &
                                 nfullrows_local, nfullcols_local, &
                                 nprow, npcol, my_prow, my_pcol, &
                                 local_rows, local_cols, proc_row_dist, proc_col_dist, &
                                 row_blk_size, col_blk_size, distribution, name)

      TYPE(block_tas_type), INTENT(IN)                                :: matrix
      INTEGER(KIND=int_8), INTENT(OUT), OPTIONAL                      :: nblkrows_total, nblkcols_total, nfullrows_total, &
                                                                         nfullcols_total
      INTEGER, INTENT(OUT), OPTIONAL                                  :: nblkrows_local, nblkcols_local, nfullrows_local, &
                                                                         nfullcols_local, nprow, npcol, my_prow, my_pcol
      INTEGER(KIND=int_8), DIMENSION(:), OPTIONAL, ALLOCATABLE        :: local_rows, local_cols
      CLASS(block_tas_distribution), ALLOCATABLE, OPTIONAL, &
         INTENT(OUT)                                                  :: proc_row_dist, proc_col_dist
      CLASS(block_tas_rowcol_data), ALLOCATABLE, OPTIONAL, &
         INTENT(OUT)                                                  :: row_blk_size, col_blk_size
      TYPE(block_tas_distribution_type), OPTIONAL                     :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL                         :: name

      TYPE(block_tas_split_info)                                      :: info
      INTEGER                                                         :: numnodes, irow, icol
      INTEGER, DIMENSION(2)                                           :: pdims, pcoord
      INTEGER, DIMENSION(:), POINTER                                  :: local_rows_local, local_cols_local

      NULLIFY (local_rows_local, local_cols_local)
      CALL block_matrix_get_info(matrix%matrix, nblkrows_local=nblkrows_local, nblkcols_local=nblkcols_local, &
                          nfullrows_local=nfullrows_local, nfullcols_local=nfullcols_local, &
                          local_rows=local_rows_local, local_cols=local_cols_local, &
                          name=name)

      IF (PRESENT(nblkrows_total)) nblkrows_total = block_tas_nblkrows_total(matrix)
      IF (PRESENT(nblkcols_total)) nblkcols_total = block_tas_nblkcols_total(matrix)

      IF (PRESENT(nfullrows_total)) nfullrows_total = block_tas_nfullrows_total(matrix)
      IF (PRESENT(nfullcols_total)) nfullcols_total = block_tas_nfullcols_total(matrix)

      info = block_tas_info(matrix)

      CALL mp_environ(numnodes, pdims, pcoord, info%mp_comm)
      IF (PRESENT(my_prow)) my_prow = pcoord(1)
      IF (PRESENT(my_pcol)) my_pcol = pcoord(2)
      IF (PRESENT(nprow)) nprow = pdims(1)
      IF (PRESENT(npcol)) npcol = pdims(2)

      IF (PRESENT(local_rows)) THEN
         ALLOCATE (local_rows(SIZE(local_rows_local)))
         DO irow = 1, SIZE(local_rows_local)
            CALL block_index_local_to_global(info, matrix%dist, row_group=local_rows_local(irow), row=local_rows(irow))
         ENDDO
      ENDIF

      IF (PRESENT(local_cols)) THEN
         ALLOCATE (local_cols(SIZE(local_cols_local)))
         DO icol = 1, SIZE(local_cols_local)
            CALL block_index_local_to_global(info, matrix%dist, column_group=local_cols_local(icol), column=local_cols(icol))
         ENDDO
      ENDIF

      IF (PRESENT(proc_row_dist)) ALLOCATE (proc_row_dist, SOURCE=matrix%dist%row_dist)
      IF (PRESENT(proc_col_dist)) ALLOCATE (proc_col_dist, SOURCE=matrix%dist%col_dist)

      IF (PRESENT(row_blk_size)) ALLOCATE (row_blk_size, SOURCE=matrix%row_blk_size)
      IF (PRESENT(col_blk_size)) ALLOCATE (col_blk_size, SOURCE=matrix%col_blk_size)
      IF (PRESENT(distribution)) distribution = matrix%dist

   END SUBROUTINE

   SUBROUTINE block_tas_iterator_next_block_d (iterator, row, column, block, &
      !! As block_matrix_iterator_next_block
                                                         row_size, col_size)
      TYPE(block_tas_iterator), INTENT(INOUT) :: iterator
      INTEGER(KIND=int_8), INTENT(OUT)      :: row, column
      REAL(dp), DIMENSION(:, :), POINTER   :: block
      INTEGER, INTENT(OUT), OPTIONAL        :: row_size, col_size

      INTEGER                               :: row_group, column_group

      CALL block_matrix_iterator_next_block(iterator%iter, row_group, column_group, block, &
                                     row_size=row_size, col_size=col_size)

      CALL block_index_local_to_global(iterator%dist%info, iterator%dist, row_group=row_group, column_group=column_group, &
                                       row=row, column=column)

   END SUBROUTINE

   !! As block_matrix_put_block
   SUBROUTINE block_tas_put_block(matrix, row, col, block, transposed, summation, &
                                               scale)
      TYPE(block_tas_type), INTENT(INOUT)     :: matrix
      REAL(dp), DIMENSION(:, :), INTENT(IN) :: block
      LOGICAL, INTENT(IN), OPTIONAL         :: transposed, summation
      INTEGER(KIND=int_8), INTENT(IN)       :: row, col
      REAL(dp), INTENT(IN), OPTIONAL       :: scale

      INTEGER                     :: col_group, row_group

      CALL block_index_global_to_local(matrix%dist%info, matrix%dist, row=row, column=col, &
                                       row_group=row_group, column_group=col_group)

      CALL block_matrix_put_block(matrix%matrix, row_group, col_group, block, &
                                  transposed=transposed, summation=summation, scale=scale)

   END SUBROUTINE

   SUBROUTINE block_tas_get_block_p(matrix, row, col, block, transposed, found, row_size, col_size)
      !! As block_matrix_get_block_p
      TYPE(block_tas_type), INTENT(INOUT)  :: matrix
      INTEGER(KIND=int_8), INTENT(IN)    :: row, col
      REAL(dp), DIMENSION(:, :), POINTER :: block
      LOGICAL, INTENT(OUT)               :: transposed
      LOGICAL, INTENT(OUT)               :: found
      INTEGER, INTENT(OUT), OPTIONAL     :: row_size, col_size

      INTEGER                     :: col_group, row_group

      CALL block_index_global_to_local(matrix%dist%info, matrix%dist, row=row, column=col, &
                                       row_group=row_group, column_group=col_group)

      CALL block_matrix_get_block_p(matrix%matrix, row_group, col_group, block, transposed, &
                                    found, row_size=row_size, col_size=col_size)

   END SUBROUTINE

   SUBROUTINE block_tas_set(matrix, alpha)
      !! As block_matrix_set
      TYPE(block_tas_type), INTENT(INOUT)    :: matrix
      REAL(dp), INTENT(IN)                :: alpha
      CALL block_matrix_set(matrix%matrix, alpha)
   END SUBROUTINE

   SUBROUTINE block_tas_filter (matrix, eps)
      !! As block_matrix_filter
      TYPE(block_tas_type), INTENT(INOUT)    :: matrix
      REAL(dp), INTENT(IN)                :: eps

      CALL block_matrix_filter(matrix%matrix, eps)

   END SUBROUTINE

END MODULE

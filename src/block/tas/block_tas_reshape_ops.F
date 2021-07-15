!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE block_tas_reshape_ops
   !! communication routines to reshape / replicate / merge tall-and-skinny matrices.

   USE block_matrix_api, ONLY: &
      block_matrix_put_block, block_matrix_reserve_blocks, &
      block_matrix_distribution_col_dist, block_matrix_distribution_row_dist, &
      block_matrix_get_stored_coordinates, &
      block_matrix_iterator_blocks_left, block_matrix_iterator_next_block, &
      block_matrix_iterator_start, block_matrix_iterator_stop, &
      block_matrix_row_block_sizes, block_matrix_col_block_sizes, &
      block_matrix_get_info, block_matrix_clear, &
      block_matrix_distribution_obj, block_matrix_iterator, block_matrix_type, &
      block_matrix_finalize, block_matrix_get_name
   USE block_tas_base, ONLY: &
      block_tas_blk_sizes, block_tas_create, block_tas_distribution_new, block_tas_finalize, &
      block_tas_get_stored_coordinates, block_tas_info, &
      block_tas_iterator_blocks_left, block_tas_iterator_next_block, block_tas_iterator_start, &
      block_tas_iterator_stop, block_tas_put_block, block_tas_reserve_blocks, &
      block_repl_get_stored_coordinates, block_tas_clear
   USE block_tas_types, ONLY: &
      block_tas_distribution_type, block_tas_iterator, block_tas_split_info, block_tas_type
   USE block_tas_global, ONLY: &
      block_tas_blk_size_arb, block_tas_blk_size_repl, block_tas_dist_arb, block_tas_dist_repl, &
      block_tas_distribution, block_tas_rowcol_data
   USE block_tas_split, ONLY: &
      colsplit, block_tas_get_split_info, rowsplit
   USE kinds, ONLY: &
      int_8, real_8, real_4, dp
   USE message_passing, ONLY: &
      mp_alltoall, mp_environ, mp_isend, mp_irecv, mp_waitall
   USE block_tas_util, ONLY: &
      swap, index_unique
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_tas_reshape_ops'

   PUBLIC :: &
      block_tas_merge, &
      block_tas_replicate, &
      block_tas_reshape

   TYPE block_buffer_type
      INTEGER :: nblock = -1
      INTEGER(KIND=int_8), DIMENSION(:, :), ALLOCATABLE :: indx
      REAL(dp), DIMENSION(:), ALLOCATABLE :: msg
      INTEGER :: endpos = -1
   END TYPE

CONTAINS

   RECURSIVE SUBROUTINE block_tas_reshape(matrix_in, matrix_out, summation, transposed, move_data)
      !! copy data (involves reshape)

      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix_in, matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation
         !! whether matrix_out = matrix_out + matrix_in
      LOGICAL, INTENT(IN), OPTIONAL                      :: transposed
      LOGICAL, INTENT(IN), OPTIONAL                      :: move_data
         !! memory optimization: move data to matrix_out such that matrix_in is empty on return

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_reshape'

      INTEGER                                            :: handle, handle2, iproc, mp_comm, mynode, nblock, &
                                                            ndata, numnodes, bcount, nblk
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: index_recv, blks_to_allocate
      INTEGER(KIND=int_8), DIMENSION(2)                  :: blk_index
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array
      INTEGER, DIMENSION(2)                              :: blk_size
      LOGICAL                                            :: tr, tr_in, move_prv
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(block_tas_iterator)                           :: iter
      TYPE(block_tas_split_info)                         :: info

      CALL timeset(routineN, handle)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL block_matrix_clear(matrix_out%matrix)
      ELSE
         CALL block_matrix_clear(matrix_out%matrix)
      ENDIF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      ENDIF

      IF (PRESENT(transposed)) THEN
         tr_in = transposed
      ELSE
         tr_in = .FALSE.
      ENDIF

      IF (.NOT. matrix_out%valid) THEN
         CPABORT("can not reshape into invalid matrix")
      ENDIF

      info = block_tas_info(matrix_in)
      mp_comm = info%mp_comm
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (buffer_send(0:numnodes - 1))
      ALLOCATE (buffer_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_send(0:numnodes - 1))
      ALLOCATE (num_entries_recv(0:numnodes - 1))
      ALLOCATE (num_entries_send(0:numnodes - 1))
      ALLOCATE (num_rec(0:2*numnodes - 1))
      ALLOCATE (num_send(0:2*numnodes - 1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))
      CALL block_tas_iterator_start(iter, matrix_in)

      CALL timeset(routineN//"_get_coord", handle2)
      DO WHILE (block_tas_iterator_blocks_left(iter))
         CALL block_tas_iterator_next_block(iter, blk_index(1), blk_index(2), nblock, transposed=tr, &
                                            row_size=blk_size(1), col_size=blk_size(2))

         IF (tr_in) THEN
            CALL block_tas_get_stored_coordinates(matrix_out, blk_index(2), blk_index(1), iproc)
         ELSE
            CALL block_tas_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         ENDIF

         num_send(2*iproc) = num_send(2*iproc) + PRODUCT(blk_size)
         num_send(2*iproc + 1) = num_send(2*iproc + 1) + 1
      ENDDO
      CALL block_tas_iterator_stop(iter)
      CALL timestop(handle2)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      CALL timeset(routineN//"_buffer_fill", handle2)
      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc))

         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc))

      ENDDO

      CALL block_tas_iterator_start(iter, matrix_in)
      DO WHILE (block_tas_iterator_blocks_left(iter))
         CALL block_tas_iterator_next_block(iter, blk_index(1), blk_index(2), block, tr, &
                                            row_size=blk_size(1), col_size=blk_size(2))

         CPASSERT(tr .EQV. .FALSE.)

         IF (tr_in) THEN
            CALL block_tas_get_stored_coordinates(matrix_out, blk_index(2), blk_index(1), iproc)
         ELSE
            CALL block_tas_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         ENDIF
         CALL block_buffer_add_block(buffer_send(iproc), blk_index, block, transposed=tr_in)
      ENDDO
      CALL block_tas_iterator_stop(iter)

      IF (move_prv) CALL block_tas_clear(matrix_in)

      CALL timestop(handle2)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL block_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

      DO iproc = 0, numnodes - 1
         CALL block_buffer_destroy(buffer_send(iproc))
      ENDDO

      CALL timestop(handle2)

      CALL timeset(routineN//"_buffer_obtain", handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = index_recv(:, :)
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      ENDDO

      CALL block_tas_reserve_blocks(matrix_out, blks_to_allocate(:, 1), blks_to_allocate(:, 2))
      DEALLOCATE (blks_to_allocate)

      DO iproc = 0, numnodes - 1
         ! First, we need to get the index to create block
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index)
            CALL block_tas_blk_sizes(matrix_out, blk_index(1), blk_index(2), blk_size(1), blk_size(2))
            ALLOCATE(block(blk_size(1), blk_size(2)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index, block)
            CALL block_tas_put_block(matrix_out, blk_index(1), blk_index(2), block, summation=summation)
            DEALLOCATE(block)
         ENDDO
         CALL block_buffer_destroy(buffer_recv(iproc))
      ENDDO

      CALL timestop(handle2)

      CALL block_tas_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_tas_replicate(matrix_in, info, matrix_out, nodata, move_data)
      !! Replicate matrix_in such that each submatrix of matrix_out is an exact copy of matrix_in

      TYPE(block_matrix_type), INTENT(INOUT)                    :: matrix_in
      TYPE(block_tas_split_info), INTENT(IN)             :: info
      TYPE(block_tas_type), INTENT(OUT)                  :: matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: nodata
         !! Don't copy data but create matrix_out
      LOGICAL, INTENT(IN), OPTIONAL                      :: move_data
         !! memory optimization: move data to matrix_out such that matrix_in is empty on return

      INTEGER                                            :: nblkcols, nblkrows
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, col_dist, &
                                                            row_blk_size, row_dist
      TYPE(block_matrix_distribution_obj)                       :: block_matrix_dist
      TYPE(block_tas_dist_arb), TARGET                     :: dir_dist
      TYPE(block_tas_dist_repl), TARGET                    :: repl_dist

      CLASS(block_tas_distribution), ALLOCATABLE :: col_dist_obj, row_dist_obj
      CLASS(block_tas_rowcol_data), ALLOCATABLE :: row_bsize_obj, col_bsize_obj
      TYPE(block_tas_blk_size_repl), TARGET :: repl_blksize
      TYPE(block_tas_blk_size_arb), TARGET :: dir_blksize
      TYPE(block_tas_distribution_type) :: dist
      INTEGER :: mp_comm, numnodes, mynode
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: num_blocks_recv, num_blocks_send, &
                                                            num_entries_recv, num_entries_send, &
                                                            num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: req_array, blks_to_allocate
      INTEGER, DIMENSION(2) :: blk_size
      INTEGER, DIMENSION(2) :: blk_index
      INTEGER(KIND=int_8), DIMENSION(2) :: blk_index_i8
      TYPE(block_matrix_iterator) :: iter
      INTEGER :: nblock, i, iproc, bcount, nblk
      INTEGER, DIMENSION(:), ALLOCATABLE :: iprocs
      LOGICAL :: tr, nodata_prv, move_prv
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :) :: index_recv
      INTEGER :: ndata

      REAL(KIND=dp), DIMENSION(:, :), POINTER :: block

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_replicate'

      INTEGER :: handle, handle2

      NULLIFY (col_blk_size, row_blk_size)

      CALL timeset(routineN, handle)

      IF (PRESENT(nodata)) THEN
         nodata_prv = nodata
      ELSE
         nodata_prv = .FALSE.
      ENDIF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      ENDIF

      CALL block_matrix_get_info(matrix_in, distribution=block_matrix_dist,  &
                          nblkrows_total=nblkrows, nblkcols_total=nblkcols, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size)
      row_dist => block_matrix_distribution_row_dist(block_matrix_dist)
      col_dist => block_matrix_distribution_col_dist(block_matrix_dist)

      mp_comm = info%mp_comm

      CALL mp_environ(numnodes, mynode, mp_comm)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm)

      SELECT CASE (info%split_rowcol)
      CASE (rowsplit)
         repl_dist = block_tas_dist_repl(row_dist, pdims(1), nblkrows, info%ngroup, info%pgrid_split_size)
         dir_dist = block_tas_dist_arb(col_dist, pdims(2), INT(nblkcols, KIND=int_8))
         repl_blksize = block_tas_blk_size_repl(row_blk_size, info%ngroup)
         dir_blksize = block_tas_blk_size_arb(col_blk_size)
         ALLOCATE (row_dist_obj, source=repl_dist)
         ALLOCATE (col_dist_obj, source=dir_dist)
         ALLOCATE (row_bsize_obj, source=repl_blksize)
         ALLOCATE (col_bsize_obj, source=dir_blksize)
      CASE (colsplit)
         dir_dist = block_tas_dist_arb(row_dist, pdims(1), INT(nblkrows, KIND=int_8))
         repl_dist = block_tas_dist_repl(col_dist, pdims(2), nblkcols, info%ngroup, info%pgrid_split_size)
         dir_blksize = block_tas_blk_size_arb(row_blk_size)
         repl_blksize = block_tas_blk_size_repl(col_blk_size, info%ngroup)
         ALLOCATE (row_dist_obj, source=dir_dist)
         ALLOCATE (col_dist_obj, source=repl_dist)
         ALLOCATE (row_bsize_obj, source=dir_blksize)
         ALLOCATE (col_bsize_obj, source=repl_blksize)
      END SELECT

      CALL block_tas_distribution_new(dist, mp_comm, row_dist_obj, col_dist_obj, split_info=info)
      CALL block_tas_create(matrix_out, TRIM(block_matrix_get_name(matrix_in))//" replicated", &
                            dist, row_bsize_obj, col_bsize_obj, own_dist=.TRUE.)

      IF (nodata_prv) THEN
         CALL block_tas_finalize(matrix_out)
         CALL timestop(handle)
         RETURN
      ENDIF

      ALLOCATE (buffer_send(0:numnodes - 1))
      ALLOCATE (buffer_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_send(0:numnodes - 1))
      ALLOCATE (num_entries_recv(0:numnodes - 1))
      ALLOCATE (num_entries_send(0:numnodes - 1))
      ALLOCATE (num_rec(0:2*numnodes - 1))
      ALLOCATE (num_send(0:2*numnodes - 1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))

      ALLOCATE (iprocs(info%ngroup))
      CALL block_matrix_iterator_start(iter, matrix_in)
      DO WHILE (block_matrix_iterator_blocks_left(iter))
         CALL block_matrix_iterator_next_block(iter, blk_index(1), blk_index(2), nblock, transposed=tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL block_repl_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iprocs)
         DO i = 1, SIZE(iprocs)
            num_send(2*iprocs(i)) = num_send(2*iprocs(i)) + PRODUCT(blk_size)
            num_send(2*iprocs(i) + 1) = num_send(2*iprocs(i) + 1) + 1
         ENDDO
      ENDDO
      CALL block_matrix_iterator_stop(iter)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc))

         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc))

      ENDDO

      CALL block_matrix_iterator_start(iter, matrix_in)
      DO WHILE (block_matrix_iterator_blocks_left(iter))
         CALL block_matrix_iterator_next_block(iter, blk_index(1), blk_index(2), block, tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL block_repl_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iprocs)
         DO i = 1, SIZE(iprocs)
            CALL block_buffer_add_block(buffer_send(iprocs(i)), INT(blk_index, KIND=int_8), block)
         ENDDO
      ENDDO
      CALL block_matrix_iterator_stop(iter)

      IF (move_prv) CALL block_matrix_clear(matrix_in)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL block_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

      DO iproc = 0, numnodes - 1
         CALL block_buffer_destroy(buffer_send(iproc))
      ENDDO

      CALL timestop(handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = INT(index_recv(:, :))
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      ENDDO

      CALL block_matrix_reserve_blocks(matrix_out%matrix, blks_to_allocate(:, 1), blks_to_allocate(:, 2))
      DEALLOCATE (blks_to_allocate)

      DO iproc = 0, numnodes - 1
         ! First, we need to get the index to create block
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8)
            CALL block_tas_blk_sizes(matrix_out, blk_index_i8(1), blk_index_i8(2), blk_size(1), blk_size(2))
            ALLOCATE(block(blk_size(1), blk_size(2)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8, block)
            CALL block_matrix_put_block(matrix_out%matrix, INT(blk_index_i8(1)), INT(blk_index_i8(2)), block)
            DEALLOCATE(block)
         ENDDO

         CALL block_buffer_destroy(buffer_recv(iproc))
      ENDDO

      CALL block_tas_finalize(matrix_out)

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE block_tas_merge(matrix_out, matrix_in, summation, move_data)
      !! Merge submatrices of matrix_in to matrix_out by sum

      TYPE(block_matrix_type), INTENT(INOUT)                    :: matrix_out
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix_in
      LOGICAL, INTENT(IN), OPTIONAL                      :: summation
      LOGICAL, INTENT(IN), OPTIONAL                      :: move_data
         !! memory optimization: move data to matrix_out such that matrix_in is empty on return

      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_merge'
      INTEGER                                            :: handle, handle2, iproc, mp_comm, &
                                                            mynode, nblock, ndata, numnodes, nblk, bcount
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:, :)  :: index_recv
      INTEGER(KIND=int_8), DIMENSION(2)                  :: blk_index_i8
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iprocs, num_blocks_recv, &
                                                            num_blocks_send, num_entries_recv, &
                                                            num_entries_send, num_rec, num_send
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: req_array, blks_to_allocate, blks_to_allocate_u
      INTEGER, DIMENSION(2)                              :: blk_index, blk_size
      LOGICAL                                            :: tr, move_prv
      TYPE(block_buffer_type), ALLOCATABLE, DIMENSION(:) :: buffer_recv, buffer_send
      REAL(dp), DIMENSION(:, :), POINTER                 :: block
      TYPE(block_matrix_iterator)                               :: iter
      TYPE(block_tas_split_info)                         :: info
      INTEGER, DIMENSION(:), POINTER                     :: row_block_sizes, col_block_sizes

      CALL timeset(routineN, handle)

      IF (PRESENT(summation)) THEN
         IF (.NOT. summation) CALL block_matrix_clear(matrix_out)
      ELSE
         CALL block_matrix_clear(matrix_out)
      ENDIF

      IF (PRESENT(move_data)) THEN
         move_prv = move_data
      ELSE
         move_prv = .FALSE.
      ENDIF

      info = block_tas_info(matrix_in)
      CALL block_tas_get_split_info(info, mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)

      ALLOCATE (buffer_send(0:numnodes - 1))
      ALLOCATE (buffer_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_recv(0:numnodes - 1))
      ALLOCATE (num_blocks_send(0:numnodes - 1))
      ALLOCATE (num_entries_recv(0:numnodes - 1))
      ALLOCATE (num_entries_send(0:numnodes - 1))
      ALLOCATE (num_rec(0:2*numnodes - 1))
      ALLOCATE (num_send(0:2*numnodes - 1))
      num_send(:) = 0
      ALLOCATE (req_array(1:numnodes, 4))

      ALLOCATE (iprocs(info%ngroup))

      CALL block_matrix_iterator_start(iter, matrix_in%matrix)
      DO WHILE (block_matrix_iterator_blocks_left(iter))
         CALL block_matrix_iterator_next_block(iter, blk_index(1), blk_index(2), nblock, transposed=tr, &
                                        row_size=blk_size(1), col_size=blk_size(2))
         CALL block_matrix_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         num_send(2*iproc) = num_send(2*iproc) + PRODUCT(blk_size)
         num_send(2*iproc + 1) = num_send(2*iproc + 1) + 1
      ENDDO
      CALL block_matrix_iterator_stop(iter)

      CALL timeset(routineN//"_alltoall", handle2)
      CALL mp_alltoall(num_send, num_rec, 2, mp_comm)
      CALL timestop(handle2)

      DO iproc = 0, numnodes - 1
         num_entries_recv(iproc) = num_rec(2*iproc)
         num_blocks_recv(iproc) = num_rec(2*iproc + 1)
         num_entries_send(iproc) = num_send(2*iproc)
         num_blocks_send(iproc) = num_send(2*iproc + 1)

         CALL block_buffer_create(buffer_send(iproc), num_blocks_send(iproc), num_entries_send(iproc))

         CALL block_buffer_create(buffer_recv(iproc), num_blocks_recv(iproc), num_entries_recv(iproc))

      ENDDO

      CALL block_matrix_iterator_start(iter, matrix_in%matrix)
      DO WHILE (block_matrix_iterator_blocks_left(iter))
         CALL block_matrix_iterator_next_block(iter, blk_index(1), blk_index(2), block, tr, &
                                               row_size=blk_size(1), col_size=blk_size(2))
         CALL block_matrix_get_stored_coordinates(matrix_out, blk_index(1), blk_index(2), iproc)
         CALL block_buffer_add_block(buffer_send(iproc), INT(blk_index, KIND=int_8), block)
      ENDDO

      CALL block_matrix_iterator_stop(iter)

      IF (move_prv) CALL block_tas_clear(matrix_in)

      CALL timeset(routineN//"_communicate_buffer", handle2)
      CALL block_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)

      DO iproc = 0, numnodes - 1
         CALL block_buffer_destroy(buffer_send(iproc))
      ENDDO

      CALL timestop(handle2)

      nblk = SUM(num_blocks_recv)
      ALLOCATE (blks_to_allocate(nblk, 2))

      bcount = 0
      DO iproc = 0, numnodes - 1
         CALL block_buffer_get_index(buffer_recv(iproc), index_recv)
         blks_to_allocate(bcount + 1:bcount + SIZE(index_recv, 1), :) = INT(index_recv(:, :))
         bcount = bcount + SIZE(index_recv, 1)
         DEALLOCATE (index_recv)
      ENDDO

      CALL index_unique(blks_to_allocate, blks_to_allocate_u)

      CALL block_matrix_reserve_blocks(matrix_out, blks_to_allocate_u(:, 1), blks_to_allocate_u(:, 2))
      DEALLOCATE (blks_to_allocate, blks_to_allocate_u)

      DO iproc = 0, numnodes - 1
         ! First, we need to get the index to create block
         DO WHILE (block_buffer_blocks_left(buffer_recv(iproc)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8)
            row_block_sizes => block_matrix_row_block_sizes(matrix_out)
            col_block_sizes => block_matrix_col_block_sizes(matrix_out)
            blk_size(1) = row_block_sizes(INT(blk_index_i8(1)))
            blk_size(2) = col_block_sizes(INT(blk_index_i8(2)))
            ALLOCATE(block(blk_size(1), blk_size(2)))
            CALL block_buffer_get_next_block(buffer_recv(iproc), ndata, blk_index_i8, block)
            CALL block_matrix_put_block(matrix_out, INT(blk_index_i8(1)), INT(blk_index_i8(2)), block, summation=.TRUE.)
            DEALLOCATE(block)
         ENDDO
         CALL block_buffer_destroy(buffer_recv(iproc))
      ENDDO

      CALL block_matrix_finalize(matrix_out)

      CALL timestop(handle)
   END SUBROUTINE

   SUBROUTINE block_buffer_get_index(buffer, index)
      !! get all indices from buffer
      TYPE(block_buffer_type), INTENT(IN)                :: buffer
      INTEGER(KIND=int_8), ALLOCATABLE, &
         DIMENSION(:, :), INTENT(OUT)                    :: index

      INTEGER, DIMENSION(2)                              :: indx_shape
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_buffer_get_index'
      INTEGER :: handle

      CALL timeset(routineN, handle)

      indx_shape = SHAPE(buffer%indx) - [0, 1]
      ALLOCATE (INDEX(indx_shape(1), indx_shape(2)))
      INDEX(:, :) = buffer%indx(1:indx_shape(1), 1:indx_shape(2))
      CALL timestop(handle)
   END SUBROUTINE

   PURE FUNCTION block_buffer_blocks_left(buffer)
      !! how many blocks left in iterator
      TYPE(block_buffer_type), INTENT(IN)                :: buffer
      LOGICAL                                            :: block_buffer_blocks_left

      block_buffer_blocks_left = buffer%endpos .LT. buffer%nblock
   END FUNCTION

   SUBROUTINE block_buffer_create(buffer, nblock, ndata)
      !! Create block buffer for MPI communication.

      TYPE(block_buffer_type), INTENT(OUT) :: buffer
         !! block buffer
      INTEGER, INTENT(IN) :: nblock, ndata
         !! number of blocks
         !! total number of block entries

      buffer%nblock = nblock
      buffer%endpos = 0
      ALLOCATE (buffer%msg(ndata))
      ALLOCATE (buffer%indx(nblock, 3))
   END SUBROUTINE

   SUBROUTINE block_buffer_destroy(buffer)
      TYPE(block_buffer_type), INTENT(INOUT) :: buffer

      DEALLOCATE (buffer%msg)
      DEALLOCATE (buffer%indx)
      buffer%nblock = -1
      buffer%endpos = -1
   END SUBROUTINE block_buffer_destroy

   SUBROUTINE block_buffer_add_block(buffer, index, block, transposed)
      !! insert a block into block buffer (at current iterator position)

      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER(KIND=int_8), DIMENSION(2), &
         INTENT(IN)                               :: index
         !! index of block
      REAL(dp), DIMENSION(:, :), INTENT(IN)       :: block
         !! block
      LOGICAL, INTENT(IN), OPTIONAL :: transposed
      INTEGER(KIND=int_8), DIMENSION(2)   :: index_prv

      LOGICAL :: tr
      INTEGER :: p
      INTEGER :: ndata
      INTEGER :: p_data

      IF (PRESENT(transposed)) THEN
         tr = transposed
      ELSE
         tr = .FALSE.
      ENDIF

      index_prv(:) = index(:)
      IF (tr) THEN
         CALL swap(index_prv)
      ENDIF
      ndata = PRODUCT(SHAPE(block))

      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = INT(buffer%indx(p, 3))
      ENDIF

      IF (tr) THEN
         buffer%msg (p_data + 1:p_data + ndata) = RESHAPE(TRANSPOSE(block), [ndata])
      ELSE
         buffer%msg (p_data + 1:p_data + ndata) = RESHAPE(block, [ndata])
      ENDIF

      buffer%indx(p + 1, 1:2) = index_prv(:)
      IF (p > 0) THEN
         buffer%indx(p + 1, 3) = buffer%indx(p, 3) + INT(ndata, KIND=int_8)
      ELSE
         buffer%indx(p + 1, 3) = INT(ndata, KIND=int_8)
      ENDIF
      buffer%endpos = buffer%endpos + 1
   END SUBROUTINE

   SUBROUTINE block_buffer_get_next_block(buffer, ndata, index, block, advance_iter)
      !! get next block from buffer. Iterator is advanced only if block is retrieved or advance_iter.
      TYPE(block_buffer_type), INTENT(INOUT)      :: buffer
      INTEGER, INTENT(OUT)                        :: ndata
      INTEGER(KIND=int_8), DIMENSION(2), &
         INTENT(OUT)                              :: index
      REAL(dp), DIMENSION(:, :), INTENT(OUT), OPTIONAL     :: block
      LOGICAL, INTENT(IN), OPTIONAL               :: advance_iter
      INTEGER :: p, p_data
      LOGICAL :: do_advance

      do_advance = .FALSE.
      IF (PRESENT(advance_iter)) THEN
         do_advance = advance_iter
      ELSE IF (PRESENT(block)) THEN
         do_advance = .TRUE.
      ENDIF

      p = buffer%endpos
      IF (p .EQ. 0) THEN
         p_data = 0
      ELSE
         p_data = INT(buffer%indx(p, 3))
      ENDIF

      IF (p > 0) THEN
         ndata = INT(buffer%indx(p + 1, 3) - buffer%indx(p, 3))
      ELSE
         ndata = INT(buffer%indx(p + 1, 3))
      ENDIF
      index(:) = buffer%indx(p + 1, 1:2)

      IF (PRESENT(block)) THEN
         block(:, :) = RESHAPE(buffer%msg(p_data + 1:p_data + ndata), SHAPE(block))
      ENDIF

      IF (do_advance) buffer%endpos = buffer%endpos + 1
   END SUBROUTINE

   SUBROUTINE block_tas_communicate_buffer(mp_comm, buffer_recv, buffer_send, req_array)
      !! communicate buffer
      INTEGER, INTENT(IN)                    :: mp_comm
      TYPE(block_buffer_type), DIMENSION(0:), INTENT(INOUT) :: buffer_recv, buffer_send
      INTEGER, DIMENSION(:, :), INTENT(OUT)               :: req_array

      INTEGER                                :: iproc, mynode, numnodes, rec_counter, &
                                                send_counter
      INTEGER                                   :: handle
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_communicate_buffer'

      CALL timeset(routineN, handle)
      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (numnodes > 1) THEN

         send_counter = 0
         rec_counter = 0

         DO iproc = 0, numnodes - 1
            IF (buffer_recv(iproc)%nblock > 0) THEN
               rec_counter = rec_counter + 1
               CALL mp_irecv(buffer_recv(iproc)%indx, iproc, mp_comm, req_array(rec_counter, 3), tag=4)
               CALL mp_irecv(buffer_recv(iproc)%msg, iproc, mp_comm, req_array(rec_counter, 4), tag=7)
            END IF
         END DO

         DO iproc = 0, numnodes - 1
            IF (buffer_send(iproc)%nblock > 0) THEN
               send_counter = send_counter + 1
               CALL mp_isend(buffer_send(iproc)%indx, iproc, mp_comm, req_array(send_counter, 1), tag=4)
               CALL mp_isend(buffer_send(iproc)%msg, iproc, mp_comm, req_array(send_counter, 2), tag=7)
            END IF
         END DO

         IF (send_counter > 0) THEN
            CALL mp_waitall(req_array(1:send_counter, 1:2))
         ENDIF
         IF (rec_counter > 0) THEN
            CALL mp_waitall(req_array(1:rec_counter, 3:4))
         ENDIF

      ELSE
         IF (buffer_recv(0)%nblock > 0) THEN
            buffer_recv(0)%indx(:, :) = buffer_send(0)%indx(:, :)
            buffer_recv(0)%msg(:) = buffer_send(0)%msg(:)
         ENDIF
      ENDIF
      CALL timestop(handle)
   END SUBROUTINE

END MODULE

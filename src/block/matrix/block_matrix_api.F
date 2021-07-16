!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE block_matrix_api

 USE dbcsr_types, ONLY: &
     dbcsr_type

 USE dbcsr_iterator_operations, ONLY: &
    dbcsr_iterator_blocks_left, &
    dbcsr_iterator_next_block, &
    dbcsr_iterator_start, &
    dbcsr_iterator_stop

 USE dbcsr_transformations, ONLY: &
    dbcsr_redistribute, &
    dbcsr_complete_redistribute

 USE dbcsr_operations, ONLY: &
      dbcsr_set, &
      dbcsr_filter, &
      dbcsr_scale, &
      dbcsr_get_info, &
      dbcsr_clear, &
      dbcsr_copy, &
      dbcsr_add, &
      dbcsr_maxabs, &
      dbcsr_zero

 USE dbcsr_block_access, ONLY: &
      dbcsr_get_block_p, &
      dbcsr_reserve_blocks, &
      dbcsr_put_block

 USE dbcsr_methods, ONLY: &
      dbcsr_mp_release, &
      dbcsr_get_num_blocks, &
      dbcsr_get_nze, &
      dbcsr_row_block_sizes, &
      dbcsr_col_block_sizes, &
      dbcsr_release

 USE dbcsr_dist_methods, ONLY: &
      dbcsr_distribution_col_dist, &
      dbcsr_distribution_new, &
      dbcsr_distribution_row_dist, &
      dbcsr_distribution_hold, &
      dbcsr_distribution_release

 USE dbcsr_dist_operations, ONLY: &
      dbcsr_get_stored_coordinates

 USE dbcsr_types, ONLY: &
      block_matrix_transpose => dbcsr_transpose, &
      block_matrix_no_transpose => dbcsr_no_transpose, &
      dbcsr_mp_obj, &
      dbcsr_type_real_8, &
      dbcsr_type_no_symmetry, &
      dbcsr_iterator, &
      dbcsr_distribution_obj

 USE dbcsr_multiply_api, ONLY: dbcsr_multiply

 USE dbcsr_work_operations, ONLY: &
      dbcsr_create, &
      dbcsr_finalize

 USE kinds, ONLY: &
      default_string_length, &
      int_8, &
      dp

 USE dbcsr_mp_methods, ONLY: &
      dbcsr_mp_new

 USE dbcsr_config, ONLY: block_matrix_cfg => dbcsr_cfg

 USE dbcsr_dist_util, ONLY: &
      dbcsr_checksum

 USE message_passing, ONLY: mp_cart_rank, mp_environ


#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_matrix_api'

   PUBLIC :: block_matrix_type
   PUBLIC :: block_matrix_release
   PUBLIC :: block_matrix_create

   PUBLIC :: block_matrix_iterator
   PUBLIC :: block_matrix_iterator_start
   PUBLIC :: block_matrix_iterator_stop
   PUBLIC :: block_matrix_iterator_blocks_left
   PUBLIC :: block_matrix_iterator_next_block

   PUBLIC :: block_matrix_clear
   PUBLIC :: block_matrix_scale
   PUBLIC :: block_matrix_get_block_p
   PUBLIC :: block_matrix_put_block
   PUBLIC :: block_matrix_reserve_blocks
   PUBLIC :: block_matrix_set
   PUBLIC :: block_matrix_filter
   PUBLIC :: block_matrix_finalize
   PUBLIC :: block_matrix_multiply
   PUBLIC :: block_matrix_redistribute
   PUBLIC :: block_matrix_copy
   PUBLIC :: block_matrix_add
   PUBLIC :: block_matrix_maxabs
   PUBLIC :: block_matrix_zero
   PUBLIC :: block_matrix_checksum
   PUBLIC :: block_matrix_transpose
   PUBLIC :: block_matrix_no_transpose

   PUBLIC :: block_matrix_distribution_obj
   PUBLIC :: block_matrix_distribution_col_dist
   PUBLIC :: block_matrix_distribution_new
   PUBLIC :: block_matrix_distribution_row_dist
   PUBLIC :: block_matrix_distribution_hold
   PUBLIC :: block_matrix_distribution_release

   PUBLIC :: block_matrix_get_name
   PUBLIC :: block_matrix_get_num_blocks
   PUBLIC :: block_matrix_get_nze
   PUBLIC :: block_matrix_get_stored_coordinates
   PUBLIC :: block_matrix_row_block_sizes
   PUBLIC :: block_matrix_col_block_sizes

   PUBLIC :: generate_larnv_seed

   ! TODO: remove these
   PUBLIC :: block_matrix_complete_redistribute
   PUBLIC :: block_matrix_get_info
   PUBLIC :: block_matrix_cfg

   TYPE block_matrix_distribution_obj
      TYPE(dbcsr_distribution_obj), PRIVATE        :: dbcsr
   END TYPE block_matrix_distribution_obj

   TYPE block_matrix_type
      TYPE(dbcsr_type), PRIVATE                    :: dbcsr
   END TYPE block_matrix_type

   TYPE block_matrix_iterator
      TYPE(dbcsr_iterator), PRIVATE        :: dbcsr
   END TYPE block_matrix_iterator

   INTERFACE block_matrix_create
      MODULE PROCEDURE block_matrix_create_new
      MODULE PROCEDURE block_matrix_create_template
   END INTERFACE block_matrix_create

   INTERFACE block_matrix_iterator_next_block
      MODULE PROCEDURE block_matrix_iterator_next_block_d
      MODULE PROCEDURE block_matrix_iterator_next_block_index
   END INTERFACE block_matrix_iterator_next_block

   CONTAINS

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_create_new(matrix, name, dist, &
                                      row_blk_size, col_blk_size)

      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix
      CHARACTER(len=*), INTENT(IN)                       :: name
      TYPE(block_matrix_distribution_obj), INTENT(IN)    :: dist
      INTEGER, DIMENSION(:), INTENT(IN), POINTER, CONTIGUOUS :: row_blk_size, col_blk_size

      CALL dbcsr_create(matrix%dbcsr, name=name, dist=dist%dbcsr, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                        data_type=dbcsr_type_real_8)
   END SUBROUTINE block_matrix_create_new

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_create_template(matrix, template, name)

      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix
      TYPE(block_matrix_type), INTENT(IN)                :: template
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name

      CALL dbcsr_create(matrix%dbcsr, template=template%dbcsr, name=name, &
                        matrix_type=dbcsr_type_no_symmetry, &
                        data_type=dbcsr_type_real_8)
   END SUBROUTINE block_matrix_create_template

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_iterator_blocks_left(iterator) RESULT(blocks_left)
      TYPE(block_matrix_iterator), INTENT(IN)            :: iterator
      LOGICAL                                            :: blocks_left

      blocks_left = dbcsr_iterator_blocks_left(iterator%dbcsr)
   END FUNCTION block_matrix_iterator_blocks_left

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_iterator_stop(iterator)
      TYPE(block_matrix_iterator), INTENT(INOUT)           :: iterator

      CALL dbcsr_iterator_stop(iterator%dbcsr)
   END SUBROUTINE block_matrix_iterator_stop

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_iterator_start(iterator, matrix)
      TYPE(block_matrix_iterator), INTENT(OUT)           :: iterator
      TYPE(block_matrix_type), INTENT(IN)                :: matrix

      CALL dbcsr_iterator_start(iterator%dbcsr, matrix%dbcsr)
   END SUBROUTINE block_matrix_iterator_start

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_iterator_next_block_d(iterator, row, column, block, row_size, col_size)
      TYPE(block_matrix_iterator), INTENT(INOUT) :: iterator
      INTEGER, INTENT(OUT)                       :: row, column
      REAL(dp), DIMENSION(:, :), POINTER         :: block
      INTEGER, INTENT(OUT), OPTIONAL             :: row_size, col_size

      LOGICAL :: transposed
      CALL dbcsr_iterator_next_block(iterator%dbcsr, row, column, block, transposed, &
                                     row_size=row_size, col_size=col_size)
      CPASSERT(.NOT.transposed)
   END SUBROUTINE block_matrix_iterator_next_block_d

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_iterator_next_block_index(iterator, row, column, blk_p, row_size, col_size)
      TYPE(block_matrix_iterator), INTENT(INOUT)         :: iterator
      INTEGER, INTENT(OUT)                               :: row, column
      INTEGER, INTENT(OUT), OPTIONAL                     :: blk_p, row_size, col_size

      INTEGER :: blk ! only needed to select right procedure from interface
      LOGICAL :: transposed

      CALL dbcsr_iterator_next_block(iterator%dbcsr, row=row, column=column, blk=blk,&
                                     transposed=transposed, blk_p=blk_p, &
                                     row_size=row_size, col_size=col_size)
      CPASSERT(.NOT.transposed)
   END SUBROUTINE block_matrix_iterator_next_block_index

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_filter(matrix, eps)
      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix
      REAL(dp), INTENT(IN)                               :: eps

      CALL dbcsr_filter(matrix%dbcsr, eps)
   END SUBROUTINE block_matrix_filter

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_set(matrix, alpha)
      TYPE(block_matrix_type), INTENT(INOUT)     :: matrix
      REAL(dp), INTENT(IN)                       :: alpha

      CALL dbcsr_set(matrix%dbcsr, alpha)
   END SUBROUTINE block_matrix_set

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_get_block_p(matrix, row, col, block, tr, found, row_size, col_size)
      TYPE(block_matrix_type), INTENT(INOUT)   :: matrix
      INTEGER, INTENT(IN)                      :: row, col
      REAL(dp), DIMENSION(:, :), POINTER       :: block
      LOGICAL, INTENT(OUT)                     :: tr
      LOGICAL, INTENT(OUT)                     :: found
      INTEGER, INTENT(OUT), OPTIONAL           :: row_size, col_size

      CALL dbcsr_get_block_p(matrix%dbcsr, row, col, block, tr, found, row_size, col_size)
   END SUBROUTINE block_matrix_get_block_p

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_put_block(matrix, row, col, block, transposed, summation, scale)
      TYPE(block_matrix_type), INTENT(INOUT)   :: matrix
      INTEGER, INTENT(IN)                      :: row, col
      REAL(dp), DIMENSION(:, :), INTENT(IN)    :: block
      LOGICAL, INTENT(IN), OPTIONAL            :: transposed, summation
      REAL(dp), INTENT(IN), OPTIONAL           :: scale

      CALL dbcsr_put_block(matrix%dbcsr, row, col, block, &
                           transposed=transposed, summation=summation, scale=scale)
   END SUBROUTINE block_matrix_put_block

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_get_info(matrix, nblkrows_total, nblkcols_total, &
                                    nfullrows_total, nfullcols_total, &
                                    nblkrows_local, nblkcols_local, &
                                    nfullrows_local, nfullcols_local, &
                                    my_prow, my_pcol, &
                                    local_rows, local_cols, proc_row_dist, proc_col_dist, &
                                    row_blk_size, col_blk_size, row_blk_offset, col_blk_offset, &
                                    distribution, name, &
                                    group)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      INTEGER, INTENT(OUT), OPTIONAL :: nblkrows_total, nblkcols_total, nfullrows_total, &
                                        nfullcols_total, nblkrows_local, nblkcols_local, nfullrows_local, nfullcols_local, &
                                        my_prow, my_pcol
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: local_rows, local_cols, proc_row_dist, &
                                                            proc_col_dist, &
                                                            row_blk_size, col_blk_size, &
                                                            row_blk_offset, col_blk_offset
      TYPE(block_matrix_distribution_obj), INTENT(OUT), &
         OPTIONAL                                        :: distribution
      CHARACTER(len=*), INTENT(OUT), OPTIONAL            :: name
      INTEGER, INTENT(OUT), OPTIONAL                     :: group

      TYPE(dbcsr_distribution_obj)                       :: dist

      CALL dbcsr_get_info(matrix=matrix%dbcsr, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nfullrows_total, &
                          nfullcols_total=nfullcols_total, &
                          nblkrows_local=nblkrows_local, &
                          nblkcols_local=nblkcols_local, &
                          nfullrows_local=nfullrows_local, &
                          nfullcols_local=nfullcols_local, &
                          my_prow=my_prow, &
                          my_pcol=my_pcol, &
                          local_rows=local_rows, &
                          local_cols=local_cols, &
                          proc_row_dist=proc_row_dist, &
                          proc_col_dist=proc_col_dist, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          row_blk_offset=row_blk_offset, &
                          col_blk_offset=col_blk_offset, &
                          distribution=dist, &
                          name=name, &
                          group=group)

      IF (PRESENT(distribution)) THEN
         distribution%dbcsr = dist
      ENDIF

   END SUBROUTINE block_matrix_get_info

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_clear(matrix)
      TYPE(block_matrix_type), INTENT(INOUT) :: matrix

      CALL dbcsr_clear(matrix%dbcsr)
   END SUBROUTINE block_matrix_clear

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_get_name(matrix) RESULT(res)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      CHARACTER(len=default_string_length)               :: res

      CALL dbcsr_get_info(matrix%dbcsr, name=res)
   END FUNCTION block_matrix_get_name

   !================================================================================================
   !================================================================================================
   PURE FUNCTION block_matrix_get_nze(matrix) RESULT(res)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      INTEGER                                            :: res

      res = dbcsr_get_nze(matrix%dbcsr)
   END FUNCTION block_matrix_get_nze

   !================================================================================================
   !================================================================================================
   PURE FUNCTION block_matrix_get_num_blocks(matrix) RESULT(num_blocks)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      INTEGER                                            :: num_blocks

      num_blocks = dbcsr_get_num_blocks(matrix%dbcsr)
   END FUNCTION block_matrix_get_num_blocks

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_row_block_sizes(matrix) RESULT(res)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: res

      res => dbcsr_row_block_sizes(matrix%dbcsr)
   END FUNCTION block_matrix_row_block_sizes

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_col_block_sizes(matrix) RESULT(res)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: res

      res => dbcsr_col_block_sizes(matrix%dbcsr)
   END FUNCTION block_matrix_col_block_sizes

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_distribution_hold(dist)
      TYPE(block_matrix_distribution_obj)                :: dist

      CALL dbcsr_distribution_hold(dist%dbcsr)
   END SUBROUTINE block_matrix_distribution_hold

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_reserve_blocks(matrix, rows, cols)
      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix
      INTEGER, DIMENSION(:), INTENT(IN)                  :: rows, cols

      CALL dbcsr_reserve_blocks(matrix%dbcsr, rows, cols)
   END SUBROUTINE block_matrix_reserve_blocks

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_distribution_row_dist(dist) RESULT(res)
      TYPE(block_matrix_distribution_obj), INTENT(IN)    :: dist
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: res

      res => dbcsr_distribution_row_dist(dist%dbcsr)
   END FUNCTION block_matrix_distribution_row_dist

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_distribution_col_dist(dist) RESULT(res)
      TYPE(block_matrix_distribution_obj), INTENT(IN)    :: dist
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: res

      res => dbcsr_distribution_col_dist(dist%dbcsr)
   END FUNCTION block_matrix_distribution_col_dist

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_distribution_new(dist, mp_comm, row_dist_block, col_dist_block)
      TYPE(block_matrix_distribution_obj), INTENT(OUT)       :: dist
      INTEGER, INTENT(IN)                                    :: mp_comm
      INTEGER, DIMENSION(:), INTENT(IN), POINTER, CONTIGUOUS :: row_dist_block, col_dist_block

      TYPE(dbcsr_mp_obj)                                 :: mp_env
      INTEGER                                            :: mynode, numnodes, pcol, prow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: pgrid
      INTEGER, DIMENSION(2)                              :: coord, mycoord, npdims

      ! Create a dbcsr mp environment from communicator
      CALL mp_environ(numnodes, npdims, mycoord, mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      ALLOCATE (pgrid(0:npdims(1) - 1, 0:npdims(2) - 1))
      DO prow = 0, npdims(1) - 1
         DO pcol = 0, npdims(2) - 1
            coord = (/prow, pcol/)
            CALL mp_cart_rank(mp_comm, coord, pgrid(prow, pcol))
         ENDDO
      ENDDO
      CPASSERT(mynode == pgrid(mycoord(1), mycoord(2)))

      CALL dbcsr_mp_new(mp_env, mp_comm, pgrid, mynode, numnodes, mycoord(1), mycoord(2))
      CALL dbcsr_distribution_new(dist=dist%dbcsr, mp_env=mp_env, &
                                  row_dist_block=row_dist_block, col_dist_block=col_dist_block)
      CALL dbcsr_mp_release(mp_env)
   END SUBROUTINE block_matrix_distribution_new

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_finalize(matrix)
      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix

      CALL dbcsr_finalize(matrix%dbcsr)
   END SUBROUTINE block_matrix_finalize

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_distribution_release(dist)
      TYPE(block_matrix_distribution_obj)                :: dist

      CALL dbcsr_distribution_release(dist%dbcsr)
   END SUBROUTINE block_matrix_distribution_release

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_release(matrix)
      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix

      CALL dbcsr_release(matrix%dbcsr)
   END SUBROUTINE block_matrix_release

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_get_stored_coordinates(matrix, row, column, processor)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      INTEGER, INTENT(IN)                                :: row, column
      INTEGER, INTENT(OUT), OPTIONAL                     :: processor

      CALL dbcsr_get_stored_coordinates(matrix%dbcsr, row, column, processor)
   END SUBROUTINE block_matrix_get_stored_coordinates

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_scale(matrix, alpha)
      TYPE(block_matrix_type), INTENT(INOUT)    :: matrix
       REAL(dp), INTENT(IN)                     :: alpha

      CALL dbcsr_scale(matrix%dbcsr, alpha)
   END SUBROUTINE block_matrix_scale

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_copy(matrix_b, matrix_a, name)
      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix_b
      TYPE(block_matrix_type), INTENT(IN)                :: matrix_a
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: name

      CALL dbcsr_copy(matrix_b%dbcsr, matrix_a%dbcsr, name)

   END SUBROUTINE block_matrix_copy

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_redistribute(matrix, redist)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      TYPE(block_matrix_type), INTENT(INOUT)             :: redist

      CALL dbcsr_redistribute(matrix%dbcsr, redist%dbcsr)
   END SUBROUTINE block_matrix_redistribute

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_complete_redistribute(matrix, redist)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      TYPE(block_matrix_type), INTENT(INOUT)             :: redist

      CALL dbcsr_redistribute(matrix%dbcsr, redist%dbcsr)
      CALL dbcsr_complete_redistribute(matrix%dbcsr, redist%dbcsr)
   END SUBROUTINE block_matrix_complete_redistribute

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_zero(matrix)
      TYPE(block_matrix_type), INTENT(INOUT)             :: matrix

      CALL dbcsr_zero(matrix%dbcsr)
   END SUBROUTINE block_matrix_zero

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_multiply(transa, transb, &
                                      alpha, matrix_a, matrix_b, beta, matrix_c, &
                                      first_row, last_row, &
                                      first_column, last_column, first_k, last_k, &
                                      retain_sparsity, filter_eps, flop)
      CHARACTER(LEN=1), INTENT(IN)                   :: transa, transb
      REAL(kind=dp), INTENT(IN)                      :: alpha
      TYPE(block_matrix_type), INTENT(IN)            :: matrix_a, matrix_b
      REAL(kind=dp), INTENT(IN)                      :: beta
      TYPE(block_matrix_type), INTENT(INOUT)         :: matrix_c
      INTEGER, INTENT(IN), OPTIONAL                  :: first_row, last_row, &
                                                        first_column, last_column, &
                                                        first_k, last_k
      LOGICAL, INTENT(IN), OPTIONAL                  :: retain_sparsity
      REAL(kind=dp), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER(int_8), INTENT(OUT), OPTIONAL          :: flop

      CALL dbcsr_multiply(transa=transa, transb=transb, &
                          alpha=alpha, matrix_a=matrix_a%dbcsr, &
                          matrix_b=matrix_b%dbcsr, beta=beta, matrix_c=matrix_c%dbcsr, &
                          first_row=first_row, last_row=last_row, &
                          first_column=first_column, last_column=last_column, &
                          first_k=first_k, last_k=last_k, &
                          retain_sparsity=retain_sparsity, filter_eps=filter_eps, flop=flop)
   END SUBROUTINE block_matrix_multiply

   !================================================================================================
   !================================================================================================
   SUBROUTINE block_matrix_add(matrix_a, matrix_b, alpha, beta)
      TYPE(block_matrix_type), INTENT(INOUT)   :: matrix_a
      TYPE(block_matrix_type), INTENT(IN)      :: matrix_b
      REAL(kind=dp), INTENT(IN)                 :: alpha, beta

      CALL dbcsr_add(matrix_a%dbcsr, matrix_b%dbcsr, alpha_scalar=alpha, beta_scalar=beta)

   END SUBROUTINE block_matrix_add

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_checksum(matrix, local, pos) RESULT(checksum)
      TYPE(block_matrix_type), INTENT(IN)                :: matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: local, pos
      REAL(KIND=dp)                                      :: checksum

      checksum = dbcsr_checksum(matrix%dbcsr, local=local, pos=pos)
   END FUNCTION block_matrix_checksum

   !================================================================================================
   !================================================================================================
   FUNCTION block_matrix_maxabs(matrix) RESULT(norm)
      TYPE(block_matrix_type), INTENT(INOUT)         :: matrix
      REAL(KIND=dp)                                  :: norm

      norm = dbcsr_maxabs(matrix%dbcsr)
   END FUNCTION block_matrix_maxabs

   !================================================================================================
   !TODO: Move to some utility module.
   !================================================================================================
   FUNCTION generate_larnv_seed(irow, nrow, icol, ncol, ival) RESULT(iseed)
      !! generate a seed respecting the lapack constraints,
      !! - values between 0..4095 (2**12-1)
      !! - iseed(4) odd
      !! also try to avoid iseed that are zero.
      !! Have but with a unique 2D mapping (irow,icol), and a 'mini-seed' ival

      INTEGER, INTENT(IN)                                :: irow, nrow, icol, ncol, ival
         !! 1..nrow
         !! nrow
         !! 1..ncol
         !! ncol
         !! mini-seed
      INTEGER                                            :: iseed(4)
         !! a lapack compatible seed

      INTEGER(KIND=int_8)                                :: map

      map = ((irow - 1 + icol*INT(nrow, int_8))*(1 + MODULO(ival, 2**16)))*2 + 1 + 0*ncol ! ncol used
      iseed(4) = INT(MODULO(map, 2_int_8**12)); map = map/2_int_8**12; ! keep odd
      iseed(3) = INT(MODULO(IEOR(map, 3541_int_8), 2_int_8**12)); map = map/2_int_8**12
      iseed(2) = INT(MODULO(IEOR(map, 1153_int_8), 2_int_8**12)); map = map/2_int_8**12
      iseed(1) = INT(MODULO(IEOR(map, 2029_int_8), 2_int_8**12)); map = map/2_int_8**12
   END FUNCTION generate_larnv_seed

END MODULE block_matrix_api

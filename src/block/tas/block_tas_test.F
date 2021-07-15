!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE block_tas_test

   !! testing infrastructure for tall-and-skinny matrices
   USE block_matrix_api, ONLY: &
      block_matrix_release, block_matrix_nblkcols_total, block_matrix_nblkrows_total, &
      block_matrix_row_block_sizes, block_matrix_col_block_sizes, &
      block_matrix_distribution_release, &
      block_matrix_multiply, &
      block_matrix_type, block_matrix_distribution_obj, &
      block_matrix_no_transpose, block_matrix_transpose, &
      block_matrix_distribution_new, &
      block_matrix_create, &
      block_matrix_finalize, &
      block_matrix_checksum, &
      block_matrix_maxabs, &
      block_matrix_add, &
      block_matrix_complete_redistribute, &
      block_matrix_lapack_larnv, block_matrix_set_larnv_seed, block_matrix_get_name, &
      block_matrix_row_blk_size_obj, block_matrix_col_blk_size_obj, &
      block_matrix_set_larnv_seed
   USE block_tas_base, ONLY: &
      block_tas_convert_to_dbcsr, block_tas_create, block_tas_distribution_new, &
      block_tas_finalize, block_tas_get_stored_coordinates, block_tas_nblkcols_total, &
      block_tas_nblkrows_total, block_tas_put_block, block_tas_info
   USE block_tas_types, ONLY: block_tas_distribution_type, &
                              block_tas_type
   USE block_tas_global, ONLY: block_tas_blk_size_arb, &
                               block_tas_dist_cyclic, &
                               block_tas_default_distvec
   USE block_tas_mm, ONLY: block_tas_multiply
   USE block_tas_split, ONLY: block_tas_mp_comm, &
                              block_tas_get_split_info
   USE block_tas_util, ONLY: invert_transpose_flag
   USE kinds, ONLY: int_8, real_8, dp
   USE message_passing, ONLY: mp_environ, &
                            mp_cart_create, &
                            mp_comm_free

#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: &
      block_tas_benchmark_mm, &
      block_tas_checksum, &
      block_tas_random_bsizes, &
      block_tas_setup_test_matrix, &
      block_tas_test_mm, &
      block_tas_reset_randmat_seed

   INTEGER, SAVE :: randmat_counter = 0
   INTEGER, PARAMETER, PRIVATE :: rand_seed_init = 12341313

CONTAINS
   SUBROUTINE block_tas_setup_test_matrix(matrix, mp_comm_out, mp_comm, nrows, ncols, rbsizes, cbsizes, &
      !! Setup tall-and-skinny matrix for testing
                                          dist_splitsize, name, sparsity, reuse_comm)

      TYPE(block_tas_type), INTENT(OUT)                    :: matrix
      INTEGER, INTENT(OUT)                               :: mp_comm_out
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER(KIND=int_8), INTENT(IN)                    :: nrows, ncols
      INTEGER, DIMENSION(nrows), INTENT(IN)              :: rbsizes
      INTEGER, DIMENSION(ncols), INTENT(IN)              :: cbsizes
      INTEGER, DIMENSION(2), INTENT(IN)                  :: dist_splitsize
      CHARACTER(len=*), INTENT(IN)                       :: name
      REAL(KIND=real_8), INTENT(IN)                      :: sparsity
      LOGICAL, INTENT(IN), OPTIONAL                      :: reuse_comm

      INTEGER                                            :: col_size, max_col_size, max_nze, &
                                                            max_row_size, mynode, node_holds_blk, &
                                                            numnodes, nze, row_size
      INTEGER(KIND=int_8)                                :: col, col_s, row, row_s, nrow, ncol
      INTEGER, DIMENSION(2)                              :: pcoord, pdims
      LOGICAL                                            :: reuse_comm_prv, tr
      REAL(KIND=real_8), DIMENSION(1)                    :: rn
      REAL(KIND=real_8), ALLOCATABLE, DIMENSION(:, :)    :: values
      TYPE(block_tas_blk_size_arb)                         :: cbsize_obj, rbsize_obj
      TYPE(block_tas_dist_cyclic)                          :: col_dist_obj, row_dist_obj
      TYPE(block_tas_distribution_type)                    :: dist
      CHARACTER(LEN=*), PARAMETER :: routineN = 'block_tas_setup_test_matrix'
      INTEGER :: handle
      INTEGER, DIMENSION(4)                              :: iseed, jseed

      ! we don't reserve blocks prior to putting them, so this time is meaningless and should not
      ! be considered in benchmark!
      CALL timeset(routineN, handle)

      ! Check that the counter was initialised (or has not overflowed)
      CPASSERT(randmat_counter .NE. 0)
      ! the counter goes into the seed. Every new call gives a new random matrix
      randmat_counter = randmat_counter + 1

      IF (PRESENT(reuse_comm)) THEN
         reuse_comm_prv = reuse_comm
      ELSE
         reuse_comm_prv = .FALSE.
      ENDIF

      IF (reuse_comm_prv) THEN
         mp_comm_out = mp_comm
      ELSE
         mp_comm_out = block_tas_mp_comm(mp_comm, nrows, ncols)
      ENDIF

      CALL mp_environ(numnodes, mynode, mp_comm_out)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm_out)

      row_dist_obj = block_tas_dist_cyclic(dist_splitsize(1), pdims(1), nrows)
      col_dist_obj = block_tas_dist_cyclic(dist_splitsize(2), pdims(2), ncols)

      rbsize_obj = block_tas_blk_size_arb(rbsizes)
      cbsize_obj = block_tas_blk_size_arb(cbsizes)

      CALL block_tas_distribution_new(dist, mp_comm_out, row_dist_obj, col_dist_obj)
      CALL block_tas_create(matrix, name, dist=dist, &
                            row_blk_size=rbsize_obj, col_blk_size=cbsize_obj, own_dist=.TRUE.)

      max_row_size = MAXVAL(rbsizes)
      max_col_size = MAXVAL(cbsizes)
      max_nze = max_row_size*max_col_size

      nrow = block_tas_nblkrows_total(matrix)
      ncol = block_tas_nblkcols_total(matrix)

      ALLOCATE (values(max_row_size, max_col_size))

      CALL block_matrix_set_larnv_seed(7, 42, 3, 42, randmat_counter, jseed)

      DO row = 1, block_tas_nblkrows_total(matrix)
         DO col = 1, block_tas_nblkcols_total(matrix)
            CALL dlarnv(1, jseed, 1, rn)
            IF (rn(1) .LT. sparsity) THEN
               tr = .FALSE.
               row_s = row; col_s = col
               CALL block_tas_get_stored_coordinates(matrix, row_s, col_s, node_holds_blk)

               IF (node_holds_blk .EQ. mynode) THEN
                  row_size = rbsize_obj%data(row_s)
                  col_size = cbsize_obj%data(col_s)
                  nze = row_size*col_size
                  CALL block_matrix_set_larnv_seed(INT(row_s), INT(nrow), INT(col_s), INT(ncol), randmat_counter, iseed)
                  CALL dlarnv(1, iseed, max_nze, values)
                  CALL block_tas_put_block(matrix, row_s, col_s, values(1:row_size, 1:col_size))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      CALL block_tas_finalize(matrix)

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE block_tas_benchmark_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, compare_dbcsr, filter_eps, io_unit)
      !! Benchmark routine. Due to random sparsity (as opposed to structured sparsity pattern), this
      !! may not be representative for actual applications.

      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c

      REAL(KIND=real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: io_unit
      LOGICAL, INTENT(IN) :: compare_dbcsr

      INTEGER                                            :: handle1, handle2
      TYPE(block_matrix_type)                                   :: block_matrix_a, block_matrix_b, block_matrix_c, &
                                                            block_matrix_a_mm, block_matrix_b_mm, block_matrix_c_mm
      INTEGER                                            :: mp_comm, comm_dbcsr
      TYPE(block_matrix_distribution_obj)                       :: dist_a, dist_b, dist_c
      INTEGER, DIMENSION(2)                              :: npdims, myploc
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: cd_a, cd_b, cd_c, &
                                                            rd_a, rd_b, rd_c

      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: row_blk_size, col_blk_size

      IF (PRESENT(io_unit)) THEN
      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A)") "starting tall-and-skinny benchmark"
      ENDIF
      ENDIF
      CALL timeset("benchmark_tas_mm", handle1)
      CALL block_tas_multiply(transa, transb, transc, 1.0_dp, matrix_a, matrix_b, &
                              0.0_dp, matrix_c, &
                              filter_eps=filter_eps, unit_nr=io_unit)
      CALL timestop(handle1)
      IF (PRESENT(io_unit)) THEN
      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A)") "tall-and-skinny benchmark completed"
      ENDIF
      ENDIF

      IF (compare_dbcsr) THEN
         CALL block_tas_convert_to_dbcsr(matrix_a, block_matrix_a)
         CALL block_tas_convert_to_dbcsr(matrix_b, block_matrix_b)
         CALL block_tas_convert_to_dbcsr(matrix_c, block_matrix_c)

         CALL block_tas_get_split_info(block_tas_info(matrix_a), mp_comm=mp_comm)
         npdims(:) = 0
         CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbcsr)

         ALLOCATE (rd_a(block_matrix_nblkrows_total(block_matrix_a))); ALLOCATE (cd_a(block_matrix_nblkcols_total(block_matrix_a)))
         ALLOCATE (rd_b(block_matrix_nblkrows_total(block_matrix_b))); ALLOCATE (cd_b(block_matrix_nblkcols_total(block_matrix_b)))
         ALLOCATE (rd_c(block_matrix_nblkrows_total(block_matrix_c))); ALLOCATE (cd_c(block_matrix_nblkcols_total(block_matrix_c)))
         CALL block_tas_default_distvec(INT(block_matrix_nblkrows_total(block_matrix_a)), &
                                        npdims(1), block_matrix_row_block_sizes(block_matrix_a), rd_a)
         CALL block_tas_default_distvec(INT(block_matrix_nblkcols_total(block_matrix_a)), &
                                        npdims(2), block_matrix_col_block_sizes(block_matrix_a), cd_a)
         CALL block_tas_default_distvec(INT(block_matrix_nblkrows_total(block_matrix_b)), &
                                        npdims(1), block_matrix_row_block_sizes(block_matrix_b), rd_b)
         CALL block_tas_default_distvec(INT(block_matrix_nblkcols_total(block_matrix_b)), &
                                        npdims(2), block_matrix_col_block_sizes(block_matrix_b), cd_b)
         CALL block_tas_default_distvec(INT(block_matrix_nblkrows_total(block_matrix_c)), &
                                        npdims(1), block_matrix_row_block_sizes(block_matrix_c), rd_c)
         CALL block_tas_default_distvec(INT(block_matrix_nblkcols_total(block_matrix_c)), &
                                        npdims(2), block_matrix_col_block_sizes(block_matrix_c), cd_c)

         CALL block_matrix_distribution_new(dist_a, comm_dbcsr, rd_a, cd_a, reuse_arrays=.TRUE.)
         CALL block_matrix_distribution_new(dist_b, comm_dbcsr, rd_b, cd_b, reuse_arrays=.TRUE.)
         CALL block_matrix_distribution_new(dist_c, comm_dbcsr, rd_c, cd_c, reuse_arrays=.TRUE.)

         row_blk_size => block_matrix_row_block_sizes(block_matrix_a)
         col_blk_size => block_matrix_col_block_sizes(block_matrix_a)
         CALL block_matrix_create(matrix=block_matrix_a_mm, &
                           name=block_matrix_get_name(block_matrix_a), &
                           dist=dist_a, &
                           row_blk_size=row_blk_size, col_blk_size=col_blk_size)
         row_blk_size => block_matrix_row_block_sizes(block_matrix_b)
         col_blk_size => block_matrix_col_block_sizes(block_matrix_b)
         CALL block_matrix_create(matrix=block_matrix_b_mm, &
                           name=block_matrix_get_name(block_matrix_b), &
                           dist=dist_b, &
                           row_blk_size=row_blk_size, col_blk_size=col_blk_size)
         row_blk_size => block_matrix_row_block_sizes(block_matrix_c)
         col_blk_size => block_matrix_col_block_sizes(block_matrix_c)
         CALL block_matrix_create(matrix=block_matrix_c_mm, &
                           name=block_matrix_get_name(block_matrix_c), &
                           dist=dist_c, &
                           row_blk_size=row_blk_size, col_blk_size=col_blk_size)

         CALL block_matrix_finalize(block_matrix_a_mm)
         CALL block_matrix_finalize(block_matrix_b_mm)
         CALL block_matrix_finalize(block_matrix_c_mm)

         CALL block_matrix_complete_redistribute(block_matrix_a, block_matrix_a_mm)
         CALL block_matrix_complete_redistribute(block_matrix_b, block_matrix_b_mm)
         IF (PRESENT(io_unit)) THEN
         IF (io_unit > 0) THEN
            WRITE (io_unit, "(A)") "starting dbcsr benchmark"
         ENDIF
         ENDIF
         CALL timeset("benchmark_block_mm", handle2)
         CALL block_matrix_multiply(transa, transb, 1.0_dp, block_matrix_a_mm, block_matrix_b_mm, &
                                    0.0_dp, block_matrix_c_mm, filter_eps=filter_eps)
         CALL timestop(handle2)
         IF (PRESENT(io_unit)) THEN
         IF (io_unit > 0) THEN
            WRITE (io_unit, "(A)") "dbcsr benchmark completed"
         ENDIF
         ENDIF

         CALL block_matrix_release(block_matrix_a)
         CALL block_matrix_release(block_matrix_b)
         CALL block_matrix_release(block_matrix_c)
         CALL block_matrix_release(block_matrix_a_mm)
         CALL block_matrix_release(block_matrix_b_mm)
         CALL block_matrix_release(block_matrix_c_mm)
         CALL block_matrix_distribution_release(dist_a)
         CALL block_matrix_distribution_release(dist_b)
         CALL block_matrix_distribution_release(dist_c)

         CALL mp_comm_free(comm_dbcsr)
      ENDIF

   END SUBROUTINE

   SUBROUTINE block_tas_test_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, filter_eps, unit_nr, log_verbose)
      !! Test tall-and-skinny matrix multiplication for accuracy
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(block_tas_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN), OPTIONAL                      :: log_verbose
      REAL(KIND=real_8), INTENT(IN), OPTIONAL            :: filter_eps

      CHARACTER(LEN=1)                                   :: transa_prv, transb_prv
      INTEGER                                            :: comm_dbcsr, io_unit, mp_comm, mynode, &
                                                            numnodes
      INTEGER, DIMENSION(2)                              :: myploc, npdims
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: cd_a, cd_b, cd_c, &
                                                            rd_a, rd_b, rd_c
      REAL(KIND=real_8)                                  :: norm, rc_cs, sq_cs
      TYPE(block_matrix_distribution_obj)                       :: dist_a, dist_b, dist_c
      TYPE(block_matrix_type)                                   :: block_matrix_a, block_matrix_a_mm, block_matrix_b, &
                                                            block_matrix_b_mm, block_matrix_c, block_matrix_c_mm, &
                                                            block_matrix_c_mm_check
      REAL(KIND=real_8), PARAMETER :: test_tol = 1.0E-10_real_8

      CALL block_tas_get_split_info(block_tas_info(matrix_a), mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      io_unit = -1
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL block_tas_multiply(transa, transb, transc, 1.0_dp, matrix_a, matrix_b, &
                              0.0_dp, matrix_c, &
                              filter_eps=filter_eps, unit_nr=io_unit, log_verbose=log_verbose, optimize_dist=.TRUE.)

      CALL block_tas_convert_to_dbcsr(matrix_a, block_matrix_a)
      CALL block_tas_convert_to_dbcsr(matrix_b, block_matrix_b)
      CALL block_tas_convert_to_dbcsr(matrix_c, block_matrix_c)

      npdims(:) = 0
      CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbcsr)

      ALLOCATE (rd_a(block_matrix_nblkrows_total(block_matrix_a))); ALLOCATE (cd_a(block_matrix_nblkcols_total(block_matrix_a)))
      ALLOCATE (rd_b(block_matrix_nblkrows_total(block_matrix_b))); ALLOCATE (cd_b(block_matrix_nblkcols_total(block_matrix_b)))
      ALLOCATE (rd_c(block_matrix_nblkrows_total(block_matrix_c))); ALLOCATE (cd_c(block_matrix_nblkcols_total(block_matrix_c)))
      CALL block_tas_default_distvec(INT(block_matrix_nblkrows_total(block_matrix_a)), &
                                     npdims(1), block_matrix_row_block_sizes(block_matrix_a), rd_a)
      CALL block_tas_default_distvec(INT(block_matrix_nblkcols_total(block_matrix_a)), &
                                     npdims(2), block_matrix_col_block_sizes(block_matrix_a), cd_a)
      CALL block_tas_default_distvec(INT(block_matrix_nblkrows_total(block_matrix_b)), &
                                     npdims(1), block_matrix_row_block_sizes(block_matrix_b), rd_b)
      CALL block_tas_default_distvec(INT(block_matrix_nblkcols_total(block_matrix_b)), &
                                     npdims(2), block_matrix_col_block_sizes(block_matrix_b), cd_b)
      CALL block_tas_default_distvec(INT(block_matrix_nblkrows_total(block_matrix_c)), &
                                     npdims(1), block_matrix_row_block_sizes(block_matrix_c), rd_c)
      CALL block_tas_default_distvec(INT(block_matrix_nblkcols_total(block_matrix_c)), &
                                     npdims(2), block_matrix_col_block_sizes(block_matrix_c), cd_c)

      CALL block_matrix_distribution_new(dist_a, comm_dbcsr, rd_a, cd_a, reuse_arrays=.TRUE.)
      CALL block_matrix_distribution_new(dist_b, comm_dbcsr, rd_b, cd_b, reuse_arrays=.TRUE.)
      CALL block_matrix_distribution_new(dist_c, comm_dbcsr, rd_c, cd_c, reuse_arrays=.TRUE.)

      CALL block_matrix_create(matrix=block_matrix_a_mm, name="matrix a", dist=dist_a, &
                        row_blk_size_obj=block_matrix_row_blk_size_obj(block_matrix_a), &
                        col_blk_size_obj=block_matrix_col_blk_size_obj(block_matrix_a))

      CALL block_matrix_create(matrix=block_matrix_b_mm, name="matrix b", dist=dist_b, &
                        row_blk_size_obj=block_matrix_row_blk_size_obj(block_matrix_b), &
                        col_blk_size_obj=block_matrix_col_blk_size_obj(block_matrix_b))

      CALL block_matrix_create(matrix=block_matrix_c_mm, name="matrix c", dist=dist_c, &
                        row_blk_size_obj=block_matrix_row_blk_size_obj(block_matrix_c), &
                        col_blk_size_obj=block_matrix_col_blk_size_obj(block_matrix_c))

      CALL block_matrix_create(matrix=block_matrix_c_mm_check, name="matrix c check", &
                        dist=dist_c, &
                        row_blk_size_obj=block_matrix_row_blk_size_obj(block_matrix_c), &
                        col_blk_size_obj=block_matrix_col_blk_size_obj(block_matrix_c))

      CALL block_matrix_finalize(block_matrix_a_mm)
      CALL block_matrix_finalize(block_matrix_b_mm)
      CALL block_matrix_finalize(block_matrix_c_mm)
      CALL block_matrix_finalize(block_matrix_c_mm_check)

      CALL block_matrix_complete_redistribute(block_matrix_a, block_matrix_a_mm)
      CALL block_matrix_complete_redistribute(block_matrix_b, block_matrix_b_mm)
      CALL block_matrix_complete_redistribute(block_matrix_c, block_matrix_c_mm_check)

      transa_prv = transa; transb_prv = transb

      IF (transc == block_matrix_no_transpose) THEN
         CALL block_matrix_multiply(transa_prv, transb_prv, 1.0_dp, &
                             block_matrix_a_mm, block_matrix_b_mm, &
                             0.0_dp, block_matrix_c_mm, filter_eps=filter_eps)
      ELSEIF (transc == block_matrix_transpose) THEN
         CALL invert_transpose_flag(transa_prv)
         CALL invert_transpose_flag(transb_prv)
         CALL block_matrix_multiply(transb_prv, transa_prv, 1.0_dp, &
                             block_matrix_b_mm, block_matrix_a_mm, &
                             0.0_dp, block_matrix_c_mm, filter_eps=filter_eps)
      ENDIF

      sq_cs = block_matrix_checksum(block_matrix_c_mm)
      rc_cs = block_matrix_checksum(block_matrix_c_mm_check)
      CALL block_matrix_add(block_matrix_c_mm_check, block_matrix_c_mm, -1.0_real_8, 1.0_real_8)
      norm = block_matrix_maxabs(block_matrix_c_mm_check)

      IF (io_unit > 0) THEN
      IF (ABS(norm) .GT. test_tol) THEN
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(block_matrix_get_name(matrix_a%matrix)), &
                                                   transa, ' X ', &
                                                   TRIM(block_matrix_get_name(matrix_b%matrix)), &
                                                   transb, 'failed!'
         WRITE (io_unit, "(A,1X,E9.2,1X,E9.2)") "checksums", sq_cs, rc_cs
         WRITE (io_unit, "(A,1X,E9.2)") "difference norm", norm
         CPABORT("")
      ELSE
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(block_matrix_get_name(matrix_a%matrix)), &
                                                   transa, ' X ', &
                                                   TRIM(block_matrix_get_name(matrix_b%matrix)), &
                                                   transb, 'passed!'
         WRITE (io_unit, "(A,1X,E9.2,1X,E9.2)") "checksums", sq_cs, rc_cs
         WRITE (io_unit, "(A,1X,E9.2)") "difference norm", norm
      ENDIF
      ENDIF

      CALL block_matrix_release(block_matrix_a)
      CALL block_matrix_release(block_matrix_a_mm)
      CALL block_matrix_release(block_matrix_b)
      CALL block_matrix_release(block_matrix_b_mm)
      CALL block_matrix_release(block_matrix_c)
      CALL block_matrix_release(block_matrix_c_mm)
      CALL block_matrix_release(block_matrix_c_mm_check)

      CALL block_matrix_distribution_release(dist_a)
      CALL block_matrix_distribution_release(dist_b)
      CALL block_matrix_distribution_release(dist_c)

      CALL mp_comm_free(comm_dbcsr)

   END SUBROUTINE

   FUNCTION block_tas_checksum(matrix, local, pos)
      !! Calculate checksum of tall-and-skinny matrix consistent with block_matrix_checksum
      TYPE(block_tas_type), INTENT(IN) :: matrix
      LOGICAL, INTENT(IN), OPTIONAL  :: local, pos
      TYPE(block_matrix_type)               :: block_matrix_m
      REAL(KIND=real_8)              :: block_tas_checksum

      CALL block_tas_convert_to_dbcsr(matrix, block_matrix_m)
      block_tas_checksum = block_matrix_checksum(block_matrix_m, local, pos)
      CALL block_matrix_release(block_matrix_m)
   END FUNCTION

   SUBROUTINE block_tas_random_bsizes(sizes, repeat, block_sizes)
      !! Create random block sizes
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes
      INTEGER, INTENT(IN)                                :: repeat
      INTEGER, DIMENSION(:), INTENT(OUT)                 :: block_sizes

      INTEGER                                            :: d, size_i

      DO d = 1, SIZE(block_sizes)
         size_i = MOD((d - 1)/repeat, SIZE(sizes)) + 1
         block_sizes(d) = sizes(size_i)
      ENDDO
   END SUBROUTINE

   SUBROUTINE block_tas_reset_randmat_seed()
      !! Reset the seed used for generating random matrices to default value
      randmat_counter = rand_seed_init
   END SUBROUTINE

END MODULE

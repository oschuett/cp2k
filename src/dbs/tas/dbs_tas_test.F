!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE dbs_tas_test

   !! testing infrastructure for tall-and-skinny matrices
   USE dbs_types, ONLY: dbs_type_real_8
   USE dbs_data_methods, ONLY: dbs_scalar
   USE dbs_methods, ONLY: &
      dbs_release, dbs_nblkcols_total, dbs_nblkrows_total, dbs_row_block_sizes, dbs_col_block_sizes, &
      dbs_mp_release, dbs_distribution_release
   USE dbs_multiply_api, ONLY: dbs_multiply
   USE dbs_tas_base, ONLY: &
      dbs_tas_convert_to_dbcsr, dbs_tas_create, dbs_tas_distribution_new, &
      dbs_tas_finalize, dbs_tas_get_stored_coordinates, dbs_tas_nblkcols_total, &
      dbs_tas_nblkrows_total, dbs_tas_put_block, dbs_tas_info
   USE dbs_tas_types, ONLY: dbs_tas_distribution_type, &
                              dbs_tas_type
   USE dbs_tas_global, ONLY: dbs_tas_blk_size_arb, &
                               dbs_tas_dist_cyclic, &
                               dbs_tas_default_distvec
   USE dbs_tas_mm, ONLY: dbs_tas_multiply
   USE dbs_tas_split, ONLY: dbs_tas_mp_comm, &
                              dbs_tas_get_split_info
   USE dbs_tas_util, ONLY: dbs_mp_environ, &
                             invert_transpose_flag
   USE dbs_types, ONLY: &
      dbs_type, dbs_distribution_obj, dbs_mp_obj, dbs_no_transpose, dbs_transpose, &
      dbs_type_no_symmetry
   USE dbs_kinds, ONLY: int_8, &
                          real_8
   USE dbs_mpiwrap, ONLY: mp_environ, &
                            mp_cart_create, &
                            mp_comm_free
   USE dbs_dist_methods, ONLY: dbs_distribution_new
   USE dbs_work_operations, ONLY: dbs_create, &
                                    dbs_finalize
   USE dbs_dist_util, ONLY: dbs_checksum
   USE dbs_operations, ONLY: dbs_maxabs, &
                               dbs_add
   USE dbs_transformations, ONLY: dbs_complete_redistribute
   USE dbs_blas_operations, ONLY: &
      dbs_lapack_larnv, set_larnv_seed
#include "base/dbs_base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: &
      dbs_tas_benchmark_mm, &
      dbs_tas_checksum, &
      dbs_tas_random_bsizes, &
      dbs_tas_setup_test_matrix, &
      dbs_tas_test_mm, &
      dbs_tas_reset_randmat_seed

   INTEGER, SAVE :: randmat_counter = 0
   INTEGER, PARAMETER, PRIVATE :: rand_seed_init = 12341313

CONTAINS
   SUBROUTINE dbs_tas_setup_test_matrix(matrix, mp_comm_out, mp_comm, nrows, ncols, rbsizes, cbsizes, &
      !! Setup tall-and-skinny matrix for testing
                                          dist_splitsize, name, sparsity, reuse_comm)

      TYPE(dbs_tas_type), INTENT(OUT)                    :: matrix
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
      TYPE(dbs_tas_blk_size_arb)                         :: cbsize_obj, rbsize_obj
      TYPE(dbs_tas_dist_cyclic)                          :: col_dist_obj, row_dist_obj
      TYPE(dbs_tas_distribution_type)                    :: dist
      CHARACTER(LEN=*), PARAMETER :: routineN = 'dbs_tas_setup_test_matrix'
      INTEGER :: handle
      INTEGER, DIMENSION(4)                              :: iseed, jseed

      ! we don't reserve blocks prior to putting them, so this time is meaningless and should not
      ! be considered in benchmark!
      CALL timeset(routineN, handle)

      ! Check that the counter was initialised (or has not overflowed)
      DBCSR_ASSERT(randmat_counter .NE. 0)
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
         mp_comm_out = dbs_tas_mp_comm(mp_comm, nrows, ncols)
      ENDIF

      CALL mp_environ(numnodes, mynode, mp_comm_out)
      CALL mp_environ(numnodes, pdims, pcoord, mp_comm_out)

      row_dist_obj = dbs_tas_dist_cyclic(dist_splitsize(1), pdims(1), nrows)
      col_dist_obj = dbs_tas_dist_cyclic(dist_splitsize(2), pdims(2), ncols)

      rbsize_obj = dbs_tas_blk_size_arb(rbsizes)
      cbsize_obj = dbs_tas_blk_size_arb(cbsizes)

      CALL dbs_tas_distribution_new(dist, mp_comm_out, row_dist_obj, col_dist_obj)
      CALL dbs_tas_create(matrix, name, dist=dist, data_type=dbs_type_real_8, &
                            row_blk_size=rbsize_obj, col_blk_size=cbsize_obj, own_dist=.TRUE.)

      max_row_size = MAXVAL(rbsizes)
      max_col_size = MAXVAL(cbsizes)
      max_nze = max_row_size*max_col_size

      nrow = dbs_tas_nblkrows_total(matrix)
      ncol = dbs_tas_nblkcols_total(matrix)

      ALLOCATE (values(max_row_size, max_col_size))

      CALL set_larnv_seed(7, 42, 3, 42, randmat_counter, jseed)

      DO row = 1, dbs_tas_nblkrows_total(matrix)
         DO col = 1, dbs_tas_nblkcols_total(matrix)
            CALL dlarnv(1, jseed, 1, rn)
            IF (rn(1) .LT. sparsity) THEN
               tr = .FALSE.
               row_s = row; col_s = col
               CALL dbs_tas_get_stored_coordinates(matrix, row_s, col_s, node_holds_blk)

               IF (node_holds_blk .EQ. mynode) THEN
                  row_size = rbsize_obj%data(row_s)
                  col_size = cbsize_obj%data(col_s)
                  nze = row_size*col_size
                  CALL set_larnv_seed(INT(row_s), INT(nrow), INT(col_s), INT(ncol), randmat_counter, iseed)
                  CALL dlarnv(1, iseed, max_nze, values)
                  CALL dbs_tas_put_block(matrix, row_s, col_s, values(1:row_size, 1:col_size))
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      CALL dbs_tas_finalize(matrix)

      CALL timestop(handle)

   END SUBROUTINE

   SUBROUTINE dbs_tas_benchmark_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, compare_dbcsr, filter_eps, io_unit)
      !! Benchmark routine. Due to random sparsity (as opposed to structured sparsity pattern), this
      !! may not be representative for actual applications.

      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(dbs_tas_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c

      REAL(KIND=real_8), INTENT(IN), OPTIONAL            :: filter_eps
      INTEGER, INTENT(IN), OPTIONAL                      :: io_unit
      LOGICAL, INTENT(IN) :: compare_dbcsr

      INTEGER                                            :: handle1, handle2
      TYPE(dbs_type)                                   :: dbs_a, dbs_b, dbs_c, &
                                                            dbs_a_mm, dbs_b_mm, dbs_c_mm
      INTEGER                                            :: mp_comm, comm_dbcsr
      TYPE(dbs_distribution_obj)                       :: dist_a, dist_b, dist_c
      INTEGER, DIMENSION(2)                              :: npdims, myploc
      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: cd_a, cd_b, cd_c, &
                                                            rd_a, rd_b, rd_c
      TYPE(dbs_mp_obj)                                 :: mp_environ_tmp

      INTEGER, DIMENSION(:), POINTER, CONTIGUOUS         :: row_blk_size, col_blk_size

      IF (PRESENT(io_unit)) THEN
      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A)") "starting tall-and-skinny benchmark"
      ENDIF
      ENDIF
      CALL timeset("benchmark_tas_mm", handle1)
      CALL dbs_tas_multiply(transa, transb, transc, dbs_scalar(1.0_real_8), matrix_a, matrix_b, &
                              dbs_scalar(0.0_real_8), matrix_c, &
                              filter_eps=filter_eps, unit_nr=io_unit)
      CALL timestop(handle1)
      IF (PRESENT(io_unit)) THEN
      IF (io_unit > 0) THEN
         WRITE (io_unit, "(A)") "tall-and-skinny benchmark completed"
      ENDIF
      ENDIF

      IF (compare_dbcsr) THEN
         CALL dbs_tas_convert_to_dbcsr(matrix_a, dbs_a)
         CALL dbs_tas_convert_to_dbcsr(matrix_b, dbs_b)
         CALL dbs_tas_convert_to_dbcsr(matrix_c, dbs_c)

         CALL dbs_tas_get_split_info(dbs_tas_info(matrix_a), mp_comm=mp_comm)
         npdims(:) = 0
         CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbcsr)

         ALLOCATE (rd_a(dbs_nblkrows_total(dbs_a))); ALLOCATE (cd_a(dbs_nblkcols_total(dbs_a)))
         ALLOCATE (rd_b(dbs_nblkrows_total(dbs_b))); ALLOCATE (cd_b(dbs_nblkcols_total(dbs_b)))
         ALLOCATE (rd_c(dbs_nblkrows_total(dbs_c))); ALLOCATE (cd_c(dbs_nblkcols_total(dbs_c)))
         CALL dbs_tas_default_distvec(INT(dbs_nblkrows_total(dbs_a)), npdims(1), dbs_row_block_sizes(dbs_a), rd_a)
         CALL dbs_tas_default_distvec(INT(dbs_nblkcols_total(dbs_a)), npdims(2), dbs_col_block_sizes(dbs_a), cd_a)
         CALL dbs_tas_default_distvec(INT(dbs_nblkrows_total(dbs_b)), npdims(1), dbs_row_block_sizes(dbs_b), rd_b)
         CALL dbs_tas_default_distvec(INT(dbs_nblkcols_total(dbs_b)), npdims(2), dbs_col_block_sizes(dbs_b), cd_b)
         CALL dbs_tas_default_distvec(INT(dbs_nblkrows_total(dbs_c)), npdims(1), dbs_row_block_sizes(dbs_c), rd_c)
         CALL dbs_tas_default_distvec(INT(dbs_nblkcols_total(dbs_c)), npdims(2), dbs_col_block_sizes(dbs_c), cd_c)

         mp_environ_tmp = dbs_mp_environ(comm_dbcsr)
         CALL dbs_distribution_new(dist_a, mp_environ_tmp, rd_a, cd_a, reuse_arrays=.TRUE.)
         CALL dbs_distribution_new(dist_b, mp_environ_tmp, rd_b, cd_b, reuse_arrays=.TRUE.)
         CALL dbs_distribution_new(dist_c, mp_environ_tmp, rd_c, cd_c, reuse_arrays=.TRUE.)
         CALL dbs_mp_release(mp_environ_tmp)

         row_blk_size => dbs_row_block_sizes(dbs_a)
         col_blk_size => dbs_col_block_sizes(dbs_a)
         CALL dbs_create(matrix=dbs_a_mm, name=dbs_a%name, dist=dist_a, matrix_type=dbs_type_no_symmetry, &
                           row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                           data_type=dbs_type_real_8)
         row_blk_size => dbs_row_block_sizes(dbs_b)
         col_blk_size => dbs_col_block_sizes(dbs_b)
         CALL dbs_create(matrix=dbs_b_mm, name=dbs_b%name, dist=dist_b, matrix_type=dbs_type_no_symmetry, &
                           row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                           data_type=dbs_type_real_8)
         row_blk_size => dbs_row_block_sizes(dbs_c)
         col_blk_size => dbs_col_block_sizes(dbs_c)
         CALL dbs_create(matrix=dbs_c_mm, name=dbs_c%name, dist=dist_c, matrix_type=dbs_type_no_symmetry, &
                           row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                           data_type=dbs_type_real_8)

         CALL dbs_finalize(dbs_a_mm)
         CALL dbs_finalize(dbs_b_mm)
         CALL dbs_finalize(dbs_c_mm)

         CALL dbs_complete_redistribute(dbs_a, dbs_a_mm)
         CALL dbs_complete_redistribute(dbs_b, dbs_b_mm)
         IF (PRESENT(io_unit)) THEN
         IF (io_unit > 0) THEN
            WRITE (io_unit, "(A)") "starting dbcsr benchmark"
         ENDIF
         ENDIF
         CALL timeset("benchmark_dbs_mm", handle2)
         CALL dbs_multiply(transa, transb, dbs_scalar(1.0_real_8), dbs_a_mm, dbs_b_mm, &
                             dbs_scalar(0.0_real_8), dbs_c_mm, filter_eps=filter_eps)
         CALL timestop(handle2)
         IF (PRESENT(io_unit)) THEN
         IF (io_unit > 0) THEN
            WRITE (io_unit, "(A)") "dbcsr benchmark completed"
         ENDIF
         ENDIF

         CALL dbs_release(dbs_a)
         CALL dbs_release(dbs_b)
         CALL dbs_release(dbs_c)
         CALL dbs_release(dbs_a_mm)
         CALL dbs_release(dbs_b_mm)
         CALL dbs_release(dbs_c_mm)
         CALL dbs_distribution_release(dist_a)
         CALL dbs_distribution_release(dist_b)
         CALL dbs_distribution_release(dist_c)

         CALL mp_comm_free(comm_dbcsr)
      ENDIF

   END SUBROUTINE

   SUBROUTINE dbs_tas_test_mm(transa, transb, transc, matrix_a, matrix_b, matrix_c, filter_eps, unit_nr, log_verbose)
      !! Test tall-and-skinny matrix multiplication for accuracy
      CHARACTER(LEN=1), INTENT(IN)                       :: transa, transb, transc
      TYPE(dbs_tas_type), INTENT(INOUT)                  :: matrix_a, matrix_b, matrix_c
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
      TYPE(dbs_distribution_obj)                       :: dist_a, dist_b, dist_c
      TYPE(dbs_mp_obj)                                 :: mp_environ_tmp
      TYPE(dbs_type)                                   :: dbs_a, dbs_a_mm, dbs_b, &
                                                            dbs_b_mm, dbs_c, dbs_c_mm, &
                                                            dbs_c_mm_check
      REAL(KIND=real_8), PARAMETER :: test_tol = 1.0E-10_real_8

      CALL dbs_tas_get_split_info(dbs_tas_info(matrix_a), mp_comm=mp_comm)
      CALL mp_environ(numnodes, mynode, mp_comm)
      io_unit = -1
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL dbs_tas_multiply(transa, transb, transc, dbs_scalar(1.0_real_8), matrix_a, matrix_b, &
                              dbs_scalar(0.0_real_8), matrix_c, &
                              filter_eps=filter_eps, unit_nr=io_unit, log_verbose=log_verbose, optimize_dist=.TRUE.)

      CALL dbs_tas_convert_to_dbcsr(matrix_a, dbs_a)
      CALL dbs_tas_convert_to_dbcsr(matrix_b, dbs_b)
      CALL dbs_tas_convert_to_dbcsr(matrix_c, dbs_c)

      npdims(:) = 0
      CALL mp_cart_create(mp_comm, 2, npdims, myploc, comm_dbcsr)

      ALLOCATE (rd_a(dbs_nblkrows_total(dbs_a))); ALLOCATE (cd_a(dbs_nblkcols_total(dbs_a)))
      ALLOCATE (rd_b(dbs_nblkrows_total(dbs_b))); ALLOCATE (cd_b(dbs_nblkcols_total(dbs_b)))
      ALLOCATE (rd_c(dbs_nblkrows_total(dbs_c))); ALLOCATE (cd_c(dbs_nblkcols_total(dbs_c)))
      CALL dbs_tas_default_distvec(INT(dbs_nblkrows_total(dbs_a)), npdims(1), dbs_row_block_sizes(dbs_a), rd_a)
      CALL dbs_tas_default_distvec(INT(dbs_nblkcols_total(dbs_a)), npdims(2), dbs_col_block_sizes(dbs_a), cd_a)
      CALL dbs_tas_default_distvec(INT(dbs_nblkrows_total(dbs_b)), npdims(1), dbs_row_block_sizes(dbs_b), rd_b)
      CALL dbs_tas_default_distvec(INT(dbs_nblkcols_total(dbs_b)), npdims(2), dbs_col_block_sizes(dbs_b), cd_b)
      CALL dbs_tas_default_distvec(INT(dbs_nblkrows_total(dbs_c)), npdims(1), dbs_row_block_sizes(dbs_c), rd_c)
      CALL dbs_tas_default_distvec(INT(dbs_nblkcols_total(dbs_c)), npdims(2), dbs_col_block_sizes(dbs_c), cd_c)

      mp_environ_tmp = dbs_mp_environ(comm_dbcsr)
      CALL dbs_distribution_new(dist_a, mp_environ_tmp, rd_a, cd_a, reuse_arrays=.TRUE.)
      CALL dbs_distribution_new(dist_b, mp_environ_tmp, rd_b, cd_b, reuse_arrays=.TRUE.)
      CALL dbs_distribution_new(dist_c, mp_environ_tmp, rd_c, cd_c, reuse_arrays=.TRUE.)
      CALL dbs_mp_release(mp_environ_tmp)

      CALL dbs_create(matrix=dbs_a_mm, name="matrix a", dist=dist_a, matrix_type=dbs_type_no_symmetry, &
                        row_blk_size_obj=dbs_a%row_blk_size, col_blk_size_obj=dbs_a%col_blk_size, &
                        data_type=dbs_type_real_8)

      CALL dbs_create(matrix=dbs_b_mm, name="matrix b", dist=dist_b, matrix_type=dbs_type_no_symmetry, &
                        row_blk_size_obj=dbs_b%row_blk_size, col_blk_size_obj=dbs_b%col_blk_size, &
                        data_type=dbs_type_real_8)

      CALL dbs_create(matrix=dbs_c_mm, name="matrix c", dist=dist_c, matrix_type=dbs_type_no_symmetry, &
                        row_blk_size_obj=dbs_c%row_blk_size, col_blk_size_obj=dbs_c%col_blk_size, &
                        data_type=dbs_type_real_8)

      CALL dbs_create(matrix=dbs_c_mm_check, name="matrix c check", dist=dist_c, matrix_type=dbs_type_no_symmetry, &
                        row_blk_size_obj=dbs_c%row_blk_size, col_blk_size_obj=dbs_c%col_blk_size, &
                        data_type=dbs_type_real_8)

      CALL dbs_finalize(dbs_a_mm)
      CALL dbs_finalize(dbs_b_mm)
      CALL dbs_finalize(dbs_c_mm)
      CALL dbs_finalize(dbs_c_mm_check)

      CALL dbs_complete_redistribute(dbs_a, dbs_a_mm)
      CALL dbs_complete_redistribute(dbs_b, dbs_b_mm)
      CALL dbs_complete_redistribute(dbs_c, dbs_c_mm_check)

      transa_prv = transa; transb_prv = transb

      IF (transc == dbs_no_transpose) THEN
         CALL dbs_multiply(transa_prv, transb_prv, dbs_scalar(1.0_real_8), &
                             dbs_a_mm, dbs_b_mm, dbs_scalar(0.0_real_8), dbs_c_mm, filter_eps=filter_eps)
      ELSEIF (transc == dbs_transpose) THEN
         CALL invert_transpose_flag(transa_prv)
         CALL invert_transpose_flag(transb_prv)
         CALL dbs_multiply(transb_prv, transa_prv, dbs_scalar(1.0_real_8), &
                             dbs_b_mm, dbs_a_mm, dbs_scalar(0.0_real_8), dbs_c_mm, filter_eps=filter_eps)
      ENDIF

      sq_cs = dbs_checksum(dbs_c_mm)
      rc_cs = dbs_checksum(dbs_c_mm_check)
      CALL dbs_add(dbs_c_mm_check, dbs_c_mm, -1.0_real_8, 1.0_real_8)
      norm = dbs_maxabs(dbs_c_mm_check)

      IF (io_unit > 0) THEN
      IF (ABS(norm) .GT. test_tol) THEN
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(matrix_a%matrix%name), transa, ' X ', TRIM(matrix_b%matrix%name), &
            transb, 'failed!'
         WRITE (io_unit, "(A,1X,E9.2,1X,E9.2)") "checksums", sq_cs, rc_cs
         WRITE (io_unit, "(A,1X,E9.2)") "difference norm", norm
         DBCSR_ABORT("")
      ELSE
         WRITE (io_unit, '(A, A, A, A, A, 1X, A)') TRIM(matrix_a%matrix%name), transa, ' X ', TRIM(matrix_b%matrix%name), &
            transb, 'passed!'
         WRITE (io_unit, "(A,1X,E9.2,1X,E9.2)") "checksums", sq_cs, rc_cs
         WRITE (io_unit, "(A,1X,E9.2)") "difference norm", norm
      ENDIF
      ENDIF

      CALL dbs_release(dbs_a)
      CALL dbs_release(dbs_a_mm)
      CALL dbs_release(dbs_b)
      CALL dbs_release(dbs_b_mm)
      CALL dbs_release(dbs_c)
      CALL dbs_release(dbs_c_mm)
      CALL dbs_release(dbs_c_mm_check)

      CALL dbs_distribution_release(dist_a)
      CALL dbs_distribution_release(dist_b)
      CALL dbs_distribution_release(dist_c)

      CALL mp_comm_free(comm_dbcsr)

   END SUBROUTINE

   FUNCTION dbs_tas_checksum(matrix, local, pos)
      !! Calculate checksum of tall-and-skinny matrix consistent with dbs_checksum
      TYPE(dbs_tas_type), INTENT(IN) :: matrix
      LOGICAL, INTENT(IN), OPTIONAL  :: local, pos
      TYPE(dbs_type)               :: dbs_m
      REAL(KIND=real_8)              :: dbs_tas_checksum

      CALL dbs_tas_convert_to_dbcsr(matrix, dbs_m)
      dbs_tas_checksum = dbs_checksum(dbs_m, local, pos)
      CALL dbs_release(dbs_m)
   END FUNCTION

   SUBROUTINE dbs_tas_random_bsizes(sizes, repeat, block_sizes)
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

   SUBROUTINE dbs_tas_reset_randmat_seed()
      !! Reset the seed used for generating random matrices to default value
      randmat_counter = rand_seed_init
   END SUBROUTINE

END MODULE

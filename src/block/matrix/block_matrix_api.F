!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
MODULE block_matrix_api

 USE dbcsr_types, ONLY: &
     block_matrix_type => dbcsr_type

 USE dbcsr_iterator_operations, ONLY: &
    block_matrix_iterator_blocks_left => dbcsr_iterator_blocks_left, &
    block_matrix_iterator_next_block => dbcsr_iterator_next_block, &
    block_matrix_iterator_start => dbcsr_iterator_start, &
    block_matrix_iterator_stop => dbcsr_iterator_stop

 USE dbcsr_transformations, ONLY: &
    dbcsr_desymmetrize_deep, &
    block_matrix_redistribute => dbcsr_redistribute, &
    block_matrix_complete_redistribute => dbcsr_complete_redistribute

 USE dbcsr_operations, ONLY: &
      block_matrix_set    => dbcsr_set, &
      block_matrix_filter => dbcsr_filter, &
      block_matrix_scale => dbcsr_scale, &
      block_matrix_get_info => dbcsr_get_info, &
      block_matrix_clear => dbcsr_clear, &
      block_matrix_copy => dbcsr_copy, &
      block_matrix_add => dbcsr_add, &
      block_matrix_maxabs => dbcsr_maxabs, &
      block_matrix_zero => dbcsr_zero

 USE dbcsr_block_access, ONLY: &
      block_matrix_get_block_p => dbcsr_get_block_p, &
      block_matrix_put_block => dbcsr_put_block, &
      block_matrix_reserve_blocks => dbcsr_reserve_blocks, &
      block_matrix_put_block => dbcsr_put_block

 USE dbcsr_methods, ONLY: &
      block_matrix_release              => dbcsr_release, &
      block_matrix_get_data_type        => dbcsr_get_data_type, &
      block_matrix_mp_release           => dbcsr_mp_release, &
      block_matrix_nblkcols_total       => dbcsr_nblkcols_total, &
      block_matrix_nblkrows_total       => dbcsr_nblkrows_total, &
      block_matrix_nfullrows_total      => dbcsr_nfullrows_total, &
      block_matrix_nfullcols_total      => dbcsr_nfullcols_total, &
      block_matrix_get_data_size        => dbcsr_get_data_size, &
      block_matrix_get_num_blocks       => dbcsr_get_num_blocks, &
      block_matrix_get_nze              => dbcsr_get_nze, &
      block_mp_release                  => dbcsr_mp_release, &
      block_matrix_get_matrix_type      => dbcsr_get_matrix_type, &

      block_matrix_blk_column_size      => dbcsr_blk_column_size, &
      block_matrix_blk_row_size         => dbcsr_blk_row_size, &
      block_matrix_row_block_sizes      => dbcsr_row_block_sizes, &
      block_matrix_col_block_sizes      => dbcsr_col_block_sizes

 USE dbcsr_dist_methods, ONLY: &
      block_matrix_distribution_get => dbcsr_distribution_get, &
      block_matrix_distribution_col_dist=> dbcsr_distribution_col_dist, &
      block_matrix_distribution_new     => dbcsr_distribution_new, &
      block_matrix_distribution_row_dist=> dbcsr_distribution_row_dist, &
      block_matrix_distribution_hold    => dbcsr_distribution_hold, &
      block_matrix_distribution_release => dbcsr_distribution_release

 USE dbcsr_dist_operations, ONLY: &
      block_matrix_get_stored_coordinates => dbcsr_get_stored_coordinates, &
      checker_tr

 USE dbcsr_data_methods, ONLY: &
      block_scalar => dbcsr_scalar, &
      block_scalar_zero => dbcsr_scalar_zero, &
      block_matrix_data_new             => dbcsr_data_new, &
      block_matrix_data_release         => dbcsr_data_release, &
      block_matrix_type_1d_to_2d        => dbcsr_type_1d_to_2d, &
      block_matrix_data_clear_pointer   => dbcsr_data_clear_pointer, &
      block_matrix_data_init            => dbcsr_data_init, &
      block_matrix_data_get_sizes       => dbcsr_data_get_sizes, &
      block_matrix_data_set_pointer     => dbcsr_data_set_pointer

 USE dbcsr_data_methods_low, ONLY: &
      block_matrix_scalar_multiply => dbcsr_scalar_multiply, &
      internal_data_allocate, &
      internal_data_deallocate, &
      block_matrix_get_data_p_2d_d      => dbcsr_get_data_p_2d_d, &
      block_matrix_get_data_p_2d_z      => dbcsr_get_data_p_2d_z, &
      block_matrix_get_data_p_2d_c      => dbcsr_get_data_p_2d_c, &
      block_matrix_get_data_p_2d_s      => dbcsr_get_data_p_2d_s

 USE dbcsr_types, ONLY: &
      block_type_real_8 => dbcsr_type_real_8, &
      block_type_real_4 => dbcsr_type_real_4, &
      block_type_complex_8 => dbcsr_type_complex_8, &
      block_type_complex_4 => dbcsr_type_complex_4, &
      block_scalar_type => dbcsr_scalar_type, &
      block_matrix_transpose => dbcsr_transpose, &
      block_matrix_no_transpose => dbcsr_no_transpose, &
      block_matrix_type_no_symmetry => dbcsr_type_no_symmetry, &
      block_matrix_distribution_obj => dbcsr_distribution_obj, &
      block_matrix_iterator => dbcsr_iterator, &
      block_mp_obj => dbcsr_mp_obj

 USE dbcsr_multiply_api, ONLY: block_matrix_multiply => dbcsr_multiply

 USE dbcsr_work_operations, ONLY: &
      block_matrix_create => dbcsr_create, &
      block_matrix_finalize => dbcsr_finalize

 USE kinds, ONLY: &
      default_string_length, &
      int_8, &
      dp

 USE dbcsr_data_types, ONLY: &
      block_scalar_type => dbcsr_scalar_type, &
      block_matrix_data_obj             => dbcsr_data_obj

 USE dbcsr_mp_methods, ONLY: &
      block_mp_new => dbcsr_mp_new

 USE dbcsr_index_operations, ONLY: &
      block_matrix_sort_indices => dbcsr_sort_indices

 USE dbcsr_config, ONLY: block_matrix_cfg => dbcsr_cfg

 USE dbcsr_toollib, ONLY: &
      swap

 USE dbcsr_dist_util, ONLY: &
      block_matrix_checksum => dbcsr_checksum

 USE dbcsr_blas_operations, ONLY: &
     block_matrix_lapack_larnv => dbcsr_lapack_larnv, &
     block_matrix_set_larnv_seed => set_larnv_seed


#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_matrix_api'

   PUBLIC :: block_matrix_type
   PUBLIC :: block_matrix_iterator
   PUBLIC :: block_matrix_iterator_blocks_left
   PUBLIC :: block_matrix_iterator_next_block
   PUBLIC :: block_matrix_iterator_start
   PUBLIC :: block_matrix_iterator_stop
   PUBLIC :: block_matrix_transpose
   PUBLIC :: block_matrix_no_transpose
   PUBLIC :: block_matrix_put_block
   PUBLIC :: block_matrix_clear
   PUBLIC :: block_matrix_release
   PUBLIC :: block_scalar
   PUBLIC :: block_scalar_type
   PUBLIC :: block_type_real_8
   PUBLIC :: block_type_real_4
   PUBLIC :: block_type_complex_8
   PUBLIC :: block_type_complex_4
   PUBLIC :: block_matrix_distribution_get
   PUBLIC :: block_matrix_distribution_obj
   PUBLIC :: block_matrix_get_info
   PUBLIC :: block_matrix_scale
   PUBLIC :: block_mp_new
   PUBLIC :: block_mp_obj
   PUBLIC :: block_matrix_sort_indices
   PUBLIC :: block_matrix_get_block_p
   PUBLIC :: block_matrix_put_block
   PUBLIC :: block_matrix_reserve_blocks
   PUBLIC :: block_matrix_data_new
   PUBLIC :: block_matrix_data_release
   PUBLIC :: block_matrix_type_1d_to_2d
   PUBLIC :: block_matrix_data_clear_pointer
   PUBLIC :: block_matrix_data_init
   PUBLIC :: block_matrix_data_obj
   PUBLIC :: block_matrix_distribution_col_dist
   PUBLIC :: block_matrix_distribution_new
   PUBLIC :: block_matrix_distribution_row_dist
   PUBLIC :: block_matrix_distribution_hold
   PUBLIC :: block_matrix_distribution_release
   PUBLIC :: block_matrix_get_data_type
   PUBLIC :: block_matrix_mp_release
   PUBLIC :: block_matrix_nblkcols_total
   PUBLIC :: block_matrix_nblkrows_total
   PUBLIC :: block_matrix_nfullrows_total
   PUBLIC :: block_matrix_nfullcols_total
   PUBLIC :: block_matrix_get_data_size
   PUBLIC :: block_matrix_get_num_blocks
   PUBLIC :: block_matrix_get_nze
   PUBLIC :: block_matrix_set
   PUBLIC :: block_matrix_filter
   PUBLIC :: block_matrix_create
   PUBLIC :: block_matrix_finalize
   PUBLIC :: block_matrix_type_no_symmetry
   PUBLIC :: block_matrix_data_get_sizes
   PUBLIC :: block_matrix_data_set_pointer
   PUBLIC :: block_matrix_get_data_p_2d_d
   PUBLIC :: block_matrix_get_data_p_2d_z
   PUBLIC :: block_matrix_get_data_p_2d_c
   PUBLIC :: block_matrix_get_data_p_2d_s
   PUBLIC :: internal_data_allocate
   PUBLIC :: internal_data_deallocate
   PUBLIC :: block_matrix_get_stored_coordinates
   PUBLIC :: block_matrix_blk_column_size
   PUBLIC :: block_matrix_blk_row_size
   PUBLIC :: block_scalar_zero
   PUBLIC :: block_matrix_scalar_multiply
   PUBLIC :: block_matrix_multiply
   PUBLIC :: block_matrix_redistribute
   PUBLIC :: block_mp_release
   PUBLIC :: block_matrix_cfg
   PUBLIC :: block_matrix_copy
   PUBLIC :: block_matrix_add
   PUBLIC :: block_matrix_maxabs
   PUBLIC :: block_matrix_zero
   PUBLIC :: block_matrix_get_matrix_type
   PUBLIC :: checker_tr
   PUBLIC :: swap
   PUBLIC :: block_matrix_row_block_sizes
   PUBLIC :: block_matrix_col_block_sizes
   PUBLIC :: block_matrix_checksum
   PUBLIC :: block_matrix_complete_redistribute
   PUBLIC :: block_matrix_lapack_larnv
   PUBLIC :: block_matrix_set_larnv_seed

   CONTAINS

END MODULE block_matrix_api

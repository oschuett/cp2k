!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE dbs_tensor_api
   !! This is the start of a dbs_tensor_api, all publically needed functions
   !! are exported here. The others remain private to the library.
   !! Currently, this is the CP2K used set.
   !! Ultimately, a reduced subset and well defined api will remain,
   !! possibly grouped in to standard and expert api.
   !! Currently, this is work in progress.

   USE dbs_tensor, ONLY: &
      dbs_t_contract, dbs_t_get_block, dbs_t_get_stored_coordinates, dbs_t_put_block, &
      dbs_t_reserve_blocks, dbs_t_copy_matrix_to_tensor, dbs_t_copy, &
      dbs_t_copy_tensor_to_matrix, dbs_t_batched_contract_init, &
      dbs_t_batched_contract_finalize, &
      dbs_t_contract_index
   USE dbs_tensor_block, ONLY: &
      dbs_t_iterator_blocks_left, dbs_t_iterator_next_block, dbs_t_iterator_start, &
      dbs_t_iterator_stop, dbs_t_iterator_type, dbs_t_reserved_block_indices
   USE dbs_tensor_types, ONLY: &
      dbs_t_create, dbs_t_destroy, dbs_t_distribution_destroy, dbs_t_distribution_new, &
      dbs_t_distribution_type, dbs_t_nd_mp_comm, dbs_t_nd_mp_free, dbs_t_type, &
      dbs_t_pgrid_type, dbs_t_pgrid_create, dbs_t_pgrid_create_expert, dbs_t_pgrid_destroy, &
      dbs_t_set, dbs_t_filter, &
      dbs_t_mp_environ_pgrid => mp_environ_pgrid, dbs_t_blk_sizes, dbs_t_get_info, &
      dbs_t_finalize, dbs_t_scale, dbs_t_get_nze, dbs_t_get_nze_total, &
      dbs_t_get_num_blocks, dbs_t_get_num_blocks_total, dbs_t_clear, &
      dbs_t_mp_dims_create, dbs_t_pgrid_change_dims, dbs_t_ndims => ndims_tensor, &
      dbs_t_dims => dims_tensor, dbs_t_ndims_matrix_row => ndims_matrix_row, &
      dbs_t_ndims_matrix_column => ndims_matrix_column, dbs_t_blk_size, dbs_t_nblks_local, &
      dbs_t_nblks_total, dbs_t_max_nblks_local, dbs_t_default_distvec
   USE dbs_tensor_test, ONLY: &
      dbs_t_contract_test, dbs_t_checksum
   USE dbs_tensor_split, ONLY: &
      dbs_t_split_blocks
   USE dbs_tensor_index, ONLY: &
      dbs_t_get_mapping_info
   USE dbs_tensor_io, ONLY: &
      dbs_t_write_tensor_info, dbs_t_write_split_info, dbs_t_write_blocks, dbs_t_write_tensor_dist

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbs_tensor_api'

   PUBLIC :: dbs_t_contract, dbs_t_contract_test
   PUBLIC :: dbs_t_get_block
   PUBLIC :: dbs_t_get_stored_coordinates
   PUBLIC :: dbs_t_put_block
   PUBLIC :: dbs_t_reserve_blocks
   PUBLIC :: dbs_t_create
   PUBLIC :: dbs_t_destroy
   PUBLIC :: dbs_t_distribution_destroy
   PUBLIC :: dbs_t_distribution_new
   PUBLIC :: dbs_t_distribution_type
   PUBLIC :: dbs_t_nd_mp_comm
   PUBLIC :: dbs_t_nd_mp_free
   PUBLIC :: dbs_t_type
   PUBLIC :: dbs_t_iterator_next_block
   PUBLIC :: dbs_t_iterator_blocks_left
   PUBLIC :: dbs_t_iterator_stop
   PUBLIC :: dbs_t_iterator_start
   PUBLIC :: dbs_t_iterator_type
   PUBLIC :: dbs_t_split_blocks
   PUBLIC :: dbs_t_pgrid_type
   PUBLIC :: dbs_t_pgrid_create
   PUBLIC :: dbs_t_pgrid_create_expert
   PUBLIC :: dbs_t_pgrid_destroy
   PUBLIC :: dbs_t_set
   PUBLIC :: dbs_t_filter
   PUBLIC :: dbs_t_mp_environ_pgrid
   PUBLIC :: dbs_t_copy_matrix_to_tensor
   PUBLIC :: dbs_t_blk_sizes
   PUBLIC :: dbs_t_copy
   PUBLIC :: dbs_t_copy_tensor_to_matrix
   PUBLIC :: dbs_t_get_info
   PUBLIC :: dbs_t_checksum
   PUBLIC :: dbs_t_finalize
   PUBLIC :: dbs_t_scale
   PUBLIC :: dbs_t_get_num_blocks, dbs_t_get_num_blocks_total
   PUBLIC :: dbs_t_get_nze, dbs_t_get_nze_total
   PUBLIC :: dbs_t_clear
   PUBLIC :: dbs_t_get_mapping_info
   PUBLIC :: dbs_t_write_split_info
   PUBLIC :: dbs_t_write_blocks
   PUBLIC :: dbs_t_write_tensor_dist
   PUBLIC :: dbs_t_write_tensor_info
   PUBLIC :: dbs_t_mp_dims_create
   PUBLIC :: dbs_t_batched_contract_init
   PUBLIC :: dbs_t_batched_contract_finalize
   PUBLIC :: dbs_t_ndims
   PUBLIC :: dbs_t_dims
   PUBLIC :: dbs_t_pgrid_change_dims
   PUBLIC :: dbs_t_reserved_block_indices
   PUBLIC :: dbs_t_contract_index
   PUBLIC :: dbs_t_ndims_matrix_row
   PUBLIC :: dbs_t_ndims_matrix_column
   PUBLIC :: dbs_t_nblks_local
   PUBLIC :: dbs_t_nblks_total
   PUBLIC :: dbs_t_blk_size
   PUBLIC :: dbs_t_max_nblks_local
   PUBLIC :: dbs_t_default_distvec

END MODULE dbs_tensor_api

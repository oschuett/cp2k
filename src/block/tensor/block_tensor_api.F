!--------------------------------------------------------------------------------------------------!
! Copyright (C) by the DBCSR developers group - All rights reserved                                !
! This file is part of the DBCSR library.                                                          !
!                                                                                                  !
! For information on the license, see the LICENSE file.                                            !
! For further information please visit https://dbcsr.cp2k.org                                      !
! SPDX-License-Identifier: GPL-2.0+                                                                !
!--------------------------------------------------------------------------------------------------!

MODULE block_tensor_api
   !! This is the start of a block_tensor_api, all publically needed functions
   !! are exported here. The others remain private to the library.
   !! Currently, this is the CP2K used set.
   !! Ultimately, a reduced subset and well defined api will remain,
   !! possibly grouped in to standard and expert api.
   !! Currently, this is work in progress.

   USE block_tensor, ONLY: &
      block_tensor_contract, block_tensor_get_block, block_tensor_get_stored_coordinates, block_tensor_put_block, &
      block_tensor_reserve_blocks, block_tensor_copy_matrix_to_tensor, block_tensor_copy, &
      block_tensor_copy_tensor_to_matrix, block_tensor_batched_contract_init, &
      block_tensor_batched_contract_finalize, &
      block_tensor_contract_index
   USE block_tensor_block, ONLY: &
      block_tensor_iterator_blocks_left, block_tensor_iterator_next_block, block_tensor_iterator_start, &
      block_tensor_iterator_stop, block_tensor_iterator_type, block_tensor_reserved_block_indices
   USE block_tensor_types, ONLY: &
      block_tensor_create, block_tensor_destroy, block_tensor_distribution_destroy, block_tensor_distribution_new, &
      block_tensor_distribution_type, block_tensor_nd_mp_comm, block_tensor_nd_mp_free, block_tensor_type, &
      block_tensor_pgrid_type, block_tensor_pgrid_create, block_tensor_pgrid_create_expert, block_tensor_pgrid_destroy, &
      block_tensor_set, block_tensor_filter, &
      block_tensor_mp_environ_pgrid => mp_environ_pgrid, block_tensor_blk_sizes, block_tensor_get_info, &
      block_tensor_finalize, block_tensor_scale, block_tensor_get_nze, block_tensor_get_nze_total, &
      block_tensor_get_num_blocks, block_tensor_get_num_blocks_total, block_tensor_clear, &
      block_tensor_mp_dims_create, block_tensor_pgrid_change_dims, block_tensor_ndims => ndims_tensor, &
      block_tensor_dims => dims_tensor, block_tensor_ndims_matrix_row => ndims_matrix_row, &
      block_tensor_ndims_matrix_column => ndims_matrix_column, block_tensor_blk_size, block_tensor_nblks_local, &
      block_tensor_nblks_total, block_tensor_max_nblks_local, block_tensor_default_distvec
   USE block_tensor_test, ONLY: &
      block_tensor_contract_test, block_tensor_checksum
   USE block_tensor_split, ONLY: &
      block_tensor_split_blocks
   USE block_tensor_index, ONLY: &
      block_tensor_get_mapping_info
   USE block_tensor_io, ONLY: &
      block_tensor_write_tensor_info, block_tensor_write_split_info, block_tensor_write_blocks, block_tensor_write_tensor_dist

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'block_tensor_api'

   PUBLIC :: block_tensor_contract, block_tensor_contract_test
   PUBLIC :: block_tensor_get_block
   PUBLIC :: block_tensor_get_stored_coordinates
   PUBLIC :: block_tensor_put_block
   PUBLIC :: block_tensor_reserve_blocks
   PUBLIC :: block_tensor_create
   PUBLIC :: block_tensor_destroy
   PUBLIC :: block_tensor_distribution_destroy
   PUBLIC :: block_tensor_distribution_new
   PUBLIC :: block_tensor_distribution_type
   PUBLIC :: block_tensor_nd_mp_comm
   PUBLIC :: block_tensor_nd_mp_free
   PUBLIC :: block_tensor_type
   PUBLIC :: block_tensor_iterator_next_block
   PUBLIC :: block_tensor_iterator_blocks_left
   PUBLIC :: block_tensor_iterator_stop
   PUBLIC :: block_tensor_iterator_start
   PUBLIC :: block_tensor_iterator_type
   PUBLIC :: block_tensor_split_blocks
   PUBLIC :: block_tensor_pgrid_type
   PUBLIC :: block_tensor_pgrid_create
   PUBLIC :: block_tensor_pgrid_create_expert
   PUBLIC :: block_tensor_pgrid_destroy
   PUBLIC :: block_tensor_set
   PUBLIC :: block_tensor_filter
   PUBLIC :: block_tensor_mp_environ_pgrid
   PUBLIC :: block_tensor_copy_matrix_to_tensor
   PUBLIC :: block_tensor_blk_sizes
   PUBLIC :: block_tensor_copy
   PUBLIC :: block_tensor_copy_tensor_to_matrix
   PUBLIC :: block_tensor_get_info
   PUBLIC :: block_tensor_checksum
   PUBLIC :: block_tensor_finalize
   PUBLIC :: block_tensor_scale
   PUBLIC :: block_tensor_get_num_blocks, block_tensor_get_num_blocks_total
   PUBLIC :: block_tensor_get_nze, block_tensor_get_nze_total
   PUBLIC :: block_tensor_clear
   PUBLIC :: block_tensor_get_mapping_info
   PUBLIC :: block_tensor_write_split_info
   PUBLIC :: block_tensor_write_blocks
   PUBLIC :: block_tensor_write_tensor_dist
   PUBLIC :: block_tensor_write_tensor_info
   PUBLIC :: block_tensor_mp_dims_create
   PUBLIC :: block_tensor_batched_contract_init
   PUBLIC :: block_tensor_batched_contract_finalize
   PUBLIC :: block_tensor_ndims
   PUBLIC :: block_tensor_dims
   PUBLIC :: block_tensor_pgrid_change_dims
   PUBLIC :: block_tensor_reserved_block_indices
   PUBLIC :: block_tensor_contract_index
   PUBLIC :: block_tensor_ndims_matrix_row
   PUBLIC :: block_tensor_ndims_matrix_column
   PUBLIC :: block_tensor_nblks_local
   PUBLIC :: block_tensor_nblks_total
   PUBLIC :: block_tensor_blk_size
   PUBLIC :: block_tensor_max_nblks_local
   PUBLIC :: block_tensor_default_distvec

END MODULE block_tensor_api

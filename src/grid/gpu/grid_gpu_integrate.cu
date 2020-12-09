/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <algorithm>
#include <assert.h>
#include <cuda.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_DO_COLLOCATE 0
#include "../common/grid_common.h"
#include "grid_gpu_collint.h"
#include "grid_gpu_integrate.h"

/*******************************************************************************
 * \brief Decontracts the subblock, going from spherical to cartesian harmonics.
 * \author Ole Schuett
 ******************************************************************************/
__device__ static void store_hab(const kernel_params *params,
                                 const smem_task *task, const double *cab) {

  // The spherical index runs over angular momentum and then over contractions.
  // The carthesian index runs over exponents and then over angular momentum.

  // This is a double matrix product. Since the block can be quite large the
  // two products are fused to conserve shared memory.
  for (int i = threadIdx.x; i < task->nsgf_setb; i += blockDim.x) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      for (int k = threadIdx.z; k < task->ncosetb; k += blockDim.z) {
        double block_val = 0.0;
        const double sphib = task->sphib[i * task->maxcob + k];
        for (int l = 0; l < task->ncoseta; l++) {
          const double sphia = task->sphia[j * task->maxcoa + l];
          block_val += cab[k * task->ncoseta + l] * sphia * sphib;
        }
        if (task->block_transposed) {
          atomicAddDouble(&task->hab_block[j * task->nsgfb + i], block_val);
        } else {
          atomicAddDouble(&task->hab_block[i * task->nsgfa + j], block_val);
        }
      }
    }
  }
  __syncthreads(); // TODO: not really neded because of concurrent writes to
  // cab
}

/*******************************************************************************
 * \brief Cuda kernel for integrating all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void integrate_kernel(const kernel_params params) {

  // Copy task from global to shared memory and precompute some stuff.
  __shared__ smem_task task;
  fill_smem_task(&params, &task);

  // Check if radius is below the resolution of the grid.
  if (2.0 * task.radius < task.dh_max) {
    return; // nothing to do
  }

  // Allot dynamic shared memory.
  extern __shared__ double shared_memory[];
  double *smem_cab = &shared_memory[params.smem_cab_offset];
  double *smem_alpha = &shared_memory[params.smem_alpha_offset];
  double *smem_cxyz = &shared_memory[params.smem_cxyz_offset];

  memset(smem_cxyz, 0, ncoset(task.lp) * sizeof(double));
  __syncthreads();

  cxyz_to_grid(&params, &task, smem_cxyz, params.grid);

  memset(smem_cab, 0, task.n1_cab * task.n2_cab * sizeof(double));
  __syncthreads();

  compute_alpha(&params, &task, smem_alpha);
  cab_to_cxyz(&params, &task, smem_alpha, smem_cab, smem_cxyz);

  store_hab(&params, &task, smem_cab);

  // if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
  //  printf("la_min: %i %lb_min: %i, ncoset: %i, %i \n",
  //      task.la_min, task.lb_min,
  //      ncoset(task.la_min-1), ncoset(task.lb_min-1));
  //  //    for (int k = 0; k < task.ncosetb; k++) {
  //  //      for (int l = 0; l < task.ncoseta; l++) {
  //  //          printf("cab %i %i %le\n", k, l, smem_cab[k * task.ncoseta +
  //  l]);
  //  //      }
  //  //    }
  //  // printf("cxyz %i %i %le\n",0, 0, smem_cxyz[0]);
  //}
}

/*******************************************************************************
 * \brief Launches the Cuda kernel that integrates all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_one_grid_level(
    const grid_gpu_task_list *task_list, const int first_task,
    const int last_task, const bool orthorhombic, const bool compute_tau,
    const bool calculate_forces, const int npts_global[3],
    const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double dh[3][3], const double dh_inv[3][3],
    const cudaStream_t stream, const double *pab_blocks_dev,
    const double *grid_dev, double *hab_blocks_dev, double *forces_dev,
    double *virial_dev) {

  const int ntasks = last_task - first_task + 1;
  if (ntasks == 0) {
    return; // Nothing to do.
  }

  init_constant_memory();

  // Compute required shared memory.
  // TODO: Currently, cab's indicies run over 0...ncoset[lmax],
  //       however only ncoset(lmin)...ncoset(lmax) are actually needed.
  const int lmax = task_list->lmax;
  const int lp_max = 2 * lmax;
  const int cab_len = ncoset(lmax) * ncoset(lmax);
  const int alpha_len = 3 * (lmax + 1) * (lmax + 1) * (lp_max + 1);
  const int cxyz_len = ncoset(lp_max);
  const size_t smem_per_block =
      (cab_len + alpha_len + cxyz_len) * sizeof(double);

  if (smem_per_block > 48 * 1024) {
    fprintf(stderr, "ERROR: Not enough shared memory.\n");
    fprintf(stderr, "cab_len: %i, ", cab_len);
    fprintf(stderr, "alpha_len: %i, ", alpha_len);
    fprintf(stderr, "cxyz_len: %i, ", cxyz_len);
    fprintf(stderr, "total smem_per_block: %f kb\n\n", smem_per_block / 1024.0);
    abort();
  }

  // kernel parameters
  kernel_params params;
  params.smem_cab_offset = 0;
  params.smem_alpha_offset = cab_len;
  params.smem_cxyz_offset = params.smem_alpha_offset + alpha_len;
  params.first_task = first_task;
  params.orthorhombic = orthorhombic;
  params.compute_tau = compute_tau;
  params.calculate_forces = calculate_forces;
  params.grid = grid_dev;
  params.tasks = task_list->tasks_dev;
  params.atom_kinds = task_list->atom_kinds_dev;
  params.basis_sets = task_list->basis_sets_dev;
  params.block_offsets = task_list->block_offsets_dev;
  params.atom_positions = task_list->atom_positions_dev;
  params.pab_blocks = pab_blocks_dev;
  params.hab_blocks = hab_blocks_dev;
  params.forces = forces_dev;
  params.virial = virial_dev;
  memcpy(params.dh, dh, 9 * sizeof(double));
  memcpy(params.dh_inv, dh_inv, 9 * sizeof(double));
  memcpy(params.npts_global, npts_global, 3 * sizeof(int));
  memcpy(params.npts_local, npts_local, 3 * sizeof(int));
  memcpy(params.shift_local, shift_local, 3 * sizeof(int));
  memcpy(params.border_width, border_width, 3 * sizeof(int));

  // Launch !
  const int nblocks = ntasks;
  const dim3 threads_per_block(4, 8, 8);

  integrate_kernel<<<nblocks, threads_per_block, smem_per_block, stream>>>(
      params);
}

#endif // __GRID_CUDA
// EOF

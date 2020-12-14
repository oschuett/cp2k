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
#include "../common/grid_process_vab.h"
#include "grid_gpu_collint.h"
#include "grid_gpu_integrate.h"

/*******************************************************************************
 * \brief Decontracts the subblock, going from spherical to cartesian harmonics.
 * \author Ole Schuett
 ******************************************************************************/
template <bool COMPUTE_TAU>
__device__ static void store_hab(const kernel_params *params,
                                 const smem_task *task, const double *cab) {

  // The spherical index runs over angular momentum and then over contractions.
  // The carthesian index runs over exponents and then over angular momentum.

  // This is a double matrix product. Since the block can be quite large the
  // two products are fused to conserve shared memory.
  // TODO: move into smem_task
  const int ico_start =
      (task->la_min_basis > 0) ? ncoset(task->la_min_basis - 1) : 0;
  const int jco_start =
      (task->lb_min_basis > 0) ? ncoset(task->lb_min_basis - 1) : 0;

  for (int i = threadIdx.x; i < task->nsgf_setb; i += blockDim.x) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      for (int jco = jco_start + threadIdx.z; jco < ncoset(task->lb_max_basis);
           jco += blockDim.z) {
        const orbital b = coset_inv[jco];
        double block_val = 0.0;
        const double sphib = task->sphib[i * task->maxcob + jco];
        for (int ico = ico_start; ico < ncoset(task->la_max_basis); ico++) {
          const orbital a = coset_inv[ico];
          double hab = 0.0;
          if (COMPUTE_TAU) {
            // Since process_tau is a register hog we use it only when needed.
            hab = extract_hab_tau(a, b, task->zeta, task->zetb, task->n1, cab);
          } else {
            // fast path for common case
            hab = extract_hab(a, b, task->n1, cab);
          }

          const double sphia = task->sphia[j * task->maxcoa + ico];
          block_val += hab * sphia * sphib;
        }
        if (task->block_transposed) {
          atomicAddDouble(&task->hab_block[j * task->nsgfb + i], block_val);
        } else {
          atomicAddDouble(&task->hab_block[i * task->nsgfa + j], block_val);
        }
      }
    }
  }
  __syncthreads(); // TODO: Probably not needed
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
template <bool COMPUTE_TAU>
__device__ static void store_forces_and_virial(const kernel_params *params,
                                               const smem_task *task,
                                               const double *cab) {

  // The spherical index runs over angular momentum and then over contractions.
  // The carthesian index runs over exponents and then over angular momentum.

  // Decontract block, apply prepare_pab, and store in cab.
  // This is a double matrix product. Since the pab block can be quite large the
  // two products are fused to conserve shared memory.

  // TODO: move into smem_task
  const int ico_start =
      (task->la_min_basis > 0) ? ncoset(task->la_min_basis - 1) : 0;
  const int jco_start =
      (task->lb_min_basis > 0) ? ncoset(task->lb_min_basis - 1) : 0;

  for (int i = threadIdx.x; i < task->nsgf_setb; i += blockDim.x) {
    for (int j = threadIdx.y; j < task->nsgf_seta; j += blockDim.y) {
      double block_val;
      if (task->block_transposed) {
        block_val = task->pab_block[j * task->nsgfb + i];
      } else {
        block_val = task->pab_block[i * task->nsgfa + j];
      }
      for (int jco = jco_start + threadIdx.z; jco < task->ncosetb;
           jco += blockDim.z) {
        const double sphib = task->sphib[i * task->maxcob + jco];
        for (int ico = ico_start; ico < task->ncoseta; ico++) {
          const double sphia = task->sphia[j * task->maxcoa + ico];
          const double pabval = block_val * sphia * sphib;
          const orbital b = coset_inv[jco];
          const orbital a = coset_inv[ico];
          for (int k = 0; k < 3; k++) {
            const double force_a =
                pabval * task->off_diagonals_twice *
                extract_force_a(a, b, k, task->zeta, task->n1, cab);
            atomicAddDouble(&task->forces_a[k], force_a);
            const double force_b =
                pabval * task->off_diagonals_twice *
                extract_force_b(a, b, k, task->zetb, task->rab, task->n1, cab);
            atomicAddDouble(&task->forces_b[k], force_b);
          }
          if (params->virial != NULL) {
            for (int k = 0; k < 3; k++) {
              for (int l = 0; l < 3; l++) {
                const double virial_a =
                    extract_virial_a(a, b, k, l, task->zeta, task->n1, cab);
                const double virial_b = extract_virial_b(
                    a, b, k, l, task->zetb, task->rab, task->n1, cab);
                const double virial =
                    pabval * task->off_diagonals_twice * (virial_a + virial_b);
                atomicAddDouble(&params->virial[k * 3 + l], virial);
              }
            }
          }
        }
      }
    }
  }
  __syncthreads(); // TODO: Probably not needed
}

/*******************************************************************************
 * \brief Cuda kernel for integrating all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
template <bool COMPUTE_TAU, bool CALCULATE_FORCES>
__device__ static void integrate_kernel(const kernel_params *params) {

  // Copy task from global to shared memory and precompute some stuff.
  __shared__ smem_task task;
  fill_smem_task(params, &task);

  // Check if radius is below the resolution of the grid.
  if (2.0 * task.radius < task.dh_max) {
    return; // nothing to do
  }

  // Allot dynamic shared memory.
  extern __shared__ double shared_memory[];
  double *smem_cab = &shared_memory[params->smem_cab_offset];
  double *smem_alpha = &shared_memory[params->smem_alpha_offset];
  double *smem_cxyz = &shared_memory[params->smem_cxyz_offset];

  memset(smem_cxyz, 0, ncoset(task.lp) * sizeof(double));
  __syncthreads();

  cxyz_to_grid(params, &task, smem_cxyz, params->grid);

  memset(smem_cab, 0, task.n1 * task.n2 * sizeof(double));
  __syncthreads();

  compute_alpha(params, &task, smem_alpha);
  cab_to_cxyz(params, &task, smem_alpha, smem_cab, smem_cxyz);

  store_hab<COMPUTE_TAU>(params, &task, smem_cab);

  if (CALCULATE_FORCES) {
    store_forces_and_virial<COMPUTE_TAU>(params, &task, smem_cab);
  }

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
 * \brief Specialized Cuda kernel for compute_tau=false & calculate_forces=false
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void integrate_density(const kernel_params params) {
  integrate_kernel<false, false>(&params);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=true & calculate_forces=false.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void integrate_tau(const kernel_params params) {
  integrate_kernel<true, false>(&params);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=false & calculate_forces=true.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void integrate_density_forces(const kernel_params params) {
  integrate_kernel<false, true>(&params);
}

/*******************************************************************************
 * \brief Specialized Cuda kernel for compute_tau=true & calculate_forces=true.
 * \author Ole Schuett
 ******************************************************************************/
__global__ static void integrate_tau_forces(const kernel_params params) {
  integrate_kernel<true, true>(&params);
}

/*******************************************************************************
 * \brief Launches the Cuda kernel that integrates all tasks of one grid level.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_one_grid_level(
    const grid_gpu_task_list *task_list, const int first_task,
    const int last_task, const bool orthorhombic, const bool compute_tau,
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
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
  const bool calculate_forces = forces_dev != NULL;
  const bool calculate_virial = virial_dev != NULL;
  const process_ldiffs ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);
  const int la_max = task_list->lmax + ldiffs.la_max_diff;
  const int lb_max = task_list->lmax + ldiffs.lb_max_diff;
  const int lp_max = la_max + lb_max;
  const int cab_len = ncoset(lb_max) * ncoset(la_max);
  const int alpha_len = 3 * (lb_max + 1) * (la_max + 1) * (lp_max + 1);
  const int cxyz_len = ncoset(lp_max);
  const size_t smem_per_block =
      (cab_len + alpha_len + cxyz_len) * sizeof(double);

  if (smem_per_block > 48 * 1024) {
    fprintf(stderr, "ERROR: Not enough shared memory in grid_gpu_integrate.\n");
    fprintf(stderr, "cab_len: %i, ", cab_len);
    fprintf(stderr, "alpha_len: %i, ", alpha_len);
    fprintf(stderr, "cxyz_len: %i, ", cxyz_len);
    fprintf(stderr, "total smem_per_block: %f kb\n\n", smem_per_block / 1024.0);
    abort();
  }

  // assert(compute_tau == false);

  // kernel parameters
  kernel_params params;
  params.smem_cab_offset = 0;
  params.smem_alpha_offset = cab_len;
  params.smem_cxyz_offset = params.smem_alpha_offset + alpha_len;
  params.first_task = first_task;
  params.orthorhombic = orthorhombic;
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
  params.la_min_diff = ldiffs.la_min_diff;
  params.lb_min_diff = ldiffs.lb_min_diff;
  params.la_max_diff = ldiffs.la_max_diff;
  params.lb_max_diff = ldiffs.lb_max_diff;
  memcpy(params.dh, dh, 9 * sizeof(double));
  memcpy(params.dh_inv, dh_inv, 9 * sizeof(double));
  memcpy(params.npts_global, npts_global, 3 * sizeof(int));
  memcpy(params.npts_local, npts_local, 3 * sizeof(int));
  memcpy(params.shift_local, shift_local, 3 * sizeof(int));
  memcpy(params.border_width, border_width, 3 * sizeof(int));

  // Launch !
  const int nblocks = ntasks;
  const dim3 threads_per_block(4, 8, 8);

  if (!compute_tau && !calculate_forces) {
    integrate_density<<<nblocks, threads_per_block, smem_per_block, stream>>>(
        params);
  } else if (compute_tau && !calculate_forces) {
    integrate_tau<<<nblocks, threads_per_block, smem_per_block, stream>>>(
        params);
  } else if (!compute_tau && calculate_forces) {
    integrate_density_forces<<<nblocks, threads_per_block, smem_per_block,
                               stream>>>(params);
  } else if (compute_tau && calculate_forces) {
    integrate_tau_forces<<<nblocks, threads_per_block, smem_per_block,
                           stream>>>(params);
  }
}

#endif // __GRID_CUDA
// EOF

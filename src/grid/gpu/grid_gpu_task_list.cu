/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#ifdef __GRID_CUDA

#include <cuda_runtime.h>

#include <assert.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "../common/grid_constants.h"
#include "grid_gpu_collocate.h"
#include "grid_gpu_integrate.h"
#include "grid_gpu_task_list.h"

/*******************************************************************************
 * \brief Check given Cuda status and upon failure abort with a nice message.
 * \author Ole Schuett
 ******************************************************************************/
#define CHECK(status)                                                          \
  if (status != cudaSuccess) {                                                 \
    fprintf(stderr, "ERROR: %s %s %d\n", cudaGetErrorString(status), __FILE__, \
            __LINE__);                                                         \
    abort();                                                                   \
  }

/*******************************************************************************
 * \brief Allocates a task list for the GPU backend.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_create_task_list(
    const int ntasks, const int nlevels, const int natoms, const int nkinds,
    const int nblocks, const int block_offsets[],
    const double atom_positions[][3], const int atom_kinds[],
    const grid_basis_set *basis_sets[], const int level_list[],
    const int iatom_list[], const int jatom_list[], const int iset_list[],
    const int jset_list[], const int ipgf_list[], const int jpgf_list[],
    const int border_mask_list[], const int block_num_list[],
    const double radius_list[], const double rab_list[][3],
    grid_gpu_task_list **task_list_out) {

  if (*task_list_out != NULL) {
    // This is actually an opportunity to reuse some buffers.
    grid_gpu_free_task_list(*task_list_out);
  }

  grid_gpu_task_list *task_list =
      (grid_gpu_task_list *)malloc(sizeof(grid_gpu_task_list));

  task_list->ntasks = ntasks;
  task_list->nlevels = nlevels;
  task_list->natoms = natoms;
  task_list->nkinds = nkinds;
  task_list->nblocks = nblocks;

  size_t size = nblocks * sizeof(int);
  CHECK(cudaMalloc(&task_list->block_offsets_dev, size));
  CHECK(cudaMemcpy(task_list->block_offsets_dev, block_offsets, size,
                   cudaMemcpyHostToDevice));

  size = 3 * natoms * sizeof(double);
  CHECK(cudaMalloc(&task_list->atom_positions_dev, size));
  CHECK(cudaMemcpy(task_list->atom_positions_dev, atom_positions, size,
                   cudaMemcpyHostToDevice));

  size = natoms * sizeof(int);
  CHECK(cudaMalloc(&task_list->atom_kinds_dev, size));
  CHECK(cudaMemcpy(task_list->atom_kinds_dev, atom_kinds, size,
                   cudaMemcpyHostToDevice));

  // Upload basis sets to device.
  grid_basis_set basis_sets_host[nkinds];
  for (int i = 0; i < nkinds; i++) {
    const grid_basis_set *basis_set = basis_sets[i];
    grid_basis_set *basis_set_host = &basis_sets_host[i];
    basis_set_host->nset = basis_set->nset;
    basis_set_host->nsgf = basis_set->nsgf;
    basis_set_host->maxco = basis_set->maxco;
    basis_set_host->maxpgf = basis_set->maxpgf;

    size_t size = basis_set->nset * sizeof(int);
    CHECK(cudaMalloc(&basis_set_host->lmin, size));
    CHECK(cudaMemcpy(basis_set_host->lmin, basis_set->lmin, size,
                     cudaMemcpyHostToDevice));

    CHECK(cudaMalloc(&basis_set_host->lmax, size));
    CHECK(cudaMemcpy(basis_set_host->lmax, basis_set->lmax, size,
                     cudaMemcpyHostToDevice));

    CHECK(cudaMalloc(&basis_set_host->npgf, size));
    CHECK(cudaMemcpy(basis_set_host->npgf, basis_set->npgf, size,
                     cudaMemcpyHostToDevice));

    CHECK(cudaMalloc(&basis_set_host->nsgf_set, size));
    CHECK(cudaMemcpy(basis_set_host->nsgf_set, basis_set->nsgf_set, size,
                     cudaMemcpyHostToDevice));

    CHECK(cudaMalloc(&basis_set_host->first_sgf, size));
    CHECK(cudaMemcpy(basis_set_host->first_sgf, basis_set->first_sgf, size,
                     cudaMemcpyHostToDevice));

    size = basis_set->nsgf * basis_set->maxco * sizeof(double);
    CHECK(cudaMalloc(&basis_set_host->sphi, size));
    CHECK(cudaMemcpy(basis_set_host->sphi, basis_set->sphi, size,
                     cudaMemcpyHostToDevice));

    size = basis_set->nset * basis_set->maxpgf * sizeof(double);
    CHECK(cudaMalloc(&basis_set_host->zet, size));
    CHECK(cudaMemcpy(basis_set_host->zet, basis_set->zet, size,
                     cudaMemcpyHostToDevice));
  }
  size = nkinds * sizeof(grid_basis_set);
  CHECK(cudaMalloc(&task_list->basis_sets_dev, size));
  CHECK(cudaMemcpy(task_list->basis_sets_dev, basis_sets_host, size,
                   cudaMemcpyHostToDevice));

  size = ntasks * sizeof(grid_gpu_task);
  grid_gpu_task *tasks_host = (grid_gpu_task *)malloc(size);
  for (int i = 0; i < ntasks; i++) {
    tasks_host[i].level = level_list[i];
    tasks_host[i].iatom = iatom_list[i];
    tasks_host[i].jatom = jatom_list[i];
    tasks_host[i].iset = iset_list[i];
    tasks_host[i].jset = jset_list[i];
    tasks_host[i].ipgf = ipgf_list[i];
    tasks_host[i].jpgf = jpgf_list[i];
    tasks_host[i].border_mask = border_mask_list[i];
    tasks_host[i].block_num = block_num_list[i];
    tasks_host[i].radius = radius_list[i];
    tasks_host[i].rab[0] = rab_list[i][0];
    tasks_host[i].rab[1] = rab_list[i][1];
    tasks_host[i].rab[2] = rab_list[i][2];
  }
  CHECK(cudaMalloc(&task_list->tasks_dev, size));
  CHECK(cudaMemcpy(task_list->tasks_dev, tasks_host, size,
                   cudaMemcpyHostToDevice));
  free(tasks_host);

  // Count tasks per level.
  size = nlevels * sizeof(int);
  task_list->tasks_per_level = (int *)malloc(size);
  memset(task_list->tasks_per_level, 0, size);
  for (int i = 0; i < ntasks; i++) {
    task_list->tasks_per_level[level_list[i] - 1]++;
    assert(i == 0 || level_list[i] >= level_list[i - 1]); // expect ordered list
  }

  // Find largest angular momentum.
  task_list->lmax = 0;
  for (int ikind = 0; ikind < nkinds; ikind++) {
    for (int iset = 0; iset < basis_sets[ikind]->nset; iset++) {
      task_list->lmax = imax(task_list->lmax, basis_sets[ikind]->lmax[iset]);
    }
  }

  // allocate main cuda stream
  CHECK(cudaStreamCreate(&task_list->main_stream));

  // allocate one cuda stream per grid level
  size = nlevels * sizeof(cudaStream_t);
  task_list->level_streams = (cudaStream_t *)malloc(size);
  for (int i = 0; i < nlevels; i++) {
    CHECK(cudaStreamCreate(&task_list->level_streams[i]));
  }

  size = nlevels * sizeof(double *);
  task_list->grid_dev = (double **)malloc(size);
  memset(task_list->grid_dev, 0, size);

  size = nlevels * sizeof(size_t);
  task_list->grid_dev_size = (size_t *)malloc(size);
  memset(task_list->grid_dev_size, 0, size);

  // return newly created task list
  *task_list_out = task_list;
}

/*******************************************************************************
 * \brief Deallocates given task list, basis_sets have to be freed separately.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_free_task_list(grid_gpu_task_list *task_list) {

  // Download basis sets from device to get device pointers to their lists.
  const int nkinds = task_list->nkinds;
  grid_basis_set basis_sets_host[nkinds];
  size_t size = nkinds * sizeof(grid_basis_set);
  CHECK(cudaMemcpy(basis_sets_host, task_list->basis_sets_dev, size,
                   cudaMemcpyDeviceToHost));
  for (int i = 0; i < nkinds; i++) {
    CHECK(cudaFree(basis_sets_host[i].lmin));
    CHECK(cudaFree(basis_sets_host[i].lmax));
    CHECK(cudaFree(basis_sets_host[i].npgf));
    CHECK(cudaFree(basis_sets_host[i].nsgf_set));
    CHECK(cudaFree(basis_sets_host[i].first_sgf));
    CHECK(cudaFree(basis_sets_host[i].sphi));
    CHECK(cudaFree(basis_sets_host[i].zet));
  }
  CHECK(cudaFree(task_list->basis_sets_dev));

  CHECK(cudaFree(task_list->block_offsets_dev));
  CHECK(cudaFree(task_list->atom_positions_dev));
  CHECK(cudaFree(task_list->atom_kinds_dev));
  CHECK(cudaFree(task_list->tasks_dev));

  CHECK(cudaStreamDestroy(task_list->main_stream));

  for (int i = 0; i < task_list->nlevels; i++) {
    CHECK(cudaStreamDestroy(task_list->level_streams[i]));
  }
  free(task_list->level_streams);

  for (int i = 0; i < task_list->nlevels; i++) {
    if (task_list->grid_dev[i] != NULL) {
      CHECK(cudaFree(task_list->grid_dev[i]));
    }
  }
  free(task_list->grid_dev);

  free(task_list->grid_dev_size);
  free(task_list->tasks_per_level);
  free(task_list);
}

/*******************************************************************************
 * \brief Collocate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_collocate_task_list(
    const grid_gpu_task_list *task_list, const bool orthorhombic,
    const enum grid_func func, const int nlevels, const int npts_global[][3],
    const int npts_local[][3], const int shift_local[][3],
    const int border_width[][3], const double dh[][3][3],
    const double dh_inv[][3][3], const grid_buffer *pab_blocks,
    double *grid[]) {
  assert(task_list->nlevels == nlevels);

  // Upload blocks buffer using the main stream
  CHECK(cudaMemcpyAsync(pab_blocks->device_buffer, pab_blocks->host_buffer,
                        pab_blocks->size, cudaMemcpyHostToDevice,
                        task_list->main_stream));

  // record an event so the level streams can wait for the blocks to be uploaded
  cudaEvent_t input_ready_event;
  CHECK(cudaEventCreate(&input_ready_event));
  CHECK(cudaEventRecord(input_ready_event, task_list->main_stream));

  int first_task = 0;
  for (int level = 0; level < task_list->nlevels; level++) {
    const int last_task = first_task + task_list->tasks_per_level[level] - 1;
    const cudaStream_t level_stream = task_list->level_streams[level];
    const size_t grid_size = npts_local[level][0] * npts_local[level][1] *
                             npts_local[level][2] * sizeof(double);

    // reallocate device grid buffers if needed
    if (task_list->grid_dev_size[level] < grid_size) {
      if (task_list->grid_dev[level] != NULL) {
        CHECK(cudaFree(task_list->grid_dev[level]));
      }
      CHECK(cudaMalloc(&task_list->grid_dev[level], grid_size));
      task_list->grid_dev_size[level] = grid_size;
    }

    // zero device grid buffers
    CHECK(cudaMemsetAsync(task_list->grid_dev[level], 0, grid_size,
                          level_stream));

    // launch kernel, but only after blocks have arrived
    CHECK(cudaStreamWaitEvent(level_stream, input_ready_event, 0));
    grid_gpu_collocate_one_grid_level(
        task_list, first_task, last_task, orthorhombic, func,
        npts_global[level], npts_local[level], shift_local[level],
        border_width[level], dh[level], dh_inv[level], level_stream,
        pab_blocks->device_buffer, task_list->grid_dev[level]);

    first_task = last_task + 1;
  }

  // download result from device to host.
  // TODO: Make these mem copies actually async by page locking the grid buffers
  for (int level = 0; level < task_list->nlevels; level++) {
    const size_t grid_size = npts_local[level][0] * npts_local[level][1] *
                             npts_local[level][2] * sizeof(double);
    CHECK(cudaMemcpyAsync(grid[level], task_list->grid_dev[level], grid_size,
                          cudaMemcpyDeviceToHost,
                          task_list->level_streams[level]));
  }

  // clean up
  CHECK(cudaEventDestroy(input_ready_event));

  // wait for all the streams to finish
  CHECK(cudaDeviceSynchronize());
}

/*******************************************************************************
 * \brief Integrate all tasks of in given list onto given grids.
 *        See grid_task_list.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_gpu_integrate_task_list(
    const grid_gpu_task_list *task_list, const bool orthorhombic,
    const bool compute_tau, const int natoms, const int nlevels,
    const int npts_global[][3], const int npts_local[][3],
    const int shift_local[][3], const int border_width[][3],
    const double dh[][3][3], const double dh_inv[][3][3],
    const grid_buffer *pab_blocks, const double *grid[],
    grid_buffer *hab_blocks, double forces[][3], double virial[3][3]) {

  assert(task_list->nlevels == nlevels);

  // Prepare shared buffers using the main stream
  double *forces_dev = NULL;
  double *virial_dev = NULL;
  double *pab_blocks_dev = NULL;
  const size_t forces_size = 3 * natoms * sizeof(double);
  const size_t virial_size = 9 * sizeof(double);
  if (forces != NULL || virial != NULL) {
    CHECK(cudaMemcpyAsync(pab_blocks->device_buffer, pab_blocks->host_buffer,
                          pab_blocks->size, cudaMemcpyHostToDevice,
                          task_list->main_stream));
    pab_blocks_dev = pab_blocks->device_buffer;
  }
  if (forces != NULL) {
    CHECK(cudaMalloc(&forces_dev, forces_size));
    CHECK(cudaMemsetAsync(forces_dev, 0, forces_size, task_list->main_stream));
  }
  if (virial != NULL) {
    CHECK(cudaMalloc(&virial_dev, virial_size));
    CHECK(cudaMemsetAsync(virial_dev, 0, virial_size, task_list->main_stream));
  }

  // zero device hab blocks buffers
  CHECK(cudaMemsetAsync(hab_blocks->device_buffer, 0, hab_blocks->size,
                        task_list->main_stream));

  // record event so other streams can wait for hab, pab, virial etc to be ready
  cudaEvent_t input_ready_event;
  CHECK(cudaEventCreate(&input_ready_event));
  CHECK(cudaEventRecord(input_ready_event, task_list->main_stream));

  int first_task = 0;
  for (int level = 0; level < task_list->nlevels; level++) {
    const int last_task = first_task + task_list->tasks_per_level[level] - 1;
    const cudaStream_t level_stream = task_list->level_streams[level];
    const size_t grid_size = npts_local[level][0] * npts_local[level][1] *
                             npts_local[level][2] * sizeof(double);

    // reallocate device grid buffer if needed
    if (task_list->grid_dev_size[level] < grid_size) {
      if (task_list->grid_dev[level] != NULL) {
        CHECK(cudaFree(task_list->grid_dev[level]));
      }
      CHECK(cudaMalloc(&task_list->grid_dev[level], grid_size));
      task_list->grid_dev_size[level] = grid_size;
    }

    // upload grid
    CHECK(cudaMemcpyAsync(task_list->grid_dev[level], grid[level], grid_size,
                          cudaMemcpyHostToDevice, level_stream));

    // launch kernel, but only after grid has arrived
    CHECK(cudaStreamWaitEvent(level_stream, input_ready_event, 0));
    grid_gpu_integrate_one_grid_level(
        task_list, first_task, last_task, orthorhombic, compute_tau,
        npts_global[level], npts_local[level], shift_local[level],
        border_width[level], dh[level], dh_inv[level], level_stream,
        pab_blocks_dev, task_list->grid_dev[level], hab_blocks->device_buffer,
        forces_dev, virial_dev);

    // Have main stream wait for level to complete before downloading results.
    cudaEvent_t level_done_event;
    CHECK(cudaEventCreate(&level_done_event));
    CHECK(cudaEventRecord(level_done_event, level_stream));
    CHECK(cudaStreamWaitEvent(task_list->main_stream, level_done_event, 0));
    CHECK(cudaEventDestroy(level_done_event));

    first_task = last_task + 1;
  }

  // download result from device to host using main stream.
  CHECK(cudaMemcpyAsync(hab_blocks->host_buffer, hab_blocks->device_buffer,
                        hab_blocks->size, cudaMemcpyDeviceToHost,
                        task_list->main_stream));
  if (forces != NULL) {
    CHECK(cudaMemcpyAsync(forces, forces_dev, forces_size,
                          cudaMemcpyDeviceToHost, task_list->main_stream));
  }
  if (virial != NULL) {
    CHECK(cudaMemcpyAsync(virial, virial_dev, virial_size,
                          cudaMemcpyDeviceToHost, task_list->main_stream));
  }

  // wait for all the streams to finish
  CHECK(cudaDeviceSynchronize());

  // clean up
  CHECK(cudaEventDestroy(input_ready_event));
  if (forces != NULL) {
    CHECK(cudaFree(forces_dev));
  }
  if (virial != NULL) {
    CHECK(cudaFree(virial_dev));
  }
}

#endif // __GRID_CUDA
// EOF

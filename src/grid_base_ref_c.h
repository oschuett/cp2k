/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_BASE_REF_C_H
#define GRID_BASE_REF_C_H

#ifdef __cplusplus
extern "C" {
#endif

//******************************************************************************
// \brief ...
//
// \param grid_size_{x,y,z} ...
// \param grid_lbound_{x,y,z} ...
// \param compute_tau ...
// \param use_subpatch ...
// \param la_max ...
// \param la_min ...
// \param lb_max ...
// \param lb_min ...
// \param zeta ...
// \param zetb ...
// \param rscale ...
// \param rab2 ...
// \param dh ...
// \param dh_inv ...
// \param ra ...
// \param rab ...
// \param ng ...
// \param lb_grid ...
// \param perd ...
// \param lmax ...  // lmax_global  ! only for general
// \param radius ...  ! only for general
// \param lb_cube ...  ! only for ortho
// \param ub_cube ...  ! only for ortho
// \param sphere_bounds ...  ! only for ortho
// \param maxco ...
// \param o1 ...
// \param o2 ...
// \param pab ...
// \param grid ...
//******************************************************************************
void grid_collocate_pgf_product_rspace(const int grid_size_x,
                                       const int grid_size_y,
                                       const int grid_size_z,
                                       const int grid_lbound_x,
                                       const int grid_lbound_y,
                                       const int grid_lbound_z,
                                       const bool compute_tau,
                                       const bool use_subpatch,
                                       const int la_max,
                                       const int la_min,
                                       const int lb_max,
                                       const int lb_min,
                                       const double zeta,
                                       const double zetb,
                                       const double rscale,
                                       const double rab2,
                                       const double dh[3][3],
                                       const double dh_inv[3][3],
                                       const double ra[3],
                                       const double rab[3],
                                       const int ng[3],
                                       const int lb_grid[3],
                                       const int perd[3],
                                       const int lmax,
                                       const double radius,
                                       const int lb_cube[3],
                                       const int ub_cube[3],
                                       const int *sphere_bounds,
                                       const int maxco,
                                       const int o1,
                                       const int o2,
                                       const double pab[maxco][maxco],
                                       double grid[grid_size_z][grid_size_y][grid_size_x]);

#ifdef __cplusplus
}
#endif

#endif
//EOF

/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../common/grid_common.h"
#include "../common/grid_library.h"

#if (GRID_DO_COLLOCATE)
#define GRID_CONST_WHEN_COLLOCATE const
#define GRID_CONST_WHEN_INTEGRATE
#else
#define GRID_CONST_WHEN_COLLOCATE
#define GRID_CONST_WHEN_INTEGRATE const
#endif

/*******************************************************************************
 * \brief Collocates coefficients C_s onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cs_to_grid(const int k, const int k2, const int j,
                                    const int j2, const int i, const int i2,
                                    const int npts_local[3],
                                    GRID_CONST_WHEN_COLLOCATE double *cs,
                                    GRID_CONST_WHEN_INTEGRATE double *grid) {

  const int stride = npts_local[1] * npts_local[0];
  const int grid_index_0 = k * stride + j * npts_local[0] + i;
  const int grid_index_1 = k2 * stride + j * npts_local[0] + i;
  const int grid_index_2 = k * stride + j2 * npts_local[0] + i;
  const int grid_index_3 = k2 * stride + j2 * npts_local[0] + i;
  const int grid_index_4 = k * stride + j * npts_local[0] + i2;
  const int grid_index_5 = k2 * stride + j * npts_local[0] + i2;
  const int grid_index_6 = k * stride + j2 * npts_local[0] + i2;
  const int grid_index_7 = k2 * stride + j2 * npts_local[0] + i2;

#if (GRID_DO_COLLOCATE)
  // collocate
  grid[grid_index_0] += cs[0];
  grid[grid_index_1] += cs[1];
  grid[grid_index_2] += cs[2];
  grid[grid_index_3] += cs[3];
  grid[grid_index_4] += cs[4];
  grid[grid_index_5] += cs[5];
  grid[grid_index_6] += cs[6];
  grid[grid_index_7] += cs[7];
#else
  // integrate
  cs[0] += grid[grid_index_0];
  cs[1] += grid[grid_index_1];
  cs[2] += grid[grid_index_2];
  cs[3] += grid[grid_index_3];
  cs[4] += grid[grid_index_4];
  cs[5] += grid[grid_index_5];
  cs[6] += grid[grid_index_6];
  cs[7] += grid[grid_index_7];
#endif
}

/*******************************************************************************
 * \brief Transforms coefficients C_x into C_s by fixing grid index i.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cx_to_cs(const int lp, const double pol_ig[lp + 1],
                                  const double pol_ig2[lp + 1],
                                  GRID_CONST_WHEN_COLLOCATE double *cx,
                                  GRID_CONST_WHEN_INTEGRATE double *cs) {

  for (int lxp = 0; lxp <= lp; lxp++) {
    const double p1 = pol_ig[lxp];
    const double p2 = pol_ig2[lxp];

#if (GRID_DO_COLLOCATE)
    // collocate
    cs[0] += cx[lxp * 4 + 0] * p1;
    cs[1] += cx[lxp * 4 + 1] * p1;
    cs[2] += cx[lxp * 4 + 2] * p1;
    cs[3] += cx[lxp * 4 + 3] * p1;
    cs[4] += cx[lxp * 4 + 0] * p2;
    cs[5] += cx[lxp * 4 + 1] * p2;
    cs[6] += cx[lxp * 4 + 2] * p2;
    cs[7] += cx[lxp * 4 + 3] * p2;
#else
    // integrate
    cx[lxp * 4 + 0] += cs[0] * p1;
    cx[lxp * 4 + 1] += cs[1] * p1;
    cx[lxp * 4 + 2] += cs[2] * p1;
    cx[lxp * 4 + 3] += cs[3] * p1;
    cx[lxp * 4 + 0] += cs[4] * p2;
    cx[lxp * 4 + 1] += cs[5] * p2;
    cx[lxp * 4 + 2] += cs[6] * p2;
    cx[lxp * 4 + 3] += cs[7] * p2;
#endif
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_x onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
ortho_cx_to_grid(const int lp, const int k, const int k2, const int jg,
                 const int jg2, const int cmax, const double kremain,
                 const double pol[3][2 * cmax + 1][lp + 1],
                 const int map[3][2 * cmax + 1], const double dh[3][3],
                 const double dh_inv[3][3], const int npts_local[3],
                 GRID_CONST_WHEN_COLLOCATE double *cx,
                 GRID_CONST_WHEN_INTEGRATE double *grid) {

  const int j = map[1][jg + cmax];
  const int j2 = map[1][jg2 + cmax];
  const int jd = (2 * jg - 1) / 2; // distance from center in grid points
  const double jr = jd * dh[1][1]; // distance from center in a.u.
  const double jremain = kremain - jr * jr;
  const int igmin = ceil(-1e-8 - sqrt(fmax(0.0, jremain)) * dh_inv[0][0]);
  for (int ig = igmin; ig <= 0; ig++) {
    const int ig2 = 1 - ig;
    const int i = map[0][ig + cmax];
    const int i2 = map[0][ig2 + cmax];
    double cs[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#if (GRID_DO_COLLOCATE)
    // collocate
    ortho_cx_to_cs(lp, pol[0][ig + cmax], pol[0][ig2 + cmax], cx, cs);
    ortho_cs_to_grid(k, k2, j, j2, i, i2, npts_local, cs, grid);
#else
    // integrate
    ortho_cs_to_grid(k, k2, j, j2, i, i2, npts_local, cs, grid);
    ortho_cx_to_cs(lp, pol[0][ig + cmax], pol[0][ig2 + cmax], cx, cs);
#endif
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_xy into C_x by fixing grid index j.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cxy_to_cx(const int lp, const double pol_jg[lp + 1],
                                   const double pol_jg2[lp + 1],
                                   GRID_CONST_WHEN_COLLOCATE double *cxy,
                                   GRID_CONST_WHEN_INTEGRATE double *cx) {

  for (int lyp = 0; lyp <= lp; lyp++) {
    for (int lxp = 0; lxp <= lp - lyp; lxp++) {
      const int cxy_index = lyp * (lp + 1) * 2 + lxp * 2; // [lyp, lxp, 0]

#if (GRID_DO_COLLOCATE)
      // collocate
      cx[lxp * 4 + 0] += cxy[cxy_index + 0] * pol_jg[lyp];
      cx[lxp * 4 + 1] += cxy[cxy_index + 1] * pol_jg[lyp];
      cx[lxp * 4 + 2] += cxy[cxy_index + 0] * pol_jg2[lyp];
      cx[lxp * 4 + 3] += cxy[cxy_index + 1] * pol_jg2[lyp];
#else
      // integrate
      cxy[cxy_index + 0] += cx[lxp * 4 + 0] * pol_jg[lyp];
      cxy[cxy_index + 1] += cx[lxp * 4 + 1] * pol_jg[lyp];
      cxy[cxy_index + 0] += cx[lxp * 4 + 2] * pol_jg2[lyp];
      cxy[cxy_index + 1] += cx[lxp * 4 + 3] * pol_jg2[lyp];
#endif
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_xy onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cxy_to_grid(
    const int lp, const int kg, const int kg2, const int cmax,
    const double pol[3][2 * cmax + 1][lp + 1], const int map[3][2 * cmax + 1],
    const double dh[3][3], const double dh_inv[3][3], const double disr_radius,
    const int npts_local[3], GRID_CONST_WHEN_COLLOCATE double *cxy,
    GRID_CONST_WHEN_INTEGRATE double *grid) {

  // The cube contains an even number of grid points in each direction and
  // collocation is always performed on a pair of two opposing grid points.
  // Hence, the points with index 0 and 1 are both assigned distance zero via
  // the formular distance=(2*index-1)/2.

  const int k = map[2][kg + cmax];
  const int k2 = map[2][kg2 + cmax];
  const int kd = (2 * kg - 1) / 2; // distance from center in grid points
  const double kr = kd * dh[2][2]; // distance from center in a.u.
  const double kremain = disr_radius * disr_radius - kr * kr;
  const int jgmin = ceil(-1e-8 - sqrt(fmax(0.0, kremain)) * dh_inv[1][1]);
  for (int jg = jgmin; jg <= 0; jg++) {
    const int jg2 = 1 - jg;

    const size_t cx_size = (lp + 1) * 4;
    double cx[cx_size];
    memset(cx, 0, cx_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
    // collocate
    ortho_cxy_to_cx(lp, pol[1][jg + cmax], pol[1][jg2 + cmax], cxy, cx);
    ortho_cx_to_grid(lp, k, k2, jg, jg2, cmax, kremain, pol, map, dh, dh_inv,
                     npts_local, cx, grid);
#else
    // integrate
    ortho_cx_to_grid(lp, k, k2, jg, jg2, cmax, kremain, pol, map, dh, dh_inv,
                     npts_local, cx, grid);
    ortho_cxy_to_cx(lp, pol[1][jg + cmax], pol[1][jg2 + cmax], cxy, cx);
#endif
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_xyz into C_xz by fixing grid index k.
 * \author Ole Schuett
 ******************************************************************************/
static inline void ortho_cxyz_to_cxy(const int lp, const double pol_kg[lp + 1],
                                     const double pol_kg2[lp + 1],
                                     GRID_CONST_WHEN_COLLOCATE double *cxyz,
                                     GRID_CONST_WHEN_INTEGRATE double *cxy) {

  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++) {
        const int cxyz_index =
            lzp * (lp + 1) * (lp + 1) + lyp * (lp + 1) + lxp; // [lzp, lyp, lxp]
        const int cxy_index = lyp * (lp + 1) * 2 + lxp * 2;   // [lyp, lxp, 0]

#if (GRID_DO_COLLOCATE)
        // collocate
        cxy[cxy_index + 0] += cxyz[cxyz_index] * pol_kg[lzp];
        cxy[cxy_index + 1] += cxyz[cxyz_index] * pol_kg2[lzp];
#else
        // integrate
        cxyz[cxyz_index] += cxy[cxy_index + 0] * pol_kg[lzp];
        cxyz[cxyz_index] += cxy[cxy_index + 1] * pol_kg2[lzp];
#endif
      }
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for orthorhombic case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
ortho_cxyz_to_grid(const int lp, const double zetp, const double dh[3][3],
                   const double dh_inv[3][3], const double rp[3],
                   const int npts_global[3], const int npts_local[3],
                   const int shift_local[3], const double radius,
                   GRID_CONST_WHEN_COLLOCATE double *cxyz,
                   GRID_CONST_WHEN_INTEGRATE double *grid) {

  // *** position of the gaussian product
  //
  // this is the actual definition of the position on the grid
  // i.e. a point rp(:) gets here grid coordinates
  // MODULO(rp(:)/dr(:),npts_global(:))+1
  // hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid in Fortran
  // and (1,1,1) on grid here in C.

  // cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
  int cubecenter[3];
  for (int i = 0; i < 3; i++) {
    double dh_inv_rp = 0.0;
    for (int j = 0; j < 3; j++) {
      dh_inv_rp += dh_inv[j][i] * rp[j];
    }
    cubecenter[i] = floor(dh_inv_rp);
  }

  double roffset[3];
  for (int i = 0; i < 3; i++) {
    roffset[i] = rp[i] - ((double)cubecenter[i]) * dh[i][i];
  }

  // Historically, the radius gets discretized.
  const double drmin = fmin(dh[0][0], fmin(dh[1][1], dh[2][2]));
  const double disr_radius = drmin * fmax(1.0, ceil(radius / drmin));

  int lb_cube[3], ub_cube[3];
  for (int i = 0; i < 3; i++) {
    lb_cube[i] = ceil(-1e-8 - disr_radius * dh_inv[i][i]);
    ub_cube[i] = 1 - lb_cube[i];
    // If grid is not period check that cube fits without wrapping.
    if (npts_global[i] != npts_local[i]) {
      const int offset =
          modulo(cubecenter[i] + lb_cube[i] - shift_local[i], npts_global[i]) -
          lb_cube[i];
      assert(offset + ub_cube[i] < npts_local[i]);
      assert(offset + lb_cube[i] >= 0);
    }
  }

  // cmax = MAXVAL(ub_cube)
  const int cmax = imax(imax(ub_cube[0], ub_cube[1]), ub_cube[2]);

  // Precompute (x-xp)**lp*exp(..) for each direction.
  double pol_mutable[3][2 * cmax + 1][lp + 1];
  for (int idir = 0; idir < 3; idir++) {
    const double dr = dh[idir][idir];
    const double ro = roffset[idir];
    //  Reuse the result from the previous gridpoint to avoid to many exps:
    //  exp( -a*(x+d)**2) = exp(-a*x**2)*exp(-2*a*x*d)*exp(-a*d**2)
    //  exp(-2*a*(x+d)*d) = exp(-2*a*x*d)*exp(-2*a*d**2)
    const double t_exp_1 = exp(-zetp * pow(dr, 2));
    const double t_exp_2 = pow(t_exp_1, 2);
    double t_exp_min_1 = exp(-zetp * pow(+dr - ro, 2));
    double t_exp_min_2 = exp(-2 * zetp * (+dr - ro) * (-dr));
    for (int ig = 0; ig >= lb_cube[idir]; ig--) {
      const double rpg = ig * dr - ro;
      t_exp_min_1 *= t_exp_min_2 * t_exp_1;
      t_exp_min_2 *= t_exp_2;
      double pg = t_exp_min_1;
      for (int icoef = 0; icoef <= lp; icoef++) {
        pol_mutable[idir][ig + cmax][icoef] = pg; // exp(-zetp*rpg**2)
        pg *= rpg;
      }
    }
    double t_exp_plus_1 = exp(-zetp * pow(-ro, 2));
    double t_exp_plus_2 = exp(-2 * zetp * (-ro) * (+dr));
    for (int ig = 0; ig >= lb_cube[idir]; ig--) {
      const double rpg = (1 - ig) * dr - ro;
      t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
      t_exp_plus_2 *= t_exp_2;
      double pg = t_exp_plus_1;
      for (int icoef = 0; icoef <= lp; icoef++) {
        pol_mutable[idir][1 - ig + cmax][icoef] = pg; // exp(-zetp*rpg**2)
        pg *= rpg;
      }
    }
  }
  const double(*pol)[2 * cmax + 1][lp + 1] =
      (const double(*)[2 * cmax + 1][lp + 1]) pol_mutable;

  // Precompute mapping from cube to grid indices for each direction
  int map_mutable[3][2 * cmax + 1];
  for (int i = 0; i < 3; i++) {
    for (int k = -cmax; k <= +cmax; k++) {
      map_mutable[i][k + cmax] =
          modulo(cubecenter[i] + k - shift_local[i], npts_global[i]);
    }
  }
  const int(*map)[2 * cmax + 1] = (const int(*)[2 * cmax + 1]) map_mutable;

  // Loop over k dimension of the cube.
  const int kgmin = ceil(-1e-8 - disr_radius * dh_inv[2][2]);
  for (int kg = kgmin; kg <= 0; kg++) {
    const int kg2 = 1 - kg;

    const size_t cxy_size = (lp + 1) * (lp + 1) * 2;
    double cxy[cxy_size];
    memset(cxy, 0, cxy_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
    // collocate
    grid_library_gather_stats((grid_library_stats){.ref_collocate_ortho = 1});
    ortho_cxyz_to_cxy(lp, pol[2][kg + cmax], pol[2][kg2 + cmax], cxyz, cxy);
    ortho_cxy_to_grid(lp, kg, kg2, cmax, pol, map, dh, dh_inv, disr_radius,
                      npts_local, cxy, grid);
#else
    // integrate
    grid_library_gather_stats((grid_library_stats){.ref_integrate_ortho = 1});
    ortho_cxy_to_grid(lp, kg, kg2, cmax, pol, map, dh, dh_inv, disr_radius,
                      npts_local, cxy, grid);
    ortho_cxyz_to_cxy(lp, pol[2][kg + cmax], pol[2][kg2 + cmax], cxyz, cxy);
#endif
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_i onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_ci_to_grid(const int lp, const int j, const int jg, const int k,
                   const int kg, const int npts_local[3], const int bounds_i[2],
                   const int index_min[3], const int map_i[], const double zetp,
                   const double dh[3][3], const double gp[3],
                   const double radius, GRID_CONST_WHEN_COLLOCATE double *ci,
                   GRID_CONST_WHEN_INTEGRATE double *grid) {

  //--------------------------------------------------------------------
  // Find bounds for the inner loop based on a quadratic equation in i.
  //
  // The real-space vector from the center of the gaussian to the
  // grid point i,j,k is given by:
  //   r = (i-gp[0])*dh[0,:] + (j-gp[1])*dh[1,:] + (k-gp[2])*dh[2,:]
  //
  // Separating the term that depends on i:
  //   r = i*dh[0,:] - gp[0]*dh[0,:] + (j-gp[1])*dh[1,:] + (k-gp[2])*dh[2,:]
  //     = i*dh[0,:] + v
  //
  // The squared distance works out to:
  //   r**2 = dh[0,:]**2 * i**2  +  2 * v * dh[0,:] * i  +  v**2
  //        = a * i**2           +  b * i                +  c
  //
  // Solving r**2==radius**2 for i yields:
  //    d =  b**2  -  4 * a * (c - radius**2)
  //    i = (-b \pm sqrt(d)) / (2*a)
  //
  double a = 0.0, b = 0.0, c = 0.0;
  for (int i = 0; i < 3; i++) {
    const double v = (0 - gp[0]) * dh[0][i] + (j - gp[1]) * dh[1][i] +
                     (k - gp[2]) * dh[2][i];
    a += dh[0][i] * dh[0][i];
    b += 2.0 * v * dh[0][i];
    c += v * v;
  }
  const double d = b * b - 4.0 * a * (c - radius * radius);
  if (d <= 0.0) {
    return;
  }
  const double sqrt_d = sqrt(d);
  const int ismin = ceil((-b - sqrt_d) / (2.0 * a));
  const int ismax = floor((-b + sqrt_d) / (2.0 * a));

  const double exp_ab = exp(-zetp * (a + b));
  const double exp_2a = exp(-zetp * 2.0 * a);
  bool exp_is_invalid = true;
  double exp_2ai = 0.0;
  double exp_aiibic = 0.0;

  for (int i = ismin; i <= ismax; i++) {
    const int ig = map_i[i - index_min[0]];
    if (ig < bounds_i[0] || bounds_i[1] < ig) {
      exp_is_invalid = true;
      continue;
    }

    const int stride = npts_local[1] * npts_local[0];
    const int grid_index =
        kg * stride + jg * npts_local[0] + ig; // [kg, jg, ig]

    if (exp_is_invalid) {
      exp_2ai = exp(-zetp * 2.0 * a * i);
      exp_aiibic = exp(-zetp * ((a * i + b) * i + c));
      exp_is_invalid = false;
    }

    // polynomial terms
    double dip = exp_aiibic; // exp(-zetp * ((a * i + b) * i + c));
    const double di = i - gp[0];
    for (int il = 0; il <= lp; il++) {
#if (GRID_DO_COLLOCATE)
      grid[grid_index] += ci[il] * dip; // collocate
#else
      ci[il] += grid[grid_index] * dip; // integrate
#endif
      dip *= di;
    }

    // update exponential term
    exp_aiibic *= exp_2ai * exp_ab;
    exp_2ai *= exp_2a;
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_ij into C_i by fixing grid index j.
 * \author Ole Schuett
 ******************************************************************************/
static inline void general_cij_to_ci(const int lp, const double dj,
                                     GRID_CONST_WHEN_COLLOCATE double *cij,
                                     GRID_CONST_WHEN_INTEGRATE double *ci) {
  double djp = 1.0;
  for (int jl = 0; jl <= lp; jl++) {
    for (int il = 0; il <= lp - jl; il++) {
      const int cij_index = jl * (lp + 1) + il; // [jl, il]
#if (GRID_DO_COLLOCATE)
      ci[il] += cij[cij_index] * djp; // collocate
#else
      cij[cij_index] += ci[il] * djp;   // integrate
#endif
    }
    djp *= dj;
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_ij onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void general_cij_to_grid(
    const int lp, const int k, const int kg, const int npts_local[3],
    const int bounds_i[2], const int bounds_j[2], const int index_min[3],
    const int index_max[3], const int map_i[], const int map_j[],
    const double zetp, const double dh[3][3], const double gp[3],
    const double radius, GRID_CONST_WHEN_COLLOCATE double *cij,
    GRID_CONST_WHEN_INTEGRATE double *grid) {

  for (int j = index_min[1]; j <= index_max[1]; j++) {
    const int jg = map_j[j - index_min[1]];
    if (jg < bounds_j[0] || bounds_j[1] < jg) {
      continue;
    }
    const double dj = j - gp[1];

    const size_t ci_size = lp + 1;
    double ci[ci_size];
    memset(ci, 0, ci_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
    // collocate
    general_cij_to_ci(lp, dj, cij, ci);
    general_ci_to_grid(lp, j, jg, k, kg, npts_local, bounds_i, index_min, map_i,
                       zetp, dh, gp, radius, ci, grid);
#else
    // integrate
    general_ci_to_grid(lp, j, jg, k, kg, npts_local, bounds_i, index_min, map_i,
                       zetp, dh, gp, radius, ci, grid);
    general_cij_to_ci(lp, dj, cij, ci);
#endif
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_ijk into C_ij by fixing grid index k.
 * \author Ole Schuett
 ******************************************************************************/
static inline void general_cijk_to_cij(const int lp, const double dk,
                                       GRID_CONST_WHEN_COLLOCATE double *cijk,
                                       GRID_CONST_WHEN_INTEGRATE double *cij) {
  double dkp = 1.0;
  for (int kl = 0; kl <= lp; kl++) {
    for (int jl = 0; jl <= lp - kl; jl++) {
      for (int il = 0; il <= lp - kl - jl; il++) {
        const int cij_index = jl * (lp + 1) + il; // [jl, il]
        const int cijk_index =
            kl * (lp + 1) * (lp + 1) + jl * (lp + 1) + il; // [kl, jl, il]
#if (GRID_DO_COLLOCATE)
        cij[cij_index] += cijk[cijk_index] * dkp; // collocate
#else
        cijk[cijk_index] += cij[cij_index] * dkp; // integrate
#endif
      }
    }
    dkp *= dk;
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_ijk onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_cijk_to_grid(const int border_mask, const int lp, const double zetp,
                     const double dh[3][3], const double dh_inv[3][3],
                     const double rp[3], const int npts_global[3],
                     const int npts_local[3], const int shift_local[3],
                     const int border_width[3], const double radius,
                     GRID_CONST_WHEN_COLLOCATE double *cijk,
                     GRID_CONST_WHEN_INTEGRATE double *grid) {

  // Default for border_mask == 0.
  int bounds_i[2] = {0, npts_local[0] - 1};
  int bounds_j[2] = {0, npts_local[1] - 1};
  int bounds_k[2] = {0, npts_local[2] - 1};

  // See also rs_find_node() in task_list_methods.F.
  // If the bit is set then we need to exclude the border in that direction.
  if (border_mask & (1 << 0))
    bounds_i[0] += border_width[0];
  if (border_mask & (1 << 1))
    bounds_i[1] -= border_width[0];
  if (border_mask & (1 << 2))
    bounds_j[0] += border_width[1];
  if (border_mask & (1 << 3))
    bounds_j[1] -= border_width[1];
  if (border_mask & (1 << 4))
    bounds_k[0] += border_width[2];
  if (border_mask & (1 << 5))
    bounds_k[1] -= border_width[2];

  // center in grid coords
  // gp = MATMUL(dh_inv, rp)
  double gp[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gp[i] += dh_inv[j][i] * rp[j];
    }
  }

  // CALL return_cube_nonortho(cube_info, radius, index_min, index_max, rp)
  //
  // get the min max indices that contain at least the cube that contains a
  // sphere around rp of radius radius if the cell is very non-orthogonal this
  // implies that many useless points are included this estimate can be improved
  // (i.e. not box but sphere should be used)
  int index_min[3] = {INT_MAX, INT_MAX, INT_MAX};
  int index_max[3] = {INT_MIN, INT_MIN, INT_MIN};
  for (int i = -1; i <= 1; i++) {
    for (int j = -1; j <= 1; j++) {
      for (int k = -1; k <= 1; k++) {
        const double x = rp[0] + i * radius;
        const double y = rp[1] + j * radius;
        const double z = rp[2] + k * radius;
        for (int idir = 0; idir < 3; idir++) {
          const double resc =
              dh_inv[0][idir] * x + dh_inv[1][idir] * y + dh_inv[2][idir] * z;
          index_min[idir] = imin(index_min[idir], floor(resc));
          index_max[idir] = imax(index_max[idir], ceil(resc));
        }
      }
    }
  }

  // precompute modulos
  int map_k[index_max[2] - index_min[2] + 1];
  for (int k = index_min[2]; k <= index_max[2]; k++) {
    map_k[k - index_min[2]] = modulo(k - shift_local[2], npts_global[2]);
  }
  int map_j[index_max[1] - index_min[1] + 1];
  for (int j = index_min[1]; j <= index_max[1]; j++) {
    map_j[j - index_min[1]] = modulo(j - shift_local[1], npts_global[1]);
  }
  int map_i[index_max[0] - index_min[0] + 1];
  for (int i = index_min[0]; i <= index_max[0]; i++) {
    map_i[i - index_min[0]] = modulo(i - shift_local[0], npts_global[0]);
  }

  // go over the grid, but cycle if the point is not within the radius
  for (int k = index_min[2]; k <= index_max[2]; k++) {
    const int kg = map_k[k - index_min[2]];
    if (kg < bounds_k[0] || bounds_k[1] < kg) {
      continue;
    }
    const double dk = k - gp[2];

    // zero coef_xyt
    const int cij_size = (lp + 1) * (lp + 1);
    double cij[cij_size];
    memset(cij, 0, cij_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
    // collocate
    general_cijk_to_cij(lp, dk, cijk, cij);
    general_cij_to_grid(lp, k, kg, npts_local, bounds_i, bounds_j, index_min,
                        index_max, map_i, map_j, zetp, dh, gp, radius, cij,
                        grid);
#else
    // integrate
    general_cij_to_grid(lp, k, kg, npts_local, bounds_i, bounds_j, index_min,
                        index_max, map_i, map_j, zetp, dh, gp, radius, cij,
                        grid);
    general_cijk_to_cij(lp, dk, cijk, cij);
#endif
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_xyz into C_ijk.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_cxyz_to_cijk(const int lp, const double dh[3][3],
                     GRID_CONST_WHEN_COLLOCATE double *cxyz,
                     GRID_CONST_WHEN_INTEGRATE double *cijk) {

  // transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
  // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
  // sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp

  // transform using multinomials
  double hmatgridp[lp + 1][3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      hmatgridp[0][j][i] = 1.0;
      for (int k = 1; k <= lp; k++) {
        hmatgridp[k][j][i] = hmatgridp[k - 1][j][i] * dh[j][i];
      }
    }
  }

  const int lpx = lp;
  for (int klx = 0; klx <= lpx; klx++) {
    for (int jlx = 0; jlx <= lpx - klx; jlx++) {
      for (int ilx = 0; ilx <= lpx - klx - jlx; ilx++) {
        const int lx = ilx + jlx + klx;
        const int lpy = lp - lx;
        for (int kly = 0; kly <= lpy; kly++) {
          for (int jly = 0; jly <= lpy - kly; jly++) {
            for (int ily = 0; ily <= lpy - kly - jly; ily++) {
              const int ly = ily + jly + kly;
              const int lpz = lp - lx - ly;
              for (int klz = 0; klz <= lpz; klz++) {
                for (int jlz = 0; jlz <= lpz - klz; jlz++) {
                  for (int ilz = 0; ilz <= lpz - klz - jlz; ilz++) {
                    const int lz = ilz + jlz + klz;
                    const int il = ilx + ily + ilz;
                    const int jl = jlx + jly + jlz;
                    const int kl = klx + kly + klz;
                    const int lp1 = lp + 1;
                    const int cijk_index =
                        kl * lp1 * lp1 + jl * lp1 + il; // [kl,jl,il]
                    const int cxyz_index =
                        lz * lp1 * lp1 + ly * lp1 + lx; // [lz,ly,lx]
                    const double p =
                        hmatgridp[ilx][0][0] * hmatgridp[jlx][1][0] *
                        hmatgridp[klx][2][0] * hmatgridp[ily][0][1] *
                        hmatgridp[jly][1][1] * hmatgridp[kly][2][1] *
                        hmatgridp[ilz][0][2] * hmatgridp[jlz][1][2] *
                        hmatgridp[klz][2][2] * fac[lx] * fac[ly] * fac[lz] /
                        (fac[ilx] * fac[ily] * fac[ilz] * fac[jlx] * fac[jly] *
                         fac[jlz] * fac[klx] * fac[kly] * fac[klz]);
#if (GRID_DO_COLLOCATE)
                    cijk[cijk_index] += cxyz[cxyz_index] * p; // collocate
#else
                    cxyz[cxyz_index] += cijk[cijk_index] * p; // integrate
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid for general case.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
general_cxyz_to_grid(const int border_mask, const int lp, const double zetp,
                     const double dh[3][3], const double dh_inv[3][3],
                     const double rp[3], const int npts_global[3],
                     const int npts_local[3], const int shift_local[3],
                     const int border_width[3], const double radius,
                     GRID_CONST_WHEN_COLLOCATE double *cxyz,
                     GRID_CONST_WHEN_INTEGRATE double *grid) {

  const size_t cijk_size = (lp + 1) * (lp + 1) * (lp + 1);
  double cijk[cijk_size];
  memset(cijk, 0, cijk_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
  // collocate
  grid_library_gather_stats((grid_library_stats){.ref_collocate_general = 1});
  general_cxyz_to_cijk(lp, dh, cxyz, cijk);
  general_cijk_to_grid(border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
                       npts_local, shift_local, border_width, radius, cijk,
                       grid);
#else
  // integrate
  grid_library_gather_stats((grid_library_stats){.ref_integrate_general = 1});
  general_cijk_to_grid(border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
                       npts_local, shift_local, border_width, radius, cijk,
                       grid);
  general_cxyz_to_cijk(lp, dh, cxyz, cijk);
#endif
}

/*******************************************************************************
 * \brief Collocates coefficients C_xyz onto the grid.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
cxyz_to_grid(const bool orthorhombic, const int border_mask, const int lp,
             const double zetp, const double dh[3][3],
             const double dh_inv[3][3], const double rp[3],
             const int npts_global[3], const int npts_local[3],
             const int shift_local[3], const int border_width[3],
             const double radius, GRID_CONST_WHEN_COLLOCATE double *cxyz,
             GRID_CONST_WHEN_INTEGRATE double *grid) {

  if (orthorhombic && border_mask == 0) {
    ortho_cxyz_to_grid(lp, zetp, dh, dh_inv, rp, npts_global, npts_local,
                       shift_local, radius, cxyz, grid);
  } else {
    general_cxyz_to_grid(border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
                         npts_local, shift_local, border_width, radius, cxyz,
                         grid);
  }
}

/*******************************************************************************
 * \brief Transforms coefficients C_ab into C_xyz.
 * \author Ole Schuett
 ******************************************************************************/
static inline void cab_to_cxyz(const int la_max, const int la_min,
                               const int lb_max, const int lb_min,
                               const double prefactor, const double ra[3],
                               const double rb[3], const double rp[3],
                               GRID_CONST_WHEN_COLLOCATE double *cab,
                               GRID_CONST_WHEN_INTEGRATE double *cxyz) {

  // Computes the polynomial expansion coefficients:
  //     (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
  const int lp = la_max + lb_max;
  double alpha[3][lb_max + 1][la_max + 1][lp + 1];
  memset(alpha, 0, 3 * (lb_max + 1) * (la_max + 1) * (lp + 1) * sizeof(double));

  for (int i = 0; i < 3; i++) {
    const double drpa = rp[i] - ra[i];
    const double drpb = rp[i] - rb[i];
    for (int lxa = 0; lxa <= la_max; lxa++) {
      for (int lxb = 0; lxb <= lb_max; lxb++) {
        double binomial_k_lxa = 1.0;
        double a = 1.0;
        for (int k = 0; k <= lxa; k++) {
          double binomial_l_lxb = 1.0;
          double b = 1.0;
          for (int l = 0; l <= lxb; l++) {
            alpha[i][lxb][lxa][lxa - l + lxb - k] +=
                binomial_k_lxa * binomial_l_lxb * a * b;
            binomial_l_lxb *= ((double)(lxb - l)) / ((double)(l + 1));
            b *= drpb;
          }
          binomial_k_lxa *= ((double)(lxa - k)) / ((double)(k + 1));
          a *= drpa;
        }
      }
    }
  }

  //   *** initialise the coefficient matrix, we transform the sum
  //
  // sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} *
  //         (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya
  //         (z-a_z)**lza
  //
  // into
  //
  // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
  //
  // where p is center of the product gaussian, and lp = la_max + lb_max
  // (current implementation is l**7)
  //

  for (int lzb = 0; lzb <= lb_max; lzb++) {
    for (int lza = 0; lza <= la_max; lza++) {
      for (int lyb = 0; lyb <= lb_max - lzb; lyb++) {
        for (int lya = 0; lya <= la_max - lza; lya++) {
          const int lxb_min = imax(lb_min - lzb - lyb, 0);
          const int lxa_min = imax(la_min - lza - lya, 0);
          for (int lxb = lxb_min; lxb <= lb_max - lzb - lyb; lxb++) {
            for (int lxa = lxa_min; lxa <= la_max - lza - lya; lxa++) {
              const int ico = coset(lxa, lya, lza);
              const int jco = coset(lxb, lyb, lzb);
              const int cab_index = jco * ncoset[la_max] + ico; // [jco, ico]
              for (int lzp = 0; lzp <= lza + lzb; lzp++) {
                for (int lyp = 0; lyp <= lp - lza - lzb; lyp++) {
                  for (int lxp = 0; lxp <= lp - lza - lzb - lyp; lxp++) {
                    const double p = alpha[0][lxb][lxa][lxp] *
                                     alpha[1][lyb][lya][lyp] *
                                     alpha[2][lzb][lza][lzp] * prefactor;
                    const int lp1 = lp + 1;
                    const int cxyz_index =
                        lzp * lp1 * lp1 + lyp * lp1 + lxp; // [lzp, lyp, lxp]
#if (GRID_DO_COLLOCATE)
                    cxyz[cxyz_index] += cab[cab_index] * p; // collocate
#else
                    cab[cab_index] += cxyz[cxyz_index] * p; // integrate
#endif
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/*******************************************************************************
 * \brief Collocates coefficients C_ab onto the grid.
 * \author Ole Schuett
 ******************************************************************************/
static inline void
cab_to_grid(const bool orthorhombic, const int border_mask, const int la_max,
            const int la_min, const int lb_max, const int lb_min,
            const double zeta, const double zetb, const double rscale,
            const double dh[3][3], const double dh_inv[3][3],
            const double ra[3], const double rab[3], const int npts_global[3],
            const int npts_local[3], const int shift_local[3],
            const int border_width[3], const double radius,
            GRID_CONST_WHEN_COLLOCATE double *cab,
            GRID_CONST_WHEN_INTEGRATE double *grid) {

  // Check if radius is too small to be mapped onto grid of given resolution.
  double dh_max = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dh_max = fmax(dh_max, fabs(dh[i][j]));
    }
  }
  if (2.0 * radius < dh_max) {
    return;
  }

  const double zetp = zeta + zetb;
  const double f = zetb / zetp;
  const double rab2 = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
  const double prefactor = rscale * exp(-zeta * f * rab2);
  double rp[3], rb[3];
  for (int i = 0; i < 3; i++) {
    rp[i] = ra[i] + f * rab[i];
    rb[i] = ra[i] + rab[i];
  }

  const int lp = la_max + lb_max;
  const size_t cxyz_size = (lp + 1) * (lp + 1) * (lp + 1);
  double cxyz[cxyz_size];
  memset(cxyz, 0, cxyz_size * sizeof(double));

#if (GRID_DO_COLLOCATE)
  // collocate
  cab_to_cxyz(la_max, la_min, lb_max, lb_min, prefactor, ra, rb, rp, cab, cxyz);
  cxyz_to_grid(orthorhombic, border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
               npts_local, shift_local, border_width, radius, cxyz, grid);
#else
  // integrate
  cxyz_to_grid(orthorhombic, border_mask, lp, zetp, dh, dh_inv, rp, npts_global,
               npts_local, shift_local, border_width, radius, cxyz, grid);
  cab_to_cxyz(la_max, la_min, lb_max, lb_min, prefactor, ra, rb, rp, cab, cxyz);
#endif
}

// EOF

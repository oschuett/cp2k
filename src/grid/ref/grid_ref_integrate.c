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

#define GRID_DO_COLLOCATE 0
#include "../common/grid_common.h"
#include "../common/grid_process_vab.h"
#include "grid_ref_collint.h"
#include "grid_ref_integrate.h"

/*******************************************************************************
 * \brief Integrates a single task. See grid_ref_integrate.h for details.
 * \author Ole Schuett
 ******************************************************************************/
void grid_ref_integrate_pgf_product(
    const bool orthorhombic, const bool compute_tau, const int border_mask,
    const int la_max, const int la_min, const int lb_max, const int lb_min,
    const double zeta, const double zetb, const double dh[3][3],
    const double dh_inv[3][3], const double ra[3], const double rab[3],
    const int npts_global[3], const int npts_local[3], const int shift_local[3],
    const int border_width[3], const double radius, const int o1, const int o2,
    const int n1, const int n2, const double *grid, double hab[n2][n1],
    const double pab[n2][n1], double forces[2][3], double virials[2][3][3],
    double hdab[n2][n1][3], double a_hdab[n2][n1][3][3]) {

  const bool calculate_forces =
      (forces != NULL || hdab != NULL || virials != NULL || a_hdab != NULL);
  const bool calculate_virial = (virials != NULL || a_hdab != NULL);
  const process_ldiffs ldiffs =
      process_get_ldiffs(calculate_forces, calculate_virial, compute_tau);

  int la_min_cab = imax(la_min + ldiffs.la_min_diff, 0);
  int lb_min_cab = imax(lb_min + ldiffs.lb_min_diff, 0);
  int la_max_cab = la_max + ldiffs.la_max_diff;
  int lb_max_cab = lb_max + ldiffs.lb_max_diff;
  const int n1_cab = ncoset(la_max_cab);
  const int n2_cab = ncoset(lb_max_cab);

  const size_t cab_size = n2_cab * n1_cab;
  double cab_mutable[cab_size];
  memset(cab_mutable, 0, cab_size * sizeof(double));

  const double rscale = 1.0; // TODO: remove rscale from cab_to_grid
  cab_to_grid(orthorhombic, border_mask, la_max_cab, la_min_cab, lb_max_cab,
              lb_min_cab, zeta, zetb, rscale, dh, dh_inv, ra, rab, npts_global,
              npts_local, shift_local, border_width, radius, cab_mutable, grid);

  //  cab contains all the information needed to find the elements of hab
  //  and optionally of derivatives of these elements
  const double(*cab)[n1_cab] = (const double(*)[n1_cab])cab_mutable;

  const double ftza = 2.0 * zeta;
  const double ftzb = 2.0 * zetb;

  for (int la = la_min; la <= la_max; la++) {
    for (int ax = 0; ax <= la; ax++) {
      for (int ay = 0; ay <= la - ax; ay++) {
        const int az = la - ax - ay;
        const orbital a = {{ax, ay, az}};
        for (int lb = lb_min; lb <= lb_max; lb++) {
          for (int bx = 0; bx <= lb; bx++) {
            for (int by = 0; by <= lb - bx; by++) {
              const int bz = lb - bx - by;
              const orbital b = {{bx, by, bz}};

              double *habval = &hab[o2 + idx(b)][o1 + idx(a)];
              const double pabval =
                  (pab == NULL) ? 0.0 : pab[o2 + idx(b)][o1 + idx(a)];

              // Fill hab, forces, and virials.
              if (compute_tau) {
                process_tau(a, b, ftza, ftzb, rab, n1_cab, (double *)cab,
                            pabval, habval, forces, virials);
              } else {
                process_normal(a, b, 1.0, ftza, ftzb, rab, n1_cab,
                               (double *)cab, pabval, habval, forces, virials);
              }

              // Fill hdab and a_hdab.
              if (hdab != NULL) {
                assert(!compute_tau);
                process_force_a(a, b, 1.0, ftza, n1_cab, (double *)cab,
                                hdab[o2 + idx(b)][o1 + idx(a)]);
              }
              if (a_hdab != NULL) {
                assert(!compute_tau);
                process_virial_a(a, b, 1.0, ftza, n1_cab, (double *)cab,
                                 a_hdab[o2 + idx(b)][o1 + idx(a)]);
              }
            }
          }
        }
      }
    }
  }
}

// EOF

/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <string.h>

#include "grid_common.h"

// *****************************************************************************
static void grid_prepare_alpha(const double ra[3],
                               const double rb[3],
                               const double rp[3],
                               const int la_max,
                               const int lb_max,
                               double alpha[3][lb_max+1][la_max+1][la_max+lb_max+1]) {

    // Initialize with zeros.
    for (int iaxis=0; iaxis<3; iaxis++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxp=0; lxp<=la_max+lb_max; lxp++) {
        alpha[iaxis][lxb][lxa][lxp] = 0.0;
    }
    }
    }
    }

    //
    //   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls} alpha(ls,lxa,lxb,1)*(x-p)**ls
    //

    for (int iaxis=0; iaxis<3; iaxis++) {
       const double drpa = rp[iaxis] - ra[iaxis];
       const double drpb = rp[iaxis] - rb[iaxis];
       for (int lxa=0; lxa<=la_max; lxa++) {
       for (int lxb=0; lxb<=lb_max; lxb++) {
          double binomial_k_lxa = 1.0;
          double a = 1.0;
          for (int k=0; k<=lxa; k++) {
             double binomial_l_lxb = 1.0;
             double b = 1.0;
             for (int l=0; l<=lxb; l++) {
                alpha[iaxis][lxb][lxa][lxa-l+lxb-k] += binomial_k_lxa * binomial_l_lxb * a * b;
                binomial_l_lxb *= ((double)(lxb - l)) / ((double)(l + 1));
                b *= drpb;
             }
             binomial_k_lxa *= ((double)(lxa-k)) / ((double)(k+1));
             a *= drpa;
          }
       }
       }
    }
}

// *****************************************************************************
static void grid_prepare_coef(const int la_max,
                      const int la_min,
                      const int lb_max,
                      const int lb_min,
                      const double prefactor,
                      const double alpha[3][lb_max+1][la_max+1][la_max+lb_max+1],
                      const double pab[ncoset[lb_max]][ncoset[la_max]],
                      const int ncoef_xyz,
                      double coef_xyz[ncoef_xyz]) {


    const int lp = la_max + lb_max;

    double coef_xyt[((lp+1)*(lp+2))/2];
    double coef_xtt[lp+1];

    int lxyz = 0;
    for (int lzp = 0; lzp<=lp; lzp++) {
    for (int lyp = 0; lyp<=lp-lzp; lyp++) {
    for (int lxp = 0; lxp<=lp-lzp-lyp; lxp++) {
       coef_xyz[lxyz++] = 0.0;
    }
    }
    }

    for (int lzb = 0; lzb<=lb_max; lzb++) {
    for (int lza = 0; lza<=la_max; lza++) {
       int lxy = 0;
       for (int lyp = 0; lyp<=lp-lza-lzb; lyp++) {
          for (int lxp = 0; lxp<=lp-lza-lzb-lyp; lxp++) {
             coef_xyt[lxy++] = 0.0;
          }
          lxy = lxy + lza + lzb;
       }
       for (int lyb = 0; lyb<=lb_max-lzb; lyb++) {
       for (int lya = 0; lya<=la_max-lza; lya++) {
          const int lxpm = (lb_max-lzb-lyb) + (la_max-lza-lya);
          for (int i=0; i<=lxpm; i++) {
              coef_xtt[i] = 0.0;
          }
          for (int lxb = max(lb_min-lzb-lyb, 0); lxb<=lb_max-lzb-lyb; lxb++) {
          for (int lxa = max(la_min-lza-lya, 0); lxa<=la_max-lza-lya; lxa++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);
             const double p_ele = prefactor * pab[jco][ico];
             for (int lxp = 0; lxp<=lxa+lxb; lxp++) {
                coef_xtt[lxp] += p_ele * alpha[0][lxb][lxa][lxp];
             }
          }
          }
          int lxy = 0;
          for (int lyp = 0; lyp<=lya+lyb; lyp++) {
             for (int lxp = 0; lxp<=lp-lza-lzb-lya-lyb; lxp++) {
                coef_xyt[lxy++] += alpha[1][lyb][lya][lyp] * coef_xtt[lxp];
             }
             lxy += lza + lzb + lya + lyb - lyp;
          }
       }
       }
       lxyz = 0;
       for (int lzp = 0; lzp<=lza+lzb; lzp++) {
          int lxy = 0;
          for (int lyp = 0; lyp<=lp-lza-lzb; lyp++) {
             for (int lxp = 0; lxp<=lp-lza-lzb-lyp; lxp++) {
                coef_xyz[lxyz++] += alpha[2][lzb][lza][lzp] * coef_xyt[lxy++];
             }
             lxy += lza + lzb;
             lxyz += lza + lzb - lzp;
          }
          for (int lyp = lp-lza-lzb+1; lyp<=lp-lzp; lyp++) {
             for (int lxp = 0; lxp<=lp-lyp-lzp; lxp++) {
                lxyz++;
             }
          }
       }
    }
    }
}

// *****************************************************************************
static void grid_fill_map(const bool periodic,
                          const int lb_cube,
                          const int ub_cube,
                          const int cubecenter,
                          const int lb_grid,
                          const int npts,
                          const int ngrid,
                          const int cmax,
                          int map[2*cmax + 1]) {

    if (periodic) {
         int start = lb_cube;
         while (true) {
            const int offset = mod(cubecenter + start, npts)  + 1 - start;
            const int length = min(ub_cube, npts - offset) - start;
            for (int ig=start; ig<=start+length; ig++) {
               map[ig + cmax] = ig + offset;
            }
            if (start + length >= ub_cube){
                break;
            }
            start += length + 1;
         }
    } else {
         // this takes partial grid + border regions into account
         const int offset = mod(cubecenter + lb_cube + lb_grid, npts) + 1 - lb_cube;
         // check for out of bounds
         assert(ub_cube + offset <= ngrid);
         assert(lb_cube + offset >= 1);
         for (int ig=lb_cube; ig <= ub_cube; ig++) {
            map[ig + cmax] = ig + offset;
         }
    }
}



// *****************************************************************************
static void grid_fill_pol(const double dr,
                          const double roffset,
                          const int lb_cube,
                          const int lp,
                          const int cmax,
                          const double zetp,
                          double pol[lp+1][2*cmax+1]) {
//
//   compute the values of all (x-xp)**lp*exp(..)
//
//  still requires the old trick:
//  new trick to avoid to many exps (reuse the result from the previous gridpoint):
//  exp( -a*(x+d)**2)=exp(-a*x**2)*exp(-2*a*x*d)*exp(-a*d**2)
//  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)
//
      const double t_exp_1 = exp(-zetp * pow(dr, 2));
      const double t_exp_2 = pow(t_exp_1, 2);

      double t_exp_min_1 = exp(-zetp * pow(+dr - roffset, 2));
      double t_exp_min_2 = exp(-2 * zetp * (+dr - roffset) * (-dr));
      for (int ig=0; ig >= lb_cube; ig--) {
          const double rpg = ig * dr - roffset;
          t_exp_min_1 *= t_exp_min_2 * t_exp_1;
          t_exp_min_2 *= t_exp_2;
          double pg = t_exp_min_1;
          // pg  = EXP(-zetp*rpg**2)
          for (int icoef=0; icoef<=lp; icoef++) {
              pol[icoef][ig+cmax] = pg;
              pg *= rpg;
          }
      }

      double t_exp_plus_1 = exp(-zetp * pow(-roffset,2));
      double t_exp_plus_2 = exp(-2 * zetp * (-roffset) * (+dr));
      for (int ig=0; ig >= lb_cube; ig--) {
          const double rpg = (1-ig) * dr - roffset;
          t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
          t_exp_plus_2 *= t_exp_2;
          double pg = t_exp_plus_1;
          // pg  = EXP(-zetp*rpg**2)
          for (int icoef=0; icoef<=lp; icoef++) {
              pol[icoef][1-ig+cmax] = pg;
              pg *= rpg;
          }
      }
}


// *****************************************************************************
static void grid_collocate_core(const int cycles,
                                const int lp,
                                const int cmax,
                                const double coef_xyz[(lp+1)*(lp+2)*(lp+3)/6],
                                const double pol[3][lp+1][2*cmax+1],
                                const int map[3][2*cmax+1],
                                const double dh[3][3],
                                const double radius,
                                const int lb_cube[3],
                                const int ub_cube[3],
                                const int ngrid[3],
                                double grid[ngrid[2]][ngrid[1]][ngrid[0]]) {

    // Create the full cube, ignoring periodicity for now.
    const int nz = ub_cube[2] - lb_cube[2] + 1;
    const int ny = ub_cube[1] - lb_cube[1] + 1;
    const int nx = ub_cube[0] - lb_cube[0] + 1;
    double cube[nz][ny][nx];
    memset(cube, 0, sizeof(cube));

    // These loops should be vectorized as possible.
    int lxyz = 0;
    for (int lzp=0; lzp <= lp; lzp++) {
    for (int lyp=0; lyp <= lp-lzp; lyp++) {
    for (int lxp=0; lxp <= lp-lzp-lyp; lxp++) {
    const double coef = coef_xyz[lxyz++];

        for (int k=0; k < nz; k++) {
        for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
            cube[k][j][i] += coef * pol[2][lzp][k + lb_cube[2] + cmax]
                                  * pol[1][lyp][j + lb_cube[1] + cmax]
                                  * pol[0][lxp][i + lb_cube[0] + cmax];
        }
        }
        }

    }
    }
    }

    //
    // Write cube back to large grid taking periodicity and radius into account.
    //

    for (int k=0; k < nz; k++) {
    const int k2 = map[2][k + lb_cube[2] + cmax];
        for (int j=0; j < ny; j++) {
           const int j2 = map[1][j + lb_cube[1] + cmax];
            for (int i=0; i < nx; i++) {
                const int i2 = map[0][i + lb_cube[0] + cmax];
                grid[k2-1][j2-1][i2-1] += cycles * cube[k][j][i]; // TODO: here we cheate a bit.
            }
        }
    }

    // Unfortunatelly, I've not yet manage to reproduce the sphere_bounds
    // exactly with an on-the-fly formular because they use quantized distances.
    // See init_cube_info() for details.
    //
    //const double drmin = min(dh[0][0], min(dh[1][1], dh[2][2]));
    //const double eff_radius = drmin * max(1, ceil(radius/drmin));
    //const double eff_radius2 = eff_radius*eff_radius;
    //
    //for (int k=0; k < nz; k++) {
    //    const int k2 = map[2][k + lb_cube[2] + cmax];
    //    const double dk = (k + lb_cube[2]) * dh[2][2];
    //    const double dk2 = dk*dk;
    //    for (int j=0; j < ny; j++) {
    //        const int j2 = map[1][j + lb_cube[1] + cmax];
    //        const double dj = (j + lb_cube[1]) * dh[1][1];
    //        const double djk2 = dj*dj + dk2;
    //        for (int i=0; i < nx; i++) {
    //            const int i2 = map[0][i + lb_cube[0] + cmax];
    //            const double di = (i + lb_cube[0]) * dh[0][0];
    //            const double cur_radius2 = djk2 + di*di;
    //            if (cur_radius2 < eff_radius2) {
    //                grid[k2-1][j2-1][i2-1] += cube[k][j][i];
    //            }
    //        }
    //    }
    //}

}


// *****************************************************************************
void grid_collocate_pgf_product_cuda(const int cycles,
                                     const bool use_ortho,
                                     const int func,
                                     const int la_max,
                                     const int la_min,
                                     const int lb_max,
                                     const int lb_min,
                                     const double zeta,
                                     const double zetb,
                                     const double rscale,
                                     const double dh[3][3],
                                     const double dh_inv[3][3],
                                     const double ra[3],
                                     const double rab[3],
                                     const int npts[3],
                                     const int ngrid[3],
                                     const int lb_grid[3],
                                     const bool periodic[3],
                                     const double radius,
                                     const int lb_cube[3],
                                     const int ub_cube[3],
                                     const int o1,
                                     const int o2,
                                     const int n1,
                                     const int n2,
                                     const double pab[n2][n1],
                                     double grid[ngrid[2]][ngrid[1]][ngrid[0]]){

    assert(use_ortho && "Non-ortho not yet implemented for CUDA.");
    assert(func==100 && "Non-GRID_FUNC_AB not yet implemented for CUDA.");

    const double zetp = zeta + zetb;
    const double f = zetb / zetp;
    const double rab2 = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
    const double prefactor = rscale * exp(-zeta * f * rab2);
    double rp[3], rb[3];
    for (int i=0; i<3; i++) {
        rp[i] = ra[i] + f * rab[i];
        rb[i] = ra[i] + rab[i];
    }

    // grid_prepare_pab_AB
    const int n1_prep = ncoset[la_max];
    const int n2_prep = ncoset[lb_max];
    double pab_prep[n2_prep][n1_prep];
    memset(pab_prep, 0, sizeof(pab_prep)); // TODO needed?

    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
       for (int lya=0; lya<=la_max-lxa; lya++) {
       for (int lyb=0; lyb<=lb_max-lxb; lyb++) {
          for (int lza=max(la_min-lxa-lya, 0); lza<=la_max-lxa-lya; lza++) {
          for (int lzb=max(lb_min-lxb-lyb, 0); lzb<=lb_max-lxb-lyb; lzb++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);
             pab_prep[jco][ico] = pab[o2+jco][o1+ico];
          }
          }
       }
       }
    }
    }

    //   *** initialise the coefficient matrix, we transform the sum
    //
    // sum_{lxa,lya,lza,lxb,lyb,lzb} P_{lxa,lya,lza,lxb,lyb,lzb} *
    //         (x-a_x)**lxa (y-a_y)**lya (z-a_z)**lza (x-b_x)**lxb (y-a_y)**lya (z-a_z)**lza
    //
    // into
    //
    // sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-p_x)**lxp (y-p_y)**lyp (z-p_z)**lzp
    //
    // where p is center of the product gaussian, and lp = la_max + lb_max
    // (current implementation is l**7)
    //

    double alpha[3][lb_max+1][la_max+1][la_max+lb_max+1];
    grid_prepare_alpha(ra,
                       rb,
                       rp,
                       la_max,
                       lb_max,
                       alpha);

    //
    //   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
    //   use a three step procedure
    //   we don't store zeros, so counting is done using lxyz,lxy in order to have
    //   contiguous memory access in collocate_fast.F
    //

    const int lp = la_max + lb_max;
    const int ncoef_xyz = (lp+1)*(lp+2)*(lp+3)/6;
    double coef_xyz[ncoef_xyz];
    grid_prepare_coef(la_max,
                      la_min,
                      lb_max,
                      lb_min,
                      prefactor,
                      alpha,
                      pab_prep,
                      ncoef_xyz,
                      coef_xyz);

   // *** position of the gaussian product
   //
   // this is the actual definition of the position on the grid
   // i.e. a point rp(:) gets here grid coordinates
   // MODULO(rp(:)/dr(:),npts(:))+1
   // hence (0.0,0.0,0.0) in real space is rsgrid%lb on the rsgrid ((1,1,1) on grid)

    // cubecenter(:) = FLOOR(MATMUL(dh_inv, rp))
    int cubecenter[3];
    for (int i=0; i<3; i++) {
        double dh_inv_rp = 0.0;
        for (int j=0; j<3; j++) {
            dh_inv_rp += dh_inv[j][i] * rp[j];
        }
        cubecenter[i] = floor(dh_inv_rp);
    }

    double roffset[3];
    for (int i=0; i<3; i++) {
        roffset[i] = rp[i] - ((double) cubecenter[i]) * dh[i][i];
    }

    //cmax = MAXVAL(ub_cube)
    int cmax = INT_MIN;
    for (int i=0; i<3; i++) {
        cmax = max(cmax, ub_cube[i]);
    }

    // a mapping so that the ig corresponds to the right grid point
    int map[3][2*cmax+1];
    for (int i=0; i<3; i++) {
        grid_fill_map(periodic[i],
                      lb_cube[i],
                      ub_cube[i],
                      cubecenter[i],
                      lb_grid[i],
                      npts[i],
                      ngrid[i],
                      cmax,
                      map[i]);
    }

    double pol[3][lp+1][2*cmax+1];
    for (int i=0; i<3; i++) {
        grid_fill_pol(dh[i][i], roffset[i], lb_cube[i], lp, cmax, zetp, pol[i]);
    }

    grid_collocate_core(cycles,
                        lp,
                        cmax,
                        coef_xyz,
                        pol,
                        map,
                        dh,
                        radius,
                        lb_cube,
                        ub_cube,
                        ngrid,
                        grid);
}
//EOF

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
#include <stdlib.h>

#include "grid_collocate_replay.h"
#include "grid_collocate_cpu.h"
#include "grid_prepare_pab.h"
#include "grid_common.h"
#include "grid_globals.h"

// *****************************************************************************
static Array2d allocate_array_2d(int n1, int n2) {
    const size_t sizeof_array = sizeof(double) * n1 * n2;
    Array2d array;
    array.data = malloc(sizeof_array);
    array.s1 = n2;
    memset(array.data, 0, sizeof_array);
    return array;
}

// *****************************************************************************
typedef struct {
   double* data;
   int s1, s2;
} Array3d;

#define Array3dAt(array, i, j, k) array.data[i*array.s1 + j*array.s2 + k]

static Array3d allocate_array_3d(int n1, int n2, int n3) {
    const size_t sizeof_array = sizeof(double) * n1 * n2 * n3;
    Array3d array;
    array.data = malloc(sizeof_array);
    array.s1 = n2*n3;
    array.s2 = n3;
    memset(array.data, 0, sizeof_array);
    return array;
}

// *****************************************************************************
typedef struct {
   double* data;
   int s1, s2, s3;
} Array4d;

#define Array4dAt(array, i, j, k, l) \
   array.data[i*array.s1 + j*array.s2 + k*array.s3 + l]

static Array4d allocate_array_4d(int n1, int n2, int n3, int n4) {
    const size_t sizeof_array = sizeof(double) * n1 * n2 * n3 * n4;
    Array4d array;
    array.data = malloc(sizeof_array);
    array.s1 = n2*n3*n4;
    array.s2 = n3*n4;
    array.s3 = n4;
    memset(array.data, 0, sizeof_array);
    return array;
}

// *****************************************************************************
static void grid_prepare_alpha(const double ra[3],
                               const double rb[3],
                               const double rp[3],
                               const int la_max,
                               const int lb_max,
                               Array4d alpha) {

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
                Array4dAt(alpha, iaxis, lxb, lxa, lxa-l+lxb-k) += binomial_k_lxa * binomial_l_lxb * a * b;
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
                              const int lp,
                              const double prefactor,
                              const Array4d alpha,
                              const Array2d pab,
                              Array3d coef_xyz) {

    double coef_xyt[lp+1][lp+1];
    double coef_xtt[lp+1];

    for (int lzb = 0; lzb<=lb_max; lzb++) {
    for (int lza = 0; lza<=la_max; lza++) {
       for (int lyp = 0; lyp<=lp-lza-lzb; lyp++) {
          for (int lxp = 0; lxp<=lp-lza-lzb-lyp; lxp++) {
             coef_xyt[lyp][lxp] = 0.0;
          }
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
             const double p_ele = prefactor * Array2dAt(pab, jco, ico);
             for (int lxp = 0; lxp<=lxa+lxb; lxp++) {
                coef_xtt[lxp] += p_ele * Array4dAt(alpha, 0, lxb, lxa, lxp);
             }
          }
          }
          for (int lyp = 0; lyp<=lya+lyb; lyp++) {
             for (int lxp = 0; lxp<=lp-lza-lzb-lya-lyb; lxp++) {
                coef_xyt[lyp][lxp] += Array4dAt(alpha, 1, lyb, lya, lyp) * coef_xtt[lxp];
             }
          }
       }
       }
       for (int lzp = 0; lzp<=lza+lzb; lzp++) {
          for (int lyp = 0; lyp<=lp-lza-lzb; lyp++) {
             for (int lxp = 0; lxp<=lp-lza-lzb-lyp; lxp++) {
                Array3dAt(coef_xyz, lzp, lyp, lxp) += Array4dAt(alpha, 2, lzb, lza, lzp) * coef_xyt[lyp][lxp];
             }
          }
       }
    }
    }
}

// *****************************************************************************
static void grid_fill_map(const bool periodic[3],
                          const int lb_cube[3],
                          const int ub_cube[3],
                          const int cubecenter[3],
                          const int lb_grid[3],
                          const int npts[3],
                          const int ngrid[3],
                          const int cmax,
                          int map[][3]) {

    for (int idir=0; idir<3; idir++) {
        if (periodic[idir]) {
             //for (int i=0; i <= 2*cmax; i++)
             //    map[i] = mod(cubecenter + i - cmax, npts) + 1;
             int start = lb_cube[idir];
             while (true) {
                const int offset = mod(cubecenter[idir] + start, npts[idir])  + 1 - start;
                const int length = min(ub_cube[idir], npts[idir] - offset) - start;
                for (int ig=start; ig<=start+length; ig++) {
                   map[ig + cmax][idir] = ig + offset;
                }
                if (start + length >= ub_cube[idir]){
                    break;
                }
                start += length + 1;
             }
        } else {
             // this takes partial grid + border regions into account
             const int offset = mod(cubecenter[idir] + lb_cube[idir] + lb_grid[idir], npts[idir]) + 1 - lb_cube[idir];
             // check for out of bounds
             assert(ub_cube[idir] + offset <= ngrid[idir]);
             assert(lb_cube[idir] + offset >= 1);
             for (int ig=lb_cube[idir]; ig <= ub_cube[idir]; ig++) {
                map[ig + cmax][idir] = ig + offset;
             }
        }
    }
}


// *****************************************************************************
static void grid_fill_pol(const double dh[3][3],
                          const double roffset[3],
                          const int lb_cube[3],
                          const int lp,
                          const double zetp,
                          Array3d pol) {
//
//   compute the values of all (x-xp)**lp*exp(..)
//
//  still requires the old trick:
//  new trick to avoid to many exps (reuse the result from the previous gridpoint):
//  exp( -a*(x+d)**2)=exp(-a*x**2)*exp(-2*a*x*d)*exp(-a*d**2)
//  exp(-2*a*(x+d)*d)=exp(-2*a*x*d)*exp(-2*a*d**2)
//
    for (int idir=0; idir<3; idir++) {
        const double dr = dh[idir][idir];
        const double t_exp_1 = exp(-zetp * pow(dr, 2));
        const double t_exp_2 = pow(t_exp_1, 2);

        double t_exp_min_1 = exp(-zetp * pow(+dr - roffset[idir], 2));
        double t_exp_min_2 = exp(-2 * zetp * (+dr - roffset[idir]) * (-dr));
        for (int ig=0; ig >= lb_cube[idir]; ig--) {
            const double rpg = ig * dr - roffset[idir];
            t_exp_min_1 *= t_exp_min_2 * t_exp_1;
            t_exp_min_2 *= t_exp_2;
            double pg = t_exp_min_1;
            // pg  = EXP(-zetp*rpg**2)
            for (int icoef=0; icoef<=lp; icoef++) {
                Array3dAt(pol, idir, icoef, ig-lb_cube[idir]) = pg;
                pg *= rpg;
            }
        }

        double t_exp_plus_1 = exp(-zetp * pow(-roffset[idir],2));
        double t_exp_plus_2 = exp(-2 * zetp * (-roffset[idir]) * (+dr));
        for (int ig=0; ig >= lb_cube[idir]; ig--) {
            const double rpg = (1-ig) * dr - roffset[idir];
            t_exp_plus_1 *= t_exp_plus_2 * t_exp_1;
            t_exp_plus_2 *= t_exp_2;
            double pg = t_exp_plus_1;
            // pg  = EXP(-zetp*rpg**2)
            for (int icoef=0; icoef<=lp; icoef++) {
                Array3dAt(pol, idir, icoef, 1-ig-lb_cube[idir]) = pg;
                pg *= rpg;
            }
        }
    }
}

// *****************************************************************************
static void grid_collocate_core(const int lp,
                                const int cmax,
                                const Array3d coef_xyz,
                                const Array3d pol,
                                const int map[][3],
                                const int lb_cube[3],
                                const int ub_cube[3],
                                const double dh[3][3],
                                const double dh_inv[3][3],
                                const double disr_radius,
                                const int ngrid[3],
                                double* grid) {

    // Create the full cube, ignoring periodicity for now.
    const int nz = ub_cube[2] - lb_cube[2] + 1;
    const int ny = ub_cube[1] - lb_cube[1] + 1;
    const int nx = ub_cube[0] - lb_cube[0] + 1;
    double cube[nz][ny][nx];   // Can be large, run with "ulimit -s unlimited".
    memset(cube, 0, nz*ny*nx*sizeof(double));

    // These loops should be as vectorized as possible.
    for (int lzp=0; lzp <= lp; lzp++) {
    for (int lyp=0; lyp <= lp; lyp++) {
    for (int lxp=0; lxp <= lp; lxp++) {
        for (int k=0; k < nz; k++) {
        for (int j=0; j < ny; j++) {
        for (int i=0; i < nx; i++) {
            cube[k][j][i] += Array3dAt(coef_xyz, lzp, lyp, lxp)
                             * Array3dAt(pol, 2, lzp, k)
                             * Array3dAt(pol, 1, lyp, j)
                             * Array3dAt(pol, 0, lxp, i);
        }
        }
        }
    }
    }
    }

    //
    // Write cube back to large grid taking periodicity and radius into account.
    //

    // The cube contains an even number of grid points in each direction and
    // collocation is always performed on a pair of two opposing grid points.
    // Hence, the points with index 0 and 1 are both assigned distance zero via
    // the formular distance=(2*index-1)/2.

    const int kgmin = ceil(-1e-8 - disr_radius * dh_inv[2][2]);
    for (int kg=kgmin; kg <= 1-kgmin; kg++) {
        const int k = map[kg + cmax][2];   // target location on the grid
        const int kd = (2*kg - 1) / 2;     // distance from center in grid points
        const double kr = kd * dh[2][2];   // distance from center in a.u.
        const double kremain = disr_radius * disr_radius - kr * kr;
        const int jgmin = ceil(-1e-8 - sqrt(max(0.0, kremain)) * dh_inv[1][1]);
        for (int jg=jgmin; jg <= 1-jgmin; jg++) {
            const int j = map[jg + cmax][1];  // target location on the grid
            const int jd = (2*jg - 1) / 2;    // distance from center in grid points
            const double jr = jd * dh[1][1];  // distance from center in a.u.
            const double jremain = kremain - jr * jr;
            const int igmin = ceil(-1e-8 - sqrt(max(0.0, jremain)) * dh_inv[0][0]);
            for (int ig=igmin; ig<=1-igmin; ig++) {
                const int i = map[ig + cmax][0];  // target location on the grid
                const double res = cube[kg - lb_cube[2]][jg - lb_cube[1]][ig - lb_cube[0]];
                //grid[k-1][j-1][i-1] += res;
                grid[(k-1)*ngrid[1]*ngrid[0] + (j-1)*ngrid[0] + i-1] += res;
            }
        }
    }
}

// *****************************************************************************
static void grid_collocate_ortho(const int lp,
                                 const double zetp,
                                 const Array3d coef_xyz,
                                 const double dh[3][3],
                                 const double dh_inv[3][3],
                                 const double rp[3],
                                 const int npts[3],
                                 const int lb_grid[3],
                                 const bool periodic[3],
                                 const double radius,
                                 const int ngrid[3],
                                 double* grid) {

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

    // Historically, the radius gets discretized.
    const double drmin = min(dh[0][0], min(dh[1][1], dh[2][2]));
    const double disr_radius = drmin * max(1, ceil(radius/drmin));

    int lb_cube[3], ub_cube[3];
    for (int i=0; i<3; i++) {
        lb_cube[i] = ceil(-1e-8 - disr_radius * dh_inv[i][i]);
        ub_cube[i] = 1 - lb_cube[i];
    }

    //cmax = MAXVAL(ub_cube)
    int cmax = INT_MIN;
    for (int i=0; i<3; i++) {
        cmax = max(cmax, ub_cube[i]);
    }

    // a mapping so that the ig corresponds to the right grid point
    const size_t sizeof_map = (2*cmax+1) * 3 * sizeof(int);
    int (*map)[3] = malloc(sizeof_map);
    grid_fill_map(periodic,
                  lb_cube,
                  ub_cube,
                  cubecenter,
                  lb_grid,
                  npts,
                  ngrid,
                  cmax,
                  map);

    Array3d pol = allocate_array_3d(3, lp+1, 2*cmax+1);
    grid_fill_pol(dh, roffset, lb_cube, lp, zetp, pol);

    grid_collocate_core(lp,
                        cmax,
                        coef_xyz,
                        pol,
                        map,
                        lb_cube,
                        ub_cube,
                        dh,
                        dh_inv,
                        disr_radius,
                        ngrid,
                        grid);
    free(map);
    free(pol.data);
}


// *****************************************************************************
static void grid_collocate_general(const int lp,
                                   const double zetp,
                                   const Array3d coef_xyz,
                                   const double dh[3][3],
                                   const double dh_inv[3][3],
                                   const double rp[3],
                                   const int npts[3],
                                   const int lb_grid[3],
                                   const bool periodic[3],
                                   const double radius,
                                   const int ngrid[3],
                                   double* grid) {

// Translated from collocate_general_opt()
//
// transform P_{lxp,lyp,lzp} into a P_{lip,ljp,lkp} such that
// sum_{lxp,lyp,lzp} P_{lxp,lyp,lzp} (x-x_p)**lxp (y-y_p)**lyp (z-z_p)**lzp =
// sum_{lip,ljp,lkp} P_{lip,ljp,lkp} (i-i_p)**lip (j-j_p)**ljp (k-k_p)**lkp
//

    // aux mapping array to simplify life
    //TODO instead of this map we could use 3D arrays like coef_xyz.
    int coef_map[lp+1][lp+1][lp+1];

    //TODO really needed?
    //coef_map = HUGE(coef_map)
    for (int lzp=0; lzp<=lp; lzp++) {
        for (int lyp=0; lyp<=lp; lyp++) {
            for (int lxp=0; lxp<=lp; lxp++) {
                coef_map[lzp][lyp][lxp] = INT_MAX;
            }
        }
    }

    int lxyz = 0;
    for (int lzp=0; lzp<=lp; lzp++) {
        for (int lyp=0; lyp<=lp-lzp; lyp++) {
            for (int lxp=0; lxp<=lp-lzp-lyp; lxp++) {
                coef_map[lzp][lyp][lxp] = ++lxyz;
            }
        }
    }

    // center in grid coords
    // gp = MATMUL(dh_inv, rp)
    double gp[3];
    for (int i=0; i<3; i++) {
        gp[i] = 0.0;
        for (int j=0; j<3; j++) {
            gp[i] += dh_inv[j][i] * rp[j];
        }
    }

    // transform using multinomials
    double hmatgridp[lp+1][3][3];
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            hmatgridp[0][j][i] = 1.0;
            for (int k=1; k<=lp; k++) {
                hmatgridp[k][j][i] = hmatgridp[k-1][j][i] * dh[j][i];
            }
        }
    }

    // zero coef_ijk
    const int ncoef_ijk = ((lp+1)*(lp+2)*(lp+3))/6;
    double coef_ijk[ncoef_ijk];
    for (int i=0; i<ncoef_ijk; i++) {
        coef_ijk[i] = 0.0;
    }

    const int lpx = lp;
    for (int klx=0; klx<=lpx; klx++) {
    for (int jlx=0; jlx<=lpx-klx; jlx++) {
    for (int ilx=0; ilx<=lpx-klx-jlx; ilx++) {
        const int lx = ilx + jlx + klx;
        const int lpy = lp - lx;
        for (int kly=0; kly<=lpy; kly++) {
        for (int jly=0; jly<=lpy-kly; jly++) {
        for (int ily=0; ily<=lpy-kly-jly; ily++) {
            const int ly = ily + jly + kly;
            const int lpz = lp - lx - ly;
            for (int klz=0; klz<=lpz; klz++) {
            for (int jlz=0; jlz<=lpz-klz; jlz++) {
            for (int ilz=0; ilz<=lpz-klz-jlz; ilz++) {
                const int lz = ilz + jlz + klz;
                const int il = ilx + ily + ilz;
                const int jl = jlx + jly + jlz;
                const int kl = klx + kly + klz;
                const int lijk= coef_map[kl][jl][il];
                coef_ijk[lijk-1] += Array3dAt(coef_xyz, lz, ly, lx) *
                   hmatgridp[ilx][0][0] * hmatgridp[jlx][1][0] * hmatgridp[klx][2][0] *
                   hmatgridp[ily][0][1] * hmatgridp[jly][1][1] * hmatgridp[kly][2][1] *
                   hmatgridp[ilz][0][2] * hmatgridp[jlz][1][2] * hmatgridp[klz][2][2] *
                   fac[lx] * fac[ly] * fac[lz] /
                   (fac[ilx] * fac[ily] * fac[ilz] * fac[jlx] * fac[jly] * fac[jlz] * fac[klx] * fac[kly] * fac[klz]);
            }
            }
            }
        }
        }
        }
    }
    }
    }

    // CALL return_cube_nonortho(cube_info, radius, index_min, index_max, rp)
    //
    // get the min max indices that contain at least the cube that contains a sphere around rp of radius radius
    // if the cell is very non-orthogonal this implies that many useless points are included
    // this estimate can be improved (i.e. not box but sphere should be used)
    int index_min[3], index_max[3];
    for (int idir=0; idir<3; idir++) {
        index_min[idir] = INT_MAX;
        index_max[idir] = INT_MIN;
    }
    for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
    for (int k=-1; k<=1; k++) {
       const double x = rp[0] + i * radius;
       const double y = rp[1] + j * radius;
       const double z = rp[2] + k * radius;
       for (int idir=0; idir<3; idir++) {
          const double resc = dh_inv[0][idir] * x + dh_inv[1][idir] * y + dh_inv[2][idir] * z;
          index_min[idir] = min(index_min[idir], floor(resc));
          index_max[idir] = max(index_max[idir], ceil(resc));
       }
    }
    }
    }

    int offset[3];
    for (int idir=0; idir<3; idir++) {
        offset[idir] = mod(index_min[idir] + lb_grid[idir], npts[idir]) + 1;
    }

    // go over the grid, but cycle if the point is not within the radius
    for (int k=index_min[2]; k<=index_max[2]; k++) {
       const double dk = k - gp[2];
       int k_index;
       if (periodic[2]) {
          k_index = mod(k, npts[2]) + 1;
       } else {
          k_index = k - index_min[2] + offset[2];
       }

       // zero coef_xyt
       const int ncoef_xyt = ((lp+1)*(lp+2))/2;
       double coef_xyt[ncoef_xyt];
       for (int i=0; i<ncoef_xyt; i++) {
           coef_xyt[i] = 0.0;
       }

       int lxyz = 0;
       double dkp = 1.0;
       for (int kl=0; kl<=lp; kl++) {
          int lxy = 0;
          for (int jl=0; jl<=lp-kl; jl++) {
             for (int il=0; il<=lp-kl-jl; il++) {
                coef_xyt[lxy++] += coef_ijk[lxyz++] * dkp;
             }
             lxy += kl;
          }
          dkp *= dk;
       }


       for (int j=index_min[1]; j<=index_max[1]; j++) {
          const double dj = j - gp[1];
          int j_index;
          if (periodic[1]) {
             j_index = mod(j, npts[1]) + 1;
          } else {
             j_index = j - index_min[1] + offset[1];
          }

          double coef_xtt[lp+1];
          for (int i=0; i<=lp; i++) {
              coef_xtt[i] = 0.0;
          }
          int lxy = 0;
          double djp = 1.0;
          for (int jl=0; jl<=lp; jl++) {
             for (int il=0; il<=lp-jl; il++) {
                coef_xtt[il] += coef_xyt[lxy++] * djp;
             }
             djp *= dj;
          }

          // find bounds for the inner loop
          // based on a quadratic equation in i
          // a*i**2+b*i+c=radius**2

          // v = pointj-gp(1)*hmatgrid(:, 1)
          // a = DOT_PRODUCT(hmatgrid(:, 1), hmatgrid(:, 1))
          // b = 2*DOT_PRODUCT(v, hmatgrid(:, 1))
          // c = DOT_PRODUCT(v, v)
          // d = b*b-4*a*(c-radius**2)
          double a=0.0, b=0.0, c=0.0;
          for (int i=0; i<3; i++) {
             const double pointk = dh[2][i] * dk;
             const double pointj = pointk + dh[1][i] * dj;
             const double v = pointj - gp[0] * dh[0][i];
             a += dh[0][i] * dh[0][i];
             b += 2.0 * v * dh[0][i];
             c += v * v;
          }
          double d = b * b -4 * a * (c - radius * radius);
          if (d < 0.0) {
             continue;
          }

          // prepare for computing -zetp*rsq
          d = sqrt(d);
          const int ismin = ceill((-b-d)/(2.0*a));
          const int ismax = floor((-b+d)/(2.0*a));
          a *= -zetp;
          b *= -zetp;
          c *= -zetp;
          const int i = ismin - 1;

          // the recursion relation might have to be done
          // from the center of the gaussian (in both directions)
          // instead as the current implementation from an edge
          double exp2i = exp((a * i + b) * i + c);
          double exp1i = exp(2.0 * a * i + a + b);
          const double exp0i = exp(2.0 * a);

          for (int i=ismin; i<=ismax; i++) {
             const double di = i - gp[0];

             // polynomial terms
             double res = 0.0;
             double dip = 1.0;
             for (int il=0; il<=lp; il++) {
                res += coef_xtt[il] * dip;
                dip *= di;
             }

             // the exponential recursion
             exp2i *= exp1i;
             exp1i *= exp0i;
             res *= exp2i;

             int i_index;
             if (periodic[0]) {
                i_index = mod(i, npts[0]) + 1;
             } else {
                i_index = i - index_min[0] + offset[0];
             }
             //grid[k_index-1][j_index-1][i_index-1] += res;
             grid[(k_index-1)*ngrid[1]*ngrid[0] + (j_index-1)*ngrid[0] + i_index-1] += res;
          }
       }
    }
}

// *****************************************************************************
static void grid_collocate_internal(const bool use_ortho,
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
                                    const int o1,
                                    const int o2,
                                    const int n1,
                                    const int n2,
                                    const double* pab,
                                    double* grid){

    const double zetp = zeta + zetb;
    const double f = zetb / zetp;
    const double rab2 = rab[0] * rab[0] + rab[1] * rab[1] + rab[2] * rab[2];
    const double prefactor = rscale * exp(-zeta * f * rab2);
    double rp[3], rb[3];
    for (int i=0; i<3; i++) {
        rp[i] = ra[i] + f * rab[i];
        rb[i] = ra[i] + rab[i];
    }

    int la_min_diff, la_max_diff, lb_min_diff, lb_max_diff;
    grid_prepare_get_ldiffs(func,
                            &la_min_diff, &la_max_diff,
                            &lb_min_diff, &lb_max_diff);

    const int la_min_prep = max(la_min + la_min_diff, 0);
    const int lb_min_prep = max(lb_min + lb_min_diff, 0);
    const int la_max_prep = la_max + la_max_diff;
    const int lb_max_prep = lb_max + lb_max_diff;
    Array2d pab_prep = allocate_array_2d(ncoset[lb_max_prep], ncoset[la_max_prep]);

    Array2d pab_orig = {.data=pab, .s1=n2};

    grid_prepare_pab(func, o1, o2, la_max, la_min, lb_max, lb_min,
                     zeta, zetb, pab_orig, pab_prep);

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

    Array4d alpha = allocate_array_4d(3, lb_max_prep+1, la_max_prep+1, la_max_prep+lb_max_prep+1);
    grid_prepare_alpha(ra, rb, rp, la_max_prep, lb_max_prep, alpha);

    //
    //   compute P_{lxp,lyp,lzp} given P_{lxa,lya,lza,lxb,lyb,lzb} and alpha(ls,lxa,lxb,1)
    //   use a three step procedure
    //   we don't store zeros, so counting is done using lxyz,lxy in order to have
    //   contiguous memory access in collocate_fast.F
    //

    const int lp = la_max_prep + lb_max_prep;
    Array3d coef_xyz = allocate_array_3d(lp+1, lp+1, lp+1);
    grid_prepare_coef(la_max_prep,
                      la_min_prep,
                      lb_max_prep,
                      lb_min_prep,
                      lp,
                      prefactor,
                      alpha,
                      pab_prep,
                      coef_xyz);

    if (use_ortho) {
        grid_stats_add((Counters){.collocate_ortho_cpu = 1});
        grid_collocate_ortho(lp,
                             zetp,
                             coef_xyz,
                             dh,
                             dh_inv,
                             rp,
                             npts,
                             lb_grid,
                             periodic,
                             radius,
                             ngrid,
                             grid);
    } else {
        grid_stats_add((Counters){.collocate_general_cpu = 1});
        grid_collocate_general(lp,
                               zetp,
                               coef_xyz,
                               dh,
                               dh_inv,
                               rp,
                               npts,
                               lb_grid,
                               periodic,
                               radius,
                               ngrid,
                               grid);
    }

    free(alpha.data);
    free(coef_xyz.data);
}


// *****************************************************************************
void grid_collocate_pgf_product_cpu(const bool use_ortho,
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
                                    const int o1,
                                    const int o2,
                                    const int n1,
                                    const int n2,
                                    const double pab[n2][n1],
                                    double* grid){

// Uncomment this to dump all tasks to file.
//#define __GRID_DUMP_TASKS

#ifdef __GRID_DUMP_TASKS
    const size_t sizeof_grid = sizeof(double) * ngrid[0] * ngrid[1] * ngrid[2];
    double *grid_before = malloc(sizeof_grid);
    memcpy(grid_before, grid, sizeof_grid);
    memset(grid, 0, sizeof_grid);
#endif

    grid_collocate_internal(use_ortho,
                            func,
                            la_max,
                            la_min,
                            lb_max,
                            lb_min,
                            zeta,
                            zetb,
                            rscale,
                            dh,
                            dh_inv,
                            ra,
                            rab,
                            npts,
                            ngrid,
                            lb_grid,
                            periodic,
                            radius,
                            o1,
                            o2,
                            n1,
                            n2,
                            pab,
                            grid);

#ifdef __GRID_DUMP_TASKS

    grid_collocate_record(use_ortho,
                          func,
                          la_max,
                          la_min,
                          lb_max,
                          lb_min,
                          zeta,
                          zetb,
                          rscale,
                          dh,
                          dh_inv,
                          ra,
                          rab,
                          npts,
                          ngrid,
                          lb_grid,
                          periodic,
                          radius,
                          o1,
                          o2,
                          n1,
                          n2,
                          pab,
                          grid);

    for (int i=0; i < ngrid[0] * ngrid[1] * ngrid[2]; i++) {
        grid[i] += grid_before[i];
    }

    free(grid_before);
#endif

}

//EOF

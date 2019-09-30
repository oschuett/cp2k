/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

// *****************************************************************************
static int min(int x, int y) {
    return (x < y) ? x : y;
}

// *****************************************************************************
static int mod(int a, int m)
{
    return (a%m + m) % m;
}

// *****************************************************************************
int grid_fill_map(const bool periodic,
                  const int lb_cube,
                  const int ub_cube,
                  const int cubecenter,
                  const int lb_grid,
                  const int grid_lbound,
                  const int grid_ubound,
                  const int ng,
                  const int cmax,
                  int map[2*cmax + 1]) {

    if (periodic) {
         int start = lb_cube;
         while (true) {
            const int offset = mod(cubecenter + start, ng)  + 1 - start;
            const int length = min(ub_cube, ng - offset) - start;
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
         const int offset = mod(cubecenter + lb_cube + lb_grid, ng) + 1 - lb_cube;
         // check for out of bounds
         if (ub_cube + offset > grid_ubound || lb_cube + offset < grid_lbound){
             return -1;
         }
         for (int ig=lb_cube; ig <= ub_cube; ig++) {
            map[ig + cmax] = ig + offset;
         }
    }

    return 0;
}


// *****************************************************************************
int grid_fill_pol(
                  const double dr,
                  const double roffset,
                  const int lb_cube,
                  const int lp,
                  const int cmax,
                  const double zetp,
                  double pol[cmax+1][lp+1][2]
                  ) {
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
              pol[ig+cmax][icoef][0] = pg;
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
              pol[ig+cmax][icoef][1] = pg;
              pg *= rpg;
          }
      }

    return 0;
}

// *****************************************************************************
int grid_collocate_core(
                        const int grid_size_x,
                        const int grid_size_y,
                        const int grid_size_z,
                        const int grid_lbound_x,
                        const int grid_lbound_y,
                        const int grid_lbound_z,
                        const int lp,
                        const int cmax,
                        const double coef_xyz[(lp+1)*(lp+2)*(lp+3)/6],
                        const double pol[3][cmax+1][lp+1][2],
                        const int map[3][2*cmax+1],
                        const int *sphere_bounds,
                        double grid[grid_size_z][grid_size_y][grid_size_x]
                        ) {

    int sci = 0;

    const int kgmin = sphere_bounds[sci++];
    for (int kg=kgmin; kg <= 0; kg++) {
        const int kg2 = 1 - kg;
        const int k = map[2][kg + cmax];
        const int k2 = map[2][kg2 + cmax];

        // initialize coef_xy
        const int n_coef_xy = (lp+1)*(lp+2)/2;
        double coef_xy[n_coef_xy][2];
        for (int i=0; i < n_coef_xy; i++) {
            coef_xy[i][0] = 0.0;
            coef_xy[i][1] = 0.0;
        }

        int lxyz = 0;
        for (int lzp=0; lzp <= lp; lzp++) {
            int lxy = 0;
            for (int lyp=0; lyp <= lp-lzp; lyp++) {
                for (int lxp=0; lxp <= lp-lzp-lyp; lxp++) {
                    coef_xy[lxy][0] += coef_xyz[lxyz] * pol[2][kg+cmax][lzp][0];
                    coef_xy[lxy][1] += coef_xyz[lxyz] * pol[2][kg+cmax][lzp][1];
                    lxyz++;
                    lxy++;
                }
                lxy += lzp;
            }
        }

        const int jgmin = sphere_bounds[sci++];
        for (int jg=jgmin; jg <= 0; jg++) {
            const int jg2 = 1 - jg;
            const int j = map[1][jg + cmax];
            const int j2 = map[1][jg2 + cmax];

            // initialize coef_x
            double coef_x[lp+1][4];
            for (int i=0; i < lp+1; i++) {
                for (int j=0; j < 4; j++) {
                    coef_x[i][j] = 0.0;
                }
            }

            int lxy = 0;
            for (int lyp=0; lyp <= lp; lyp++) {
                for (int lxp=0; lxp <= lp-lyp; lxp++) {
                    coef_x[lxp][0] += coef_xy[lxy][0]*pol[1][jg+cmax][lyp][0];
                    coef_x[lxp][1] += coef_xy[lxy][1]*pol[1][jg+cmax][lyp][0];
                    coef_x[lxp][2] += coef_xy[lxy][0]*pol[1][jg+cmax][lyp][1];
                    coef_x[lxp][3] += coef_xy[lxy][1]*pol[1][jg+cmax][lyp][1];
                    lxy++;
                }
            }

            const int igmin = sphere_bounds[sci++];
            for (int ig=igmin; ig<=0; ig++) {
                const int ig2 = 1 - ig;
                const int i = map[0][ig + cmax];
                const int i2 = map[0][ig2 + cmax];

                double s01 = 0.0;
                double s02 = 0.0;
                double s03 = 0.0;
                double s04 = 0.0;
                double s05 = 0.0;
                double s06 = 0.0;
                double s07 = 0.0;
                double s08 = 0.0;

                for (int lxp=0; lxp <= lp; lxp++) {
                    s01 += coef_x[lxp][0]*pol[0][ig+cmax][lxp][0];
                    s02 += coef_x[lxp][1]*pol[0][ig+cmax][lxp][0];
                    s03 += coef_x[lxp][2]*pol[0][ig+cmax][lxp][0];
                    s04 += coef_x[lxp][3]*pol[0][ig+cmax][lxp][0];
                    s05 += coef_x[lxp][0]*pol[0][ig+cmax][lxp][1];
                    s06 += coef_x[lxp][1]*pol[0][ig+cmax][lxp][1];
                    s07 += coef_x[lxp][2]*pol[0][ig+cmax][lxp][1];
                    s08 += coef_x[lxp][3]*pol[0][ig+cmax][lxp][1];
                }

                grid[k-grid_lbound_z][j-grid_lbound_y][i-grid_lbound_x] += s01;
                grid[k2-grid_lbound_z][j-grid_lbound_y][i-grid_lbound_x] += s02;
                grid[k-grid_lbound_z][j2-grid_lbound_y][i-grid_lbound_x] += s03;
                grid[k2-grid_lbound_z][j2-grid_lbound_y][i-grid_lbound_x] += s04;
                grid[k-grid_lbound_z][j-grid_lbound_y][i2-grid_lbound_x] += s05;
                grid[k2-grid_lbound_z][j-grid_lbound_y][i2-grid_lbound_x] += s06;
                grid[k-grid_lbound_z][j2-grid_lbound_y][i2-grid_lbound_x] += s07;
                grid[k2-grid_lbound_z][j2-grid_lbound_y][i2-grid_lbound_x] += s08;
            }
        }
    }
    return 0;
}

//EOF

/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>

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
                        const double pol_x[2*cmax+1][lp+1],
                        const double pol_y[cmax+1][lp+1][2],
                        const double pol_z[cmax+1][lp+1][2],
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
                    //printf("OLE lxyz: %i, lxy: %i\n", lxyz, lxy);
                    coef_xy[lxy][0] += coef_xyz[lxyz] * pol_z[kg+cmax][lzp][0];
                    coef_xy[lxy][1] += coef_xyz[lxyz] * pol_z[kg+cmax][lzp][1];
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

            const int igmin = sphere_bounds[sci++];
            const int igmax = 1 - igmin;

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
                    coef_x[lxp][0] += coef_xy[lxy][0]*pol_y[jg+cmax][lyp][0];
                    coef_x[lxp][1] += coef_xy[lxy][1]*pol_y[jg+cmax][lyp][0];
                    coef_x[lxp][2] += coef_xy[lxy][0]*pol_y[jg+cmax][lyp][1];
                    coef_x[lxp][3] += coef_xy[lxy][1]*pol_y[jg+cmax][lyp][1];
                    lxy++;
                }
            }

            for (int ig=igmin; ig<=igmax; ig++) {
                const int i = map[0][ig + cmax];
                double s01 = 0.0;
                double s02 = 0.0;
                double s03 = 0.0;
                double s04 = 0.0;
                for (int lxp=0; lxp <= lp; lxp++) {
                    s01 += coef_x[lxp][0]*pol_x[ig+cmax][lxp];
                    s02 += coef_x[lxp][1]*pol_x[ig+cmax][lxp];
                    s03 += coef_x[lxp][2]*pol_x[ig+cmax][lxp];
                    s04 += coef_x[lxp][3]*pol_x[ig+cmax][lxp];
                }
                grid[k-grid_lbound_z][j-grid_lbound_y][i-grid_lbound_x] += s01;
                grid[k-grid_lbound_z][j2-grid_lbound_y][i-grid_lbound_x] += s03;
                grid[k2-grid_lbound_z][j-grid_lbound_y][i-grid_lbound_x] += s02;
                grid[k2-grid_lbound_z][j2-grid_lbound_y][i-grid_lbound_x] += s04;
            }
        }
    }
    return 0;
}

//EOF

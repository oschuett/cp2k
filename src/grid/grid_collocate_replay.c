/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#define _XOPEN_SOURCE 700   /* Enable POSIX 2008/13 */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#include "grid_collocate_replay.h"
#include "grid_collocate_cpu.h"


// *****************************************************************************
void grid_collocate_record(const bool use_ortho,
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
                           const int nspheres,
                           const int sphere_bounds[nspheres],
                           const int o1,
                           const int o2,
                           const int n1,
                           const int n2,
                           const double pab[n2][n1],
                           const double grid[ngrid[2]][ngrid[1]][ngrid[0]]){

    static int counter = 0;
    counter++;
    char filename[100];
    snprintf(filename, sizeof(filename), "grid_collocate_%05i.task", counter);

    const int D = DECIMAL_DIG;  // In C11 we could use DBL_DECIMAL_DIG.
    FILE *fp = fopen(filename, "w+");
    fprintf(fp, "#Grid collocate task v3\n");
    fprintf(fp, "use_ortho %i\n", use_ortho);
    fprintf(fp, "func %i\n", func);
    fprintf(fp, "la_max %i\n", la_max);
    fprintf(fp, "la_min %i\n", la_min);
    fprintf(fp, "lb_max %i\n", lb_max);
    fprintf(fp, "lb_min %i\n", lb_min);
    fprintf(fp, "zeta %.*e\n", D, zeta);
    fprintf(fp, "zetb %.*e\n", D, zetb);
    fprintf(fp, "rscale %.*e\n", D, rscale);
    for (int i=0; i<3; i++)
        fprintf(fp, "dh %i %.*e %.*e %.*e\n", i, D, dh[i][0], D, dh[i][1], D, dh[i][2]);
    for (int i=0; i<3; i++)
        fprintf(fp, "dh_inv %i %.*e %.*e %.*e\n", i, D, dh_inv[i][0], D, dh_inv[i][1], D, dh_inv[i][2]);
    fprintf(fp, "ra %.*e %.*e %.*e\n", D, ra[0], D, ra[1], D, ra[2]);
    fprintf(fp, "rab %.*e %.*e %.*e\n", D, rab[0], D, rab[1], D, rab[2]);
    fprintf(fp, "npts %i %i %i\n", npts[0], npts[1], npts[2]);
    fprintf(fp, "ngrid %i %i %i\n", ngrid[0], ngrid[1], ngrid[2]);
    fprintf(fp, "lb_grid %i %i %i\n", lb_grid[0], lb_grid[1], lb_grid[2]);
    fprintf(fp, "periodic %i %i %i\n", periodic[0], periodic[1], periodic[2]);
    fprintf(fp, "radius %.*e\n", D, radius);
    if (use_ortho) {
        fprintf(fp, "lb_cube %i %i %i\n", lb_cube[0], lb_cube[1], lb_cube[2]);
        fprintf(fp, "ub_cube %i %i %i\n", ub_cube[0], ub_cube[1], ub_cube[2]);
    }
    fprintf(fp, "nspheres %i\n", nspheres);

    int nspheres_nonzero = 0;
    for (int i=0; i<nspheres; i++) {
        if (sphere_bounds[i] != 0) {
            nspheres_nonzero++;
        }
    }
    fprintf(fp, "nspheres_nonzero %i\n", nspheres_nonzero);

    for (int i=0; i<nspheres; i++) {
        if (sphere_bounds[i] != 0) {
            fprintf(fp, "sphere_bounds %i %i\n", i, sphere_bounds[i]);
        }
    }

    fprintf(fp, "o1 %i\n", o1);
    fprintf(fp, "o2 %i\n", o2);
    fprintf(fp, "n1 %i\n", n1);
    fprintf(fp, "n2 %i\n", n2);

    for (int i=0; i < n2; i++) {
    for (int j=0; j < n1; j++) {
        fprintf(fp, "pab %i %i %.*e\n", i, j, D, pab[i][j]);
    }
    }

    int ngrid_nonzero = 0;
    for (int i=0; i<ngrid[2]; i++) {
    for (int j=0; j<ngrid[1]; j++) {
    for (int k=0; j<ngrid[0]; j++) {
        if (grid[i][j][k] != 0.0) {
            ngrid_nonzero++;
        }
    }
    }
    }
    fprintf(fp, "ngrid_nonzero %i\n", ngrid_nonzero);

    for (int i=0; i<ngrid[2]; i++) {
    for (int j=0; j<ngrid[1]; j++) {
    for (int k=0; j<ngrid[0]; j++) {
        if (grid[i][j][k] != 0.0) {
            fprintf(fp, "grid %i %i %i %.*e\n", i, j, k, D, grid[i][j][k]);
        }
    }
    }
    }
    fprintf(fp, "#THE_END\n");
    fclose(fp);
    printf("Wrote %s\n", filename);

}

// *****************************************************************************
double grid_collocate_replay(const char* filename, const int cycles){
    printf("Task:     %s\n", filename);
    FILE *fp = fopen(filename, "r");
    assert(fp != NULL && "Could not open task file.");

    char line[100], key[100];

    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(strcmp(line, "#Grid collocate task v3\n") == 0);

    int use_ortho_i;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &use_ortho_i) == 2);
    assert(strcmp(key, "use_ortho") == 0);
    bool use_ortho = use_ortho_i;

    int func;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &func) == 2);
    assert(strcmp(key, "func") == 0);

    int la_max;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &la_max) == 2);
    assert(strcmp(key, "la_max") == 0);

    int la_min;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &la_min) == 2);
    assert(strcmp(key, "la_min") == 0);

    int lb_max;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &lb_max) == 2);
    assert(strcmp(key, "lb_max") == 0);

    int lb_min;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &lb_min) == 2);
    assert(strcmp(key, "lb_min") == 0);

    double zeta;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &zeta) == 2);
    assert(strcmp(key, "zeta") == 0);

    double zetb;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &zetb) == 2);
    assert(strcmp(key, "zetb") == 0);

    double rscale;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &rscale) == 2);
    assert(strcmp(key, "rscale") == 0);

    double dh[3][3];
    for (int i=0; i<3; i++) {
        int j;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %le %le %le", key, &j, &dh[i][0], &dh[i][1], &dh[i][2]) == 5);
        assert(strcmp(key, "dh") == 0 && i == j);
    }

    double dh_inv[3][3];
    for (int i=0; i<3; i++) {
        int j;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %le %le %le", key, &j, &dh_inv[i][0], &dh_inv[i][1], &dh_inv[i][2]) == 5);
        assert(strcmp(key, "dh_inv") == 0 && i == j);
    }

    double ra[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le %le %le", key, &ra[0], &ra[1], &ra[2]) == 4);
    assert(strcmp(key, "ra") == 0);

    double rab[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le %le %le", key, &rab[0], &rab[1], &rab[2]) == 4);
    assert(strcmp(key, "rab") == 0);

    int npts[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &npts[0], &npts[1], &npts[2]) == 4);
    assert(strcmp(key, "npts") == 0);

    int ngrid[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &ngrid[0], &ngrid[1], &ngrid[2]) == 4);
    assert(strcmp(key, "ngrid") == 0);

    int lb_grid[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &lb_grid[0], &lb_grid[1], &lb_grid[2]) == 4);
    assert(strcmp(key, "lb_grid") == 0);

    int periodic_i[3];
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i %i %i", key, &periodic_i[0], &periodic_i[1], &periodic_i[2]) == 4);
    assert(strcmp(key, "periodic") == 0);
    bool periodic[3] = {periodic_i[0], periodic_i[1], periodic_i[2]};

    double radius;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %le", key, &radius) == 2);
    assert(strcmp(key, "radius") == 0);

    int *lb_cube, *ub_cube;
    int lb_cube_arr[3], ub_cube_arr[3];
    if (use_ortho) {
        lb_cube = lb_cube_arr;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i %i", key, &lb_cube[0], &lb_cube[1], &lb_cube[2]) == 4);
        assert(strcmp(key, "lb_cube") == 0);

        ub_cube = ub_cube_arr;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i %i", key, &ub_cube[0], &ub_cube[1], &ub_cube[2]) == 4);
        assert(strcmp(key, "ub_cube") == 0);
    } else {
        lb_cube = NULL;
        ub_cube = NULL;
    }

    int nspheres;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &nspheres) == 2);
    assert(strcmp(key, "nspheres") == 0);

    int sphere_bounds[nspheres];
    for (int i=0; i < nspheres; i++) {
        sphere_bounds[i] = 0;
    }

    int nspheres_nonzero;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &nspheres_nonzero) == 2);
    assert(strcmp(key, "nspheres_nonzero") == 0);

    for (int i=0; i < nspheres_nonzero; i++) {
        int j, value;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i", key, &j, &value) == 3);
        assert(strcmp(key, "sphere_bounds") == 0);
        sphere_bounds[j] = value;
    }

    int o1;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &o1) == 2);
    assert(strcmp(key, "o1") == 0);

    int o2;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &o2) == 2);
    assert(strcmp(key, "o2") == 0);

    int n1;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &n1) == 2);
    assert(strcmp(key, "n1") == 0);

    int n2;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &n2) == 2);
    assert(strcmp(key, "n2") == 0);

    double pab[n2][n1];
    for (int i=0; i<n2; i++) {
    for (int j=0; j<n1; j++) {
        int i2, j2;
        double value;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i %le", key, &i2, &j2, &value) == 4);
        assert(strcmp(key, "pab") == 0 && i == i2 && j==j2);
        pab[i][j] = value;
    }
    }

    int ngrid_nonzero;
    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(sscanf(line, "%99s %i", key, &ngrid_nonzero) == 2);
    assert(strcmp(key, "ngrid_nonzero") == 0);

    double grid_ref[ngrid[2]][ngrid[1]][ngrid[0]];
    for (int i=0; i<ngrid[2]; i++) {
    for (int j=0; j<ngrid[1]; j++) {
    for (int k=0; j<ngrid[0]; j++) {
        grid_ref[i][j][k] = 0.0;
    }
    }
    }

    for (int n=0; n < ngrid_nonzero; n++) {
        int i, j, k;
        double value;
        assert(fgets(line, sizeof(line), fp) != NULL);
        assert(sscanf(line, "%99s %i %i %i %le", key, &i, &j, &k, &value) == 5);
        assert(strcmp(key, "grid") == 0);
        grid_ref[i][j][k] = value;
    }

    assert(fgets(line, sizeof(line), fp) != NULL);
    assert(strcmp(line, "#THE_END\n") == 0);

    double grid_test[ngrid[2]][ngrid[1]][ngrid[0]];
    for (int i=0; i<ngrid[2]; i++) {
    for (int j=0; j<ngrid[1]; j++) {
    for (int k=0; j<ngrid[0]; j++) {
        grid_test[i][j][k] = 0.0;
    }
    }
    }

    printf("Cycles:   %e\n", (float)cycles);

    struct timespec start_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);

    for (int i=0; i < cycles ; i++) {
        grid_collocate_pgf_product_cpu(use_ortho,
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
                                       lb_cube,
                                       ub_cube,
                                       nspheres,
                                       sphere_bounds,
                                       o1,
                                       o2,
                                       n1,
                                       n2,
                                       pab,
                                       grid_test);
    }

    struct timespec end_time;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);

    double max_diff = 0.0;

    for (int i=0; i<ngrid[2]; i++) {
    for (int j=0; j<ngrid[1]; j++) {
    for (int k=0; j<ngrid[0]; j++) {
        const double diff = fabs(grid_test[i][j][k] - cycles * grid_ref[i][j][k]);
        max_diff = fmax(max_diff, diff);
    }
    }
    }

    printf("Max diff: %le\n", max_diff);

    const double delta_sec = (end_time.tv_sec - start_time.tv_sec) + 1e-9 * (end_time.tv_nsec - start_time.tv_nsec);
    printf("Time:     %le sec\n", delta_sec);

    return max_diff;
}

//EOF

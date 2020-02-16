/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/
#ifndef GRID_COLLOCATE_CPU_H
#define GRID_COLLOCATE_CPU_H

#ifdef __cplusplus
extern "C" {
#endif

//******************************************************************************
// \brief Collocates a single task. A task consists of a pair of atoms each
//        with a position, Gaussian exponent, and a range of angular momentum.
//        This function then collocates all combinations of spherical harmonics.
//
// \param use_ortho     When true use the faster ortho algorithm.
// \param func          TODO
// \param func_idir     TODO
// \param l{a,b}_max    Max angular momentum to collocate for give atom.
// \param l{a,b}_min    Lowest angular momentum to collocate for give atom.
// \param zet_{a,b}     Gaussian's exponent of given atom.
// \param rscale        Prefactor to take density matrix symmetry in account.
// \param dh            Incremental grid matrix
// \param dh_inv        Inverse incremental grid matrix
// \param ra            Position of atom a.
// \param rab           Vector difference between position of atom a and atom b.
// \param npts          Global number of grid points in each direction.
// \param ngrid         Local number of grid points in each direction.
// \param lb_grid       Lower bounds of the grid.
// \param periodic      Whether simulation box is periodic in given direction.
// \param lmax          Global maximum angular moment.
// \param radius        Radius where Gaussian becomes small than threshold eps.
// \param o{1,2}        Offsets. The sub-block to be collocated starts at pab[o2][o1]
// \param n{1,2}        Dimensions of density matrix block pab.
// \param pab           The atom-pair's density matrix block P_{ab}
//
// \param grid The output grid array to collocate into.
//******************************************************************************
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
                                    double grid[ngrid[2]][ngrid[1]][ngrid[0]]);

#ifdef __cplusplus
}
#endif

#endif
//EOF

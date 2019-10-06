/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2019  CP2K developers group                         *
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

const int ncoset[] = {1,  // l=0
                      4,  // l=1
                      10, // l=2 ...
                      20, 35, 56, 84, 120, 165, 220, 286, 364,
                      455, 560, 680, 816, 969, 1140, 1330};

// *****************************************************************************
// Returns zero based indices.
static int coset(int lx, int ly, int lz) {
    const int l = lx + ly + lz;
    if (l==0) {
        return 0;
    } else {
        return ncoset[l-1] + ((l-lx) * (l-lx+1)) /2 + lz;
    }
}

// *****************************************************************************
static int min(int x, int y) {
    return (x < y) ? x : y;
}

// *****************************************************************************
static int max(int x, int y) {
    return (x > y) ? x : y;
}

// *****************************************************************************
static int mod(int a, int m)
{
    return (a%m + m) % m;
}

// *****************************************************************************
int grid_prepare_pab_tau(const int o1,
                         const int o2,
                         const int la_max,
                         const int la_min,
                         const int lb_max,
                         const int lb_min,
                         const int maxco,
                         const double zeta,
                         const double zetb,
                         const double pab[maxco][maxco],
                         double pab_tau[ncoset[lb_max+1]][ncoset[la_max+1]]) {

    // create a new pab_tau so that mapping pab_tau with pgf_a pgf_b
    // is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
    // (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})

    const int nla = ncoset[la_max+1];
    const int nlb = ncoset[lb_max+1];
    // ALLOCATE (pab_tau(nla, nlb)) //TODO: Move from Fortran

    // Initialize with zeros.
    for (int ico=0; ico<nla; ico++) {
    for (int jco=0; jco<nlb; jco++) {
        pab_tau[jco][ico] = 0.0;
    }
    }

    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
       for (int lya=0; lya<=la_max-lxa; lya++) {
       for (int lyb=0; lyb<=lb_max-lxb; lyb++) {
          for (int lza=max(la_min-lxa-lya, 0); lza<=la_max-lxa-lya; lza++) {
          for (int lzb=max(lb_min-lxb-lyb, 0); lzb<=lb_max-lxb-lyb; lzb++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);

             int ico_l, jco_l;
             // x  (all safe if lxa = 0, as the spurious added terms have zero prefactor)

             ico_l = coset(max(lxa-1, 0), lya, lza);
             jco_l = coset(max(lxb-1, 0), lyb, lzb);
             pab_tau[jco_l][ico_l] += lxa * lxb * pab[o2+jco][o1+ico];
             ico_l = coset(max(lxa-1, 0), lya, lza);
             jco_l = coset((lxb+1), lyb, lzb);
             pab_tau[jco_l][ico_l] += -2.0 * lxa * zetb * pab[o2+jco][o1+ico];
             ico_l = coset((lxa+1), lya, lza);
             jco_l = coset(max(lxb-1, 0), lyb, lzb);
             pab_tau[jco_l][ico_l] += -2.0 * zeta * lxb * pab[o2+jco][o1+ico];
             ico_l = coset((lxa+1), lya, lza);
             jco_l = coset((lxb+1), lyb, lzb);
             pab_tau[jco_l][ico_l] += 4.0 * zeta * zetb * pab[o2+jco][o1+ico];

             // y

             ico_l = coset(lxa, max(lya-1, 0), lza);
             jco_l = coset(lxb, max(lyb-1, 0), lzb);
             pab_tau[jco_l][ico_l] += lya * lyb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, max(lya-1, 0), lza);
             jco_l = coset(lxb, (lyb+1), lzb);
             pab_tau[jco_l][ico_l] += -2.0 * lya * zetb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, (lya+1), lza);
             jco_l = coset(lxb, max(lyb-1, 0), lzb);
             pab_tau[jco_l][ico_l] += -2.0 * zeta * lyb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, (lya+1), lza);
             jco_l = coset(lxb, (lyb+1), lzb);
             pab_tau[jco_l][ico_l] += 4.0 * zeta * zetb * pab[o2+jco][o1+ico];

             // z

             ico_l = coset(lxa, lya, max(lza-1, 0));
             jco_l = coset(lxb, lyb, max(lzb-1, 0));
             pab_tau[jco_l][ico_l] += lza * lzb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, lya, max(lza-1, 0));
             jco_l = coset(lxb, lyb, (lzb+1));
             pab_tau[jco_l][ico_l] += -2.0 * lza * zetb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, lya, (lza+1));
             jco_l = coset(lxb, lyb, max(lzb-1, 0));
             pab_tau[jco_l][ico_l] += -2.0 * zeta * lzb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, lya, (lza+1));
             jco_l = coset(lxb, lyb, (lzb+1));
             pab_tau[jco_l][ico_l] += 4.0 * zeta * zetb * pab[o2+jco][o1+ico];
          }
          }
       }
       }
    }
    }

    // Divide by two.
    //TODO Maybe divide all prefactors above instead.
    for (int ico=0; ico<nla; ico++) {
    for (int jco=0; jco<nlb; jco++) {
        pab_tau[jco][ico] *= 0.5;
    }
    }

    return 0;
}

// *****************************************************************************
int grid_prepare_pab_rho(const int o1,
                         const int o2,
                         const int la_max,
                         const int la_min,
                         const int lb_max,
                         const int lb_min,
                         const int maxco,
                         const double pab[maxco][maxco],
                         double pab_rho[ncoset[lb_max]][ncoset[la_max]]) {

    const int nla = ncoset[la_max];
    const int nlb = ncoset[lb_max];
    //   ALLOCATE (pab_rho(nla, nlb)) // TODO move from Fortran.

    // Initialize with zeros.
    for (int ico=0; ico<nla; ico++) {
    for (int jco=0; jco<nlb; jco++) {
        pab_rho[jco][ico] = 0.0;
    }
    }
    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
       for (int lya=0; lya<=la_max-lxa; lya++) {
       for (int lyb=0; lyb<=lb_max-lxb; lyb++) {
          for (int lza=max(la_min-lxa-lya, 0); lza<=la_max-lxa-lya; lza++) {
          for (int lzb=max(lb_min-lxb-lyb, 0); lzb<=lb_max-lxb-lyb; lzb++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);
             pab_rho[jco][ico] = pab[o2+jco][o1+ico];
          }
          }
       }
       }
    }
    }

    return 0;
}

// *****************************************************************************
int grid_prepare_alpha(const double ra[3],
                       const double rb[3],
                       const double rp[3],
                       const int la_max,
                       const int lb_max,
                       const int lmax,
                       double alpha[3][lmax+1][lmax+1][2*lmax+1]) {

    // Initialize with zeros.
    for (int iaxis=0; iaxis<3; iaxis++) {
    for (int lxb=0; lxb<=lmax; lxb++) {
    for (int lxa=0; lxa<=lmax; lxa++) {
    for (int lxp=0; lxp<=2*lmax; lxp++) {
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

    return 0;
}

// *****************************************************************************
int grid_prepare_coef(const int la_max,
                      const int la_min,
                      const int lb_max,
                      const int lb_min,
                      const int lmax,
                      const double prefactor,
                      const double alpha[3][lmax+1][lmax+1][2*lmax+1],
                      const double pab[ncoset[lb_max]][ncoset[la_max]],
                      double coef_xyz[]) {


    const int lp = la_max + lb_max;

    const int n = ((lmax*2+1)*(lmax*2+2))/2;
    double coef_xyt[n];
    double coef_xtt[lmax*2 + 1];

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

    return 0;
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

/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2020  CP2K developers group                         *
 *****************************************************************************/

#include "grid_prepare_pab.h"
#include <stdbool.h>
#include <assert.h>

// *****************************************************************************
static void grid_prepare_pab_DADB(const int o1,
                                  const int o2,
                                  const int la_max,
                                  const int la_min,
                                  const int lb_max,
                                  const int lb_min,
                                  const double zeta,
                                  const double zetb,
                                  const int n1,
                                  const int n2,
                                  const double pab[n2][n1],
                                  const int n1_prep,
                                  const int n2_prep,
                                  double pab_prep[n2_prep][n1_prep]) {

    // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
    // is equivalent to mapping pab with 0.5 * (nabla pgf_a) . (nabla pgf_b)
    // (ddx pgf_a ) (ddx pgf_b) = (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x})*(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x})

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
             pab_prep[jco_l][ico_l] +=  0.5 * lxa * lxb * pab[o2+jco][o1+ico];
             ico_l = coset(max(lxa-1, 0), lya, lza);
             jco_l = coset((lxb+1), lyb, lzb);
             pab_prep[jco_l][ico_l] += -1.0 * lxa * zetb * pab[o2+jco][o1+ico];
             ico_l = coset((lxa+1), lya, lza);
             jco_l = coset(max(lxb-1, 0), lyb, lzb);
             pab_prep[jco_l][ico_l] += -1.0 * zeta * lxb * pab[o2+jco][o1+ico];
             ico_l = coset((lxa+1), lya, lza);
             jco_l = coset((lxb+1), lyb, lzb);
             pab_prep[jco_l][ico_l] += 2.0 * zeta * zetb * pab[o2+jco][o1+ico];

             // y

             ico_l = coset(lxa, max(lya-1, 0), lza);
             jco_l = coset(lxb, max(lyb-1, 0), lzb);
             pab_prep[jco_l][ico_l] +=  0.5 * lya * lyb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, max(lya-1, 0), lza);
             jco_l = coset(lxb, (lyb+1), lzb);
             pab_prep[jco_l][ico_l] += -1.0 * lya * zetb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, (lya+1), lza);
             jco_l = coset(lxb, max(lyb-1, 0), lzb);
             pab_prep[jco_l][ico_l] += -1.0 * zeta * lyb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, (lya+1), lza);
             jco_l = coset(lxb, (lyb+1), lzb);
             pab_prep[jco_l][ico_l] += 2.0 * zeta * zetb * pab[o2+jco][o1+ico];

             // z

             ico_l = coset(lxa, lya, max(lza-1, 0));
             jco_l = coset(lxb, lyb, max(lzb-1, 0));
             pab_prep[jco_l][ico_l] += 0.5 * lza * lzb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, lya, max(lza-1, 0));
             jco_l = coset(lxb, lyb, (lzb+1));
             pab_prep[jco_l][ico_l] += -1.0 * lza * zetb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, lya, (lza+1));
             jco_l = coset(lxb, lyb, max(lzb-1, 0));
             pab_prep[jco_l][ico_l] += -1.0 * zeta * lzb * pab[o2+jco][o1+ico];
             ico_l = coset(lxa, lya, (lza+1));
             jco_l = coset(lxb, lyb, (lzb+1));
             pab_prep[jco_l][ico_l] += 2.0 * zeta * zetb * pab[o2+jco][o1+ico];
          }
          }
       }
       }
    }
    }
}

// *****************************************************************************
static void grid_prepare_pab_ADBmDAB(const int idir,
                                     const int o1,
                                     const int o2,
                                     const int la_max,
                                     const int la_min,
                                     const int lb_max,
                                     const int lb_min,
                                     const double zeta,
                                     const double zetb,
                                     const int n1,
                                     const int n2,
                                     const double pab[n2][n1],
                                     const int n1_prep,
                                     const int n2_prep,
                                     double pab_prep[n2_prep][n1_prep]) {

     // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
     // is equivalent to mapping pab with
     //    pgf_a (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) pgf_b
     // ( pgf_a ) (ddx pgf_b) - (ddx pgf_a)( pgf_b ) =
     //          pgf_a *(lbx pgf_{b-1x} - 2*zetb*pgf_{b+1x}) -
     //                   (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_b

    assert(1 <= idir && idir <= 3);

    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
       for (int lya=0; lya<=la_max-lxa; lya++) {
       for (int lyb=0; lyb<=lb_max-lxb; lyb++) {
          for (int lza=max(la_min-lxa-lya, 0); lza<=la_max-lxa-lya; lza++) {
          for (int lzb=max(lb_min-lxb-lyb, 0); lzb<=lb_max-lxb-lyb; lzb++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);

             int ico_l, jco_l;

             // ! this element of pab results in 4 elements of pab_prep

             if (idir == 1) {  // x
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(max(lxb - 1, 0), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lxb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(max(lxa - 1, 0), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lxa*pab[o2 + jco][o1 + ico];
                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 2) {  // y
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, max(lyb - 1, 0), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lyb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, max(lya - 1, 0), lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lya*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else {  // z
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, max(lzb - 1, 0));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lzb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, max(lza - 1, 0));
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lza*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             }
          }
          }
       }
       }
    }
    }
}

// *****************************************************************************
static void grid_prepare_pab_DABpADB(const int idir,
                                     const int o1,
                                     const int o2,
                                     const int la_max,
                                     const int la_min,
                                     const int lb_max,
                                     const int lb_min,
                                     const double zeta,
                                     const double zetb,
                                     const int n1,
                                     const int n2,
                                     const double pab[n2][n1],
                                     const int n1_prep,
                                     const int n2_prep,
                                     double pab_prep[n2_prep][n1_prep]) {

    // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
    // is equivalent to mapping pab with
    //    pgf_a (nabla_{idir} pgf_b) + (nabla_{idir} pgf_a) pgf_b
    // ( pgf_a ) (ddx pgf_b) + (ddx pgf_a)( pgf_b ) =
    //          pgf_a *(lbx pgf_{b-1x} + 2*zetb*pgf_{b+1x}) +
    //                   (lax pgf_{a-1x} + 2*zeta*pgf_{a+1x}) pgf_b
    assert(1 <= idir && idir <= 3);

    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
       for (int lya=0; lya<=la_max-lxa; lya++) {
       for (int lyb=0; lyb<=lb_max-lxb; lyb++) {
          for (int lza=max(la_min-lxa-lya, 0); lza<=la_max-lxa-lya; lza++) {
          for (int lzb=max(lb_min-lxb-lyb, 0); lzb<=lb_max-lxb-lyb; lzb++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);

             int ico_l, jco_l;

             // this element of pab results in 4 elements of pab_prep

             if (idir == 1) {  // x
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(max(lxb - 1, 0), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lxb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(max(lxa - 1, 0), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lxa*pab[o2 + jco][o1 + ico];
                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 2) {  // y
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, max(lyb - 1, 0), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lyb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, max(lya - 1, 0), lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lya*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else { // z
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, max(lzb - 1, 0));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lzb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, max(lza - 1, 0));
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lza*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zeta*pab[o2 + jco][o1 + ico];
             }
          }
          }
       }
       }
    }
    }
}

// *****************************************************************************
static void grid_prepare_pab_ARDBmDARB(const int idir,
                                       const int ir,
                                       const int o1,
                                       const int o2,
                                       const int la_max,
                                       const int la_min,
                                       const int lb_max,
                                       const int lb_min,
                                       const double zeta,
                                       const double zetb,
                                       const int n1,
                                       const int n2,
                                       const double pab[n2][n1],
                                       const int n1_prep,
                                       const int n2_prep,
                                       double pab_prep[n2_prep][n1_prep]) {

    // create a new pab_prep so that mapping pab_prep with pgf_a pgf_b
    // is equivalent to mapping pab with
    // pgf_a (r-Rb)_{ir} (nabla_{idir} pgf_b) - (nabla_{idir} pgf_a) (r-Rb)_{ir}  pgf_b
    // ( pgf_a )(r-Rb)_{ir} (ddx pgf_b) - (ddx pgf_a) (r-Rb)_{ir} ( pgf_b ) =
    //                        pgf_a *(lbx pgf_{b-1x+1ir} - 2*zetb*pgf_{b+1x+1ir}) -
    //                       (lax pgf_{a-1x} - 2*zeta*pgf_{a+1x}) pgf_{b+1ir}


    assert(1 <= idir && idir <= 3);
    assert(1 <= ir && ir <= 3);

    for (int lxa=0; lxa<=la_max; lxa++) {
    for (int lxb=0; lxb<=lb_max; lxb++) {
       for (int lya=0; lya<=la_max-lxa; lya++) {
       for (int lyb=0; lyb<=lb_max-lxb; lyb++) {
          for (int lza=max(la_min-lxa-lya, 0); lza<=la_max-lxa-lya; lza++) {
          for (int lzb=max(lb_min-lxb-lyb, 0); lzb<=lb_max-lxb-lyb; lzb++) {
             const int ico = coset(lxa, lya, lza);
             const int jco = coset(lxb, lyb, lzb);

             int ico_l, jco_l;

             // this element of pab results in 4 elements of pab_prep

             if (idir == 1 && ir == 1) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lxb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 2), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(max(lxa - 1, 0), lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lxa*pab[o2 + jco][o1 + ico];
                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 1 && ir == 2) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(max(lxb - 1, 0), (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lxb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(max(lxa - 1, 0), lya, lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lxa*pab[o2 + jco][o1 + ico];
                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 1 && ir == 3) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(max(lxb - 1, 0), lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lxb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(max(lxa - 1, 0), lya, lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lxa*pab[o2 + jco][o1 + ico];
                ico_l = coset((lxa + 1), lya, lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 2 && ir == 1) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), max(lyb - 1, 0), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lyb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, max(lya - 1, 0), lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lya*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 2 && ir == 2) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lyb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 2), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, max(lya - 1, 0), lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lya*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 2 && ir == 3) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, max(lyb - 1, 0), (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lyb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 1), (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, max(lya - 1, 0), lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lya*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, (lya + 1), lza);
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 3 && ir == 1) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, max(lzb - 1, 0));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lzb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset((lxb + 1), lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, max(lza - 1, 0));
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lza*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset((lxb + 1), lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 3 && ir == 2) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 1), max(lzb - 1, 0));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lzb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, (lyb + 1), (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, max(lza - 1, 0));
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lza*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, (lyb + 1), lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             } else if (idir == 3 && ir == 3) {
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, lzb);
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + lzb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, lza);
                jco_l = coset(lxb, lyb, (lzb + 2));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - 2.0*zetb*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, max(lza - 1, 0));
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] - lza*pab[o2 + jco][o1 + ico];
                ico_l = coset(lxa, lya, (lza + 1));
                jco_l = coset(lxb, lyb, (lzb + 1));
                pab_prep[jco_l][ico_l] = pab_prep[jco_l][ico_l] + 2.0*zeta*pab[o2 + jco][o1 + ico];
             }
          }
          }
       }
       }
    }
    }
}

// *****************************************************************************
static void grid_prepare_pab_AB(const int o1,
                                const int o2,
                                const int la_max,
                                const int la_min,
                                const int lb_max,
                                const int lb_min,
                                const int n1,
                                const int n2,
                                const double pab[n2][n1],
                                const int n1_prep,
                                const int n2_prep,
                                double pab_prep[n2_prep][n1_prep]) {

    const int nla = ncoset[la_max];
    const int nlb = ncoset[lb_max];

    // Initialize with zeros.
    for (int ico=0; ico<nla; ico++) {
    for (int jco=0; jco<nlb; jco++) {
        pab_prep[jco][ico] = 0.0;
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
             pab_prep[jco][ico] = pab[o2+jco][o1+ico];
          }
          }
       }
       }
    }
    }
}

// *****************************************************************************
void grid_prepare_get_ldiffs(const int func,
                             int* la_min_diff,
                             int* la_max_diff,
                             int* lb_min_diff,
                             int* lb_max_diff) {
   switch(func) {
      case GRID_FUNC_DADB:
      case GRID_FUNC_ADBmDAB_X:
      case GRID_FUNC_ADBmDAB_Y:
      case GRID_FUNC_ADBmDAB_Z:
      case GRID_FUNC_DABpADB_X:
      case GRID_FUNC_DABpADB_Y:
      case GRID_FUNC_DABpADB_Z:
        *la_max_diff = +1;
        *la_min_diff = -1;
        *lb_max_diff = +1;
        *lb_min_diff = -1;
        break;
      case GRID_FUNC_ARDBmDARB_XX:
      case GRID_FUNC_ARDBmDARB_XY:
      case GRID_FUNC_ARDBmDARB_XZ:
      case GRID_FUNC_ARDBmDARB_YX:
      case GRID_FUNC_ARDBmDARB_YY:
      case GRID_FUNC_ARDBmDARB_YZ:
      case GRID_FUNC_ARDBmDARB_ZX:
      case GRID_FUNC_ARDBmDARB_ZY:
      case GRID_FUNC_ARDBmDARB_ZZ:
         *la_max_diff = +1;//TODO: mistake???
         *la_min_diff = -1;
         *lb_max_diff = +2;
         *lb_min_diff = -1;
         break;
      case GRID_FUNC_AB:
         *la_max_diff = 0;
         *la_min_diff = 0;
         *lb_max_diff = 0;
         *lb_min_diff = 0;
         break;
      default:
         assert(false && "Unknown ga_gb_function.");
    }
}

// *****************************************************************************
void grid_prepare_pab(const int func,
                      const int o1,
                      const int o2,
                      const int la_max,
                      const int la_min,
                      const int lb_max,
                      const int lb_min,
                      const double zeta,
                      const double zetb,
                      const int n1,
                      const int n2,
                      const double pab[n2][n1],
                      const int n1_prep,
                      const int n2_prep,
                      double pab_prep[n2_prep][n1_prep]) {

   switch(func) {
      case GRID_FUNC_DADB:
        grid_prepare_pab_DADB(o1, o2, la_max, la_min, lb_max, lb_min,
                              zeta, zetb,
                              n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ADBmDAB_X:
        grid_prepare_pab_ADBmDAB(1, o1, o2, la_max, la_min,
                                 lb_max, lb_min, zeta, zetb,
                                 n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ADBmDAB_Y:
        grid_prepare_pab_ADBmDAB(2, o1, o2, la_max, la_min,
                                 lb_max, lb_min, zeta, zetb,
                                 n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ADBmDAB_Z:
        grid_prepare_pab_ADBmDAB(3, o1, o2, la_max, la_min,
                                 lb_max, lb_min, zeta, zetb,
                                 n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_DABpADB_X:
        grid_prepare_pab_DABpADB(1, o1, o2, la_max, la_min,
                                 lb_max, lb_min, zeta, zetb,
                                 n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_DABpADB_Y:
        grid_prepare_pab_DABpADB(2, o1, o2, la_max, la_min,
                                 lb_max, lb_min, zeta, zetb,
                                 n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_DABpADB_Z:
        grid_prepare_pab_DABpADB(3, o1, o2, la_max, la_min,
                                 lb_max, lb_min, zeta, zetb,
                                 n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_XX:
        grid_prepare_pab_ARDBmDARB(1, 1, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_XY:
        grid_prepare_pab_ARDBmDARB(1, 2, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_XZ:
        grid_prepare_pab_ARDBmDARB(1, 3, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_YX:
        grid_prepare_pab_ARDBmDARB(2, 1, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_YY:
        grid_prepare_pab_ARDBmDARB(2, 2, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_YZ:
        grid_prepare_pab_ARDBmDARB(2, 3, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_ZX:
        grid_prepare_pab_ARDBmDARB(3, 1, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_ZY:
        grid_prepare_pab_ARDBmDARB(3, 2, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_ARDBmDARB_ZZ:
        grid_prepare_pab_ARDBmDARB(3, 3, o1, o2, la_max, la_min,
                                   lb_max, lb_min, zeta, zetb,
                                   n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      case GRID_FUNC_AB:
        grid_prepare_pab_AB(o1, o2, la_max, la_min, lb_max, lb_min,
                            n1, n2, pab, n1_prep, n2_prep, pab_prep);
        break;
      default:
        assert(false && "Unknown ga_gb_function.");
    }
}

//EOF

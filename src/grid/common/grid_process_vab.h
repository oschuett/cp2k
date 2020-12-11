/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__CUDACC__)
#define GRID_DEVICE __device__
#else
#define GRID_DEVICE
#endif

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the forces for atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_force_a(const orbital a, const orbital b, const double pab,
                const double ftza, const int m1, const int m2,
                const double vab[m2][m1], double force_a[3]) {

  for (int i = 0; i < 3; i++) {
    const double aip1 = vab[idx(b)][idx(up(i, a))];
    const double aim1 = vab[idx(b)][idx(down(i, a))];
    force_a[i] += pab * (ftza * aip1 - a.l[i] * aim1);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the forces for atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_force_b(const orbital a, const orbital b, const double pab,
                const double ftzb, const double rab[3], const int m1,
                const int m2, const double vab[m2][m1], double force_b[3]) {

  const double axpm0 = vab[idx(b)][idx(a)];
  for (int i = 0; i < 3; i++) {
    const double aip1 = vab[idx(b)][idx(up(i, a))];
    const double bim1 = vab[idx(down(i, b))][idx(a)];
    force_b[i] += pab * (ftzb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the virial for atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_virial_a(const orbital a, const orbital b, const double pab,
                 const double ftza, const int m1, const int m2,
                 const double vab[m2][m1], double virial_a[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      virial_a[i][j] += pab * ftza * vab[idx(b)][idx(up(i, up(j, a)))] -
                        pab * a.l[j] * vab[idx(b)][idx(up(i, down(j, a)))];
    }
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the virial for atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_virial_b(const orbital a, const orbital b, const double pab,
                 const double ftzb, const double rab[3], const int m1,
                 const int m2, const double vab[m2][m1],
                 double virial_b[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      virial_b[i][j] += pab * ftzb *
                            (vab[idx(b)][idx(up(i, up(j, a)))] -
                             vab[idx(b)][idx(up(i, a))] * rab[j] -
                             vab[idx(b)][idx(up(j, a))] * rab[i] +
                             vab[idx(b)][idx(a)] * rab[j] * rab[i]) -
                        pab * b.l[j] * vab[idx(up(i, down(j, b)))][idx(a)];
    }
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain forces and virials.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static void
process_normal(const orbital a, const orbital b, const double f,
               const double ftza, const double ftzb, const double rab[3],
               const int m1, const int m2, const double vab[m2][m1],
               const double pab, double *hab, double forces[2][3],
               double virials[2][3][3]) {

  *hab += f * vab[idx(b)][idx(a)];

  if (forces != NULL) {
    process_force_a(a, b, f * pab, ftza, m1, m2, vab, forces[0]);
    process_force_b(a, b, f * pab, ftzb, rab, m1, m2, vab, forces[1]);
  }

  if (virials != NULL) {
    process_virial_a(a, b, f * pab, ftza, m1, m2, vab, virials[0]);
    process_virial_b(a, b, f * pab, ftzb, rab, m1, m2, vab, virials[1]);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain forces and virials for tau.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static void
process_tau(const orbital a, const orbital b, const double ftza,
            const double ftzb, const double rab[3], const int m1, const int m2,
            const double vab[m2][m1], const double pab, double *hab,
            double forces[2][3], double virials[2][3][3]) {

  for (int i = 0; i < 3; i++) {
    process_normal(down(i, a), down(i, b), 0.5 * a.l[i] * b.l[i], ftza, ftzb,
                   rab, m1, m2, vab, pab, hab, forces, virials);
    process_normal(up(i, a), down(i, b), -0.5 * ftza * b.l[i], ftza, ftzb, rab,
                   m1, m2, vab, pab, hab, forces, virials);
    process_normal(down(i, a), up(i, b), -0.5 * a.l[i] * ftzb, ftza, ftzb, rab,
                   m1, m2, vab, pab, hab, forces, virials);
    process_normal(up(i, a), up(i, b), 0.5 * ftza * ftzb, ftza, ftzb, rab, m1,
                   m2, vab, pab, hab, forces, virials);
  }
}

// EOF

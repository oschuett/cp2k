/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <stdbool.h>

#if defined(__CUDACC__)
#define GRID_DEVICE __device__
#else
#define GRID_DEVICE
#endif

/*******************************************************************************
 * \brief Returns matrix element cab[idx(b)][idx(a)].
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double get_term(const orbital a, const orbital b,
                                          const int n, const double *cab) {
  return cab[idx(b) * n + idx(a)];
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the forces for atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void process_force_a(const orbital a, const orbital b,
                                               const double pab,
                                               const double ftza, const int n,
                                               const double *cab,
                                               double force_a[3]) {

  for (int i = 0; i < 3; i++) {
    const double aip1 = get_term(up(i, a), b, n, cab);
    const double aim1 = get_term(down(i, a), b, n, cab);
    force_a[i] += pab * (ftza * aip1 - a.l[i] * aim1);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the forces for atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_force_b(const orbital a, const orbital b, const double pab,
                const double ftzb, const double rab[3], const int n,
                const double *cab, double force_b[3]) {

  const double axpm0 = get_term(a, b, n, cab);
  for (int i = 0; i < 3; i++) {
    const double aip1 = get_term(up(i, a), b, n, cab);
    const double bim1 = get_term(a, down(i, b), n, cab);
    force_b[i] += pab * (ftzb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1);
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the virial for atom a.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_virial_a(const orbital a, const orbital b, const double pab,
                 const double ftza, const int n, const double *cab,
                 double virial_a[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      virial_a[i][j] += pab * ftza * get_term(up(i, up(j, a)), b, n, cab) -
                        pab * a.l[j] * get_term(up(i, down(j, a)), b, n, cab);
    }
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain the virial for atom b.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline void
process_virial_b(const orbital a, const orbital b, const double pab,
                 const double ftzb, const double rab[3], const int n,
                 const double *cab, double virial_b[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      virial_b[i][j] += pab * ftzb *
                            (get_term(up(i, up(j, a)), b, n, cab) -
                             get_term(up(i, a), b, n, cab) * rab[j] -
                             get_term(up(j, a), b, n, cab) * rab[i] +
                             get_term(a, b, n, cab) * rab[j] * rab[i]) -
                        pab * b.l[j] * get_term(a, up(i, down(j, b)), n, cab);
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
               const int n, const double *cab, const double pab, double *hab,
               double forces[2][3], double virials[2][3][3]) {

  *hab += f * get_term(a, b, n, cab);

  if (forces != NULL) {
    process_force_a(a, b, f * pab, ftza, n, cab, forces[0]);
    process_force_b(a, b, f * pab, ftzb, rab, n, cab, forces[1]);
  }

  if (virials != NULL) {
    process_virial_a(a, b, f * pab, ftza, n, cab, virials[0]);
    process_virial_b(a, b, f * pab, ftzb, rab, n, cab, virials[1]);
  }
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_force_a_normal(const orbital a, const orbital b, const int i,
                       const double zeta, const int n, const double *cab) {
  const double aip1 = get_term(up(i, a), b, n, cab);
  const double aim1 = get_term(down(i, a), b, n, cab);
  return 2.0 * zeta * aip1 - a.l[i] * aim1;
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_force_a(const orbital a, const orbital b, const int i,
                const double zeta, const double zetb, const int n,
                const double *cab, const bool compute_tau) {
  if (!compute_tau) {
    return extract_force_a_normal(a, b, i, zeta, n, cab);
  } else {
    double force = 0.0;
    for (int i = 0; i < 3; i++) {
      force += 0.5 * a.l[i] * b.l[i] *
               extract_force_a_normal(down(i, a), down(i, b), i, zeta, n, cab);
      force -= zeta * b.l[i] *
               extract_force_a_normal(up(i, a), down(i, b), i, zeta, n, cab);
      force -= a.l[i] * zetb *
               extract_force_a_normal(down(i, a), up(i, b), i, zeta, n, cab);
      force += 2.0 * zeta * zetb *
               extract_force_a_normal(up(i, a), up(i, b), i, zeta, n, cab);
    }
    return force;
  }
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_force_b_normal(const orbital a, const orbital b, const int i,
                       const double zetb, const double rab[3], const int n,
                       const double *cab) {
  const double axpm0 = get_term(a, b, n, cab);
  const double aip1 = get_term(up(i, a), b, n, cab);
  const double bim1 = get_term(a, down(i, b), n, cab);
  return 2.0 * zetb * (aip1 - rab[i] * axpm0) - b.l[i] * bim1;
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_force_b(const orbital a, const orbital b, const int i,
                const double zeta, const double zetb, const double rab[3],
                const int n, const double *cab, const bool compute_tau) {
  if (!compute_tau) {
    return extract_force_b_normal(a, b, i, zetb, rab, n, cab);
  } else {
    double force = 0.0;
    for (int i = 0; i < 3; i++) {
      force +=
          0.5 * a.l[i] * b.l[i] *
          extract_force_b_normal(down(i, a), down(i, b), i, zetb, rab, n, cab);
      force -=
          zeta * b.l[i] *
          extract_force_b_normal(up(i, a), down(i, b), i, zetb, rab, n, cab);
      force -=
          a.l[i] * zetb *
          extract_force_b_normal(down(i, a), up(i, b), i, zetb, rab, n, cab);
      force += 2.0 * zeta * zetb *
               extract_force_b_normal(up(i, a), up(i, b), i, zetb, rab, n, cab);
    }
    return force;
  }
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_virial_a_normal(const orbital a, const orbital b, const int i,
                        const int j, const double zeta, const int n,
                        const double *cab) {
  return 2.0 * zeta * get_term(up(i, up(j, a)), b, n, cab) -
         a.l[j] * get_term(up(i, down(j, a)), b, n, cab);
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_virial_a(const orbital a, const orbital b, const int i, const int j,
                 const double zeta, const double zetb, const int n,
                 const double *cab, const bool compute_tau) {

  if (!compute_tau) {
    return extract_virial_a_normal(a, b, i, j, zeta, n, cab);
  } else {
    double virial = 0.0;
    for (int i = 0; i < 3; i++) {
      virial +=
          0.5 * a.l[i] * b.l[i] *
          extract_virial_a_normal(down(i, a), down(i, b), i, j, zeta, n, cab);
      virial -=
          zeta * b.l[i] *
          extract_virial_a_normal(up(i, a), down(i, b), i, j, zeta, n, cab);
      virial -=
          a.l[i] * zetb *
          extract_virial_a_normal(down(i, a), up(i, b), i, j, zeta, n, cab);
      virial += 2.0 * zeta * zetb *
                extract_virial_a_normal(up(i, a), up(i, b), i, j, zeta, n, cab);
    }
    return virial;
  }
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_virial_b_normal(const orbital a, const orbital b, const int i,
                        const int j, const double zetb, const double rab[3],
                        const int n, const double *cab) {

  return 2.0 * zetb *
             (get_term(up(i, up(j, a)), b, n, cab) -
              get_term(up(i, a), b, n, cab) * rab[j] -
              get_term(up(j, a), b, n, cab) * rab[i] +
              get_term(a, b, n, cab) * rab[j] * rab[i]) -
         b.l[j] * get_term(a, up(i, down(j, b)), n, cab);
}

/*******************************************************************************
 * \brief TODO
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double
extract_virial_b(const orbital a, const orbital b, const int i, const int j,
                 const double zeta, const double zetb, const double rab[3],
                 const int n, const double *cab, const bool compute_tau) {

  if (!compute_tau) {
    return extract_virial_b_normal(a, b, i, j, zetb, rab, n, cab);
  } else {
    double virial = 0.0;
    for (int i = 0; i < 3; i++) {
      virial += 0.5 * a.l[i] * b.l[i] *
                extract_virial_b_normal(down(i, a), down(i, b), i, j, zetb, rab,
                                        n, cab);
      virial -= zeta * b.l[i] *
                extract_virial_b_normal(up(i, a), down(i, b), i, j, zetb, rab,
                                        n, cab);
      virial -= a.l[i] * zetb *
                extract_virial_b_normal(down(i, a), up(i, b), i, j, zetb, rab,
                                        n, cab);
      virial +=
          2.0 * zeta * zetb *
          extract_virial_b_normal(up(i, a), up(i, b), i, j, zetb, rab, n, cab);
    }
    return virial;
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain forces and virials for tau.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static inline double extract_hab(const orbital a, const orbital b,
                                             const double zeta,
                                             const double zetb, const int n,
                                             const double *cab,
                                             const bool compute_tau) {
  if (!compute_tau) {
    return get_term(a, b, n, cab);
  } else {
    double hab = 0.0;
    for (int i = 0; i < 3; i++) {
      hab += 0.5 * a.l[i] * b.l[i] * get_term(down(i, a), down(i, b), n, cab);
      hab -= zeta * b.l[i] * get_term(up(i, a), down(i, b), n, cab);
      hab -= a.l[i] * zetb * get_term(down(i, a), up(i, b), n, cab);
      hab += 2.0 * zeta * zetb * get_term(up(i, a), up(i, b), n, cab);
    }
    return hab;
  }
}

/*******************************************************************************
 * \brief Contracts given matrix elements to obtain forces and virials for tau.
 * \author Ole Schuett
 ******************************************************************************/
GRID_DEVICE static void process_tau(const orbital a, const orbital b,
                                    const double ftza, const double ftzb,
                                    const double rab[3], const int n,
                                    const double *cab, const double pab,
                                    double *hab, double forces[2][3],
                                    double virials[2][3][3]) {

  for (int i = 0; i < 3; i++) {
    process_normal(down(i, a), down(i, b), 0.5 * a.l[i] * b.l[i], ftza, ftzb,
                   rab, n, cab, pab, hab, forces, virials);
    process_normal(up(i, a), down(i, b), -0.5 * ftza * b.l[i], ftza, ftzb, rab,
                   n, cab, pab, hab, forces, virials);
    process_normal(down(i, a), up(i, b), -0.5 * a.l[i] * ftzb, ftza, ftzb, rab,
                   n, cab, pab, hab, forces, virials);
    process_normal(up(i, a), up(i, b), 0.5 * ftza * ftzb, ftza, ftzb, rab, n,
                   cab, pab, hab, forces, virials);
  }
}

/*******************************************************************************
 * \brief Differences in angular momentum.
 * \author Ole Schuett
 ******************************************************************************/
typedef struct {
  int la_max_diff;
  int la_min_diff;
  int lb_max_diff;
  int lb_min_diff;
} process_ldiffs;

/*******************************************************************************
 * \brief Returns difference in angular momentum range for given flags.
 * \author Ole Schuett
 ******************************************************************************/
static process_ldiffs process_get_ldiffs(bool calculate_forces,
                                         bool calculate_virial,
                                         bool compute_tau) {
  process_ldiffs ldiffs;

  ldiffs.la_max_diff = 0;
  ldiffs.lb_max_diff = 0;
  ldiffs.la_min_diff = 0;
  ldiffs.lb_min_diff = 0;

  if (calculate_forces || calculate_virial) {
    ldiffs.la_max_diff += 1; // for deriv. of gaussian, unimportant which one
    ldiffs.la_min_diff -= 1;
    ldiffs.lb_min_diff -= 1;
  }

  if (calculate_virial) {
    ldiffs.la_max_diff += 1;
    ldiffs.lb_max_diff += 1;
  }

  if (compute_tau) {
    ldiffs.la_max_diff += 1;
    ldiffs.lb_max_diff += 1;
    ldiffs.la_min_diff -= 1;
    ldiffs.lb_min_diff -= 1;
  }

  return ldiffs;
}

// EOF

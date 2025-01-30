#include "RungeKutta.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include <iostream>
#include <vector>

RungeKutta::RungeKutta(Hamiltonian H, Generator eta) {
  set_zero_body_flow_eq(H, eta);
  set_one_body_flow_eq(H, eta);
}

void RungeKutta::set_zero_body_flow_eq(Hamiltonian H, Generator eta) {
  double E = H.get_zero_body();
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();

  std::vector<OneBody> eta_1body = eta.get_one_body();
  std::vector<TwoBody> eta_2body = eta.get_two_body();

  // E flow equation
  dEds = 0;

  int eta_1body_size = eta_1body.size();

  int eta_2body_size = eta_2body.size();

#pragma omp parallel
  {
#pragma omp for nowait reduction(+ : dEds)
    for (int i = 1; i < eta_1body_size; ++i) {
      if (eta_1body[i].s < 4 && eta_1body[i].ss > 3) {
        for (const OneBody &ff : f) {
          if (ff.s == eta_1body[i].s && ff.ss == eta_1body[i].ss) {
            // n_a eta_ab f_ab - n_b (-eta_ab) f_ab
            dEds += 2 * eta_1body[i].me * ff.me;

            continue;
          }
        }
      }
    }
#pragma omp for nowait reduction(+ : dEds)
    for (int i = 0; i < eta_2body_size; ++i) {
      if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4 && eta_2body[i].ss1 > 3 &&
          eta_2body[i].ss2 > 3) {
        for (const TwoBody &g : Gamma) {
          if (g.s1 == eta_2body[i].s1 && g.s2 == eta_2body[i].s2 &&
              g.ss1 == eta_2body[i].ss1 && g.ss2 == eta_2body[i].ss2) {
            // 1/2 eta_abcd Gamma_cdab n_a n_b \bar{n}_c \bar{n}_d,
            // only index switch beteen a <-> b and c <-> d, else 0
            // 1/2 (eta_abcd Gamma_cdab - eta_bacd Gamma_cdba + eta_badc
            // Gamma dcba - eta_abdc Gamma_dcab) = 2 eta_abcd Gamma_cdab
            dEds += 2 * eta_2body[i].me * g.me;

            continue;
          }
        }
      }
    }
  }
}

void RungeKutta::set_one_body_flow_eq(Hamiltonian H, Generator eta) {
  int num_states = H.get_num_states();

  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();

  std::vector<OneBody> eta1 = eta.get_one_body();
  std::vector<TwoBody> eta2 = eta.get_two_body();

  int eta1_size = eta1.size();
  int eta2_size = eta2.size();

  for (int i = 0; i < eta1_size; ++i) {
    for (const OneBody &ff : f) {
      // first sum
      if (eta1[i].ss == ff.s) {
        dfds = dfds + OneBody(eta1[i].s, ff.ss, eta1[i].me * ff.me);
      }
      if (eta1[i].ss == ff.ss) {
        dfds = dfds + OneBody(eta1[i].s, ff.s, eta1[i].me * ff.me);
      }
      if (eta1[i].s == ff.s) {
        dfds = dfds + OneBody(eta1[i].ss, ff.ss, -eta1[i].me * ff.me);
      }
      if (eta1[i].s == ff.ss) {
        dfds = dfds + OneBody(eta1[i].ss, ff.s, -eta1[i].me * ff.me);
      }
    }

    // for (int i = 0; i < eta1_size; ++i) {
    //   if (eta1[i].s < 4 && eta1[i].ss > 3) {
    //     for (const Twobody &gamma : gamma) {
    //     }
    //   }
    // }
  }
}

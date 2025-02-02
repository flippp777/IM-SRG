#include "Hamiltonian.hpp"
#include <iostream>

void  Hamiltonian::set_generators() {

  int f_size = f.size();
  int gamma_size = gamma.size();

#pragma omp parallel
  {
    // one-body part
    std::vector<OneBody> local_eta_1body;

#pragma omp for
    for (int i = 0; i < f_size; ++i) {
      // f_ph
      if (f[i].s > 3 && f[i].ss < 4) {
        double delta = 0;
        for (int j = 0; j < f_size; ++j) {
          // f_pp
          if (f[j].s == f[i].s && f[j].ss == f[j].s) {
            delta += f[j].me;
            // f_hh
          } else if (f[j].s == f[i].ss && f[j].ss == f[j].s) {
            delta += -f[j].me;
          }
        }
        if (delta != 0) {
          double me = f[i].me / delta;
          if (std::abs(me) > 1e-9) {
            local_eta_1body.push_back(
                OneBody(f[i].s, f[i].ss, f[i].me / delta));
            // eta_1body.push_back(OneBody(f[i].ss, f[i].s, f[i].me /
            // delta));
          }
        } else {
          std::cout << "delta_onebody = 0" << std::endl;
        }
        // - H.c. (f_hp)
      } else if (f[i].s < 4 && f[i].ss > 3) {
        double delta = 0;
        for (int j = 0; j < f_size; ++j) {
          // f_hh
          if (f[j].s == f[i].s && f[j].ss == f[j].s) {
            delta += -f[j].me;
            // f_pp
          } else if (f[j].s == f[i].ss && f[j].ss == f[j].s) {
            delta += f[j].me;
          }
        }
        if (delta != 0) {
          double me = f[i].me / delta;
          if (std::abs(me) > 1e-9) {
            local_eta_1body.push_back(
                OneBody(f[i].s, f[i].ss, f[i].me / delta));
          }
        } else {
          std::cout << "delta_onebody = 0" << std::endl;
        }
      }
    }
#pragma omp critical
    eta_1body.insert(eta_1body.end(), local_eta_1body.begin(),
                     local_eta_1body.end());

    // two-body part
    std::vector<TwoBody> local_eta_2body;

#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      // Gamma_pp'hh'
      if (gamma[i].s1 > 3 && gamma[i].s2 > 3 && gamma[i].ss1 < 4 &&
          gamma[i].ss2 < 4) {
        double delta = 0;
        for (const OneBody &f : f) {
          // f_pp
          if (f.s == gamma[i].s1 && f.s == f.ss) {
            delta += f.me;
            // f_p'p'
          } else if (f.s == gamma[i].s2 && f.s == f.ss) {
            delta += f.me;
            // f_hh
          } else if (f.s == gamma[i].ss1 && f.s == f.ss) {
            delta += -f.me;
            // f_h'h'
          } else if (f.s == gamma[i].ss2 && f.s == f.ss) {
            delta += -f.me;
          }
        }
        if (delta != 0) {
          double me = gamma[i].me / delta;
          if (std::abs(me) > 1e-9) {
            local_eta_2body.push_back(TwoBody(gamma[i].s1, gamma[i].s2,
                                              gamma[i].ss1, gamma[i].ss2, me));
            // local_eta_2body.push_back(TwoBody(Gamma[i].ss2,
            // Gamma[i].ss1,
            //                             Gamma[i].s2, Gamma[i].s1, me));
          }
        } else {
#pragma omp critical
          std::cout << "delta_onebody = 0" << std::endl;
        }
      }
      // Gamma_hh'pp'
      if (gamma[i].s1 < 4 && gamma[i].s2 < 4 && gamma[i].ss1 > 3 &&
          gamma[i].ss2 > 3) {
        double delta = 0;
        for (const OneBody &f : f) {
          // f_hh
          if (f.s == gamma[i].s1 && f.s == f.ss) {
            delta += -f.me;
            // f_h'h'
          } else if (f.s == gamma[i].s2 && f.s == f.ss) {
            delta += -f.me;
            // f_pp
          } else if (f.s == gamma[i].ss1 && f.s == f.ss) {
            delta += f.me;
            // f_p'p'
          } else if (f.s == gamma[i].ss2 && f.s == f.ss) {
            delta += f.me;
          }
        }
        if (delta != 0) {
          double me = gamma[i].me / delta;
          if (std::abs(me) > 1e-9) {
            local_eta_2body.push_back(TwoBody(gamma[i].s1, gamma[i].s2,
                                              gamma[i].ss1, gamma[i].ss2, me));
          }
        } else {
          std::cout << "delta_onebody = 0" << std::endl;
        }
      }
    }
#pragma omp critical
    eta_2body.insert(eta_2body.end(), local_eta_2body.begin(),
                     local_eta_2body.end());
  }
}

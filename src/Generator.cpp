#include "Generator.hpp"
#include <algorithm>
#include <iostream>
#include <vector>

Generator::Generator(Hamiltonian H) {
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();
  int num_states = H.get_num_states();

#pragma omp parallel
  {
// one-body part
#pragma omp for nowait
  for (int i = 0; i < f.size(); ++i) {
    if (f[i].s < 4 && f[i].ss > 3) {
      double delta = 0;
      for (int j = 0; j < f.size(); ++j) {
        if (f[j].s == f[i].s && f[j].ss == f[i].s) {
          delta += -f[j].me;
        } else if (f[j].s == f[i].ss && f[j].ss == f[i].ss) {
          delta += f[j].me;
        }
      }
      if (delta != 0) {
        double me = f[i].me / delta;
        if (std::abs(me) > 1e-9) {
#pragma omp critical
          eta_1body.push_back(OneBody(f[i].s, f[i].ss, f[i].me / delta));
        }
      } else {
        std::cout << "delta_onebody = 0" << std::endl;
      }
    }
  }

// two-body part
#pragma omp for
  for (int i = 0; i < Gamma.size(); ++i) {
    if (Gamma[i].s1 < 4 && Gamma[i].s2 < 4 && Gamma[i].ss1 > 3 &&
        Gamma[i].ss2 > 3) {
      double delta = 0;
      for (const OneBody &f : f) {
        if (f.s == Gamma[i].s1 && f.s == f.ss) {
          delta += -f.me;
        } else if (f.s == Gamma[i].s2 && f.s == f.ss) {
          delta += -f.me;
        } else if (f.s == Gamma[i].ss1 && f.s == f.ss) {
          delta += f.me;
        } else if (f.s == Gamma[i].ss2 && f.s == f.ss) {
          delta += f.me;
        }
      }
      if (delta != 0) {
        double me = Gamma[i].me / delta;
        if (std::abs(me) > 1e-9) {
#pragma omp critical
          eta_2body.push_back(TwoBody(Gamma[i].s1, Gamma[i].s2, Gamma[i].ss1,
                                      Gamma[i].ss2, me));
        }
      } else {
        std::cout << "delta_onebody = 0" << std::endl;
      }
    }
  }
}}

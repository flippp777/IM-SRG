#include "Generator.hpp"
#include <cmath>
#include <iostream>
#include <vector>

Generator::Generator(Hamiltonian H) {
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();
  int f_size = f.size();
  int gamma_size = Gamma.size();

#pragma omp parallel
  {
// one-body part
#pragma omp for nowait
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
#pragma omp critical
            eta_1body.push_back(OneBody(f[i].s, f[i].ss, f[i].me / delta));
            eta_1body.push_back(OneBody(f[i].ss, f[i].s, f[i].me / delta));
          }
        } else {
          std::cout << "delta_onebody = 0" << std::endl;
        }
        // - H.c. (f_hp)
      } 
//       else if (f[i].s < 4 && f[i].ss > 3) {
//         double delta = 0;
//         for (int j = 0; j < f_size; ++j) {
//           // f_hh
//           if (f[j].s == f[i].s && f[j].ss == f[j].s) {
//             delta += -f[j].me;
//             // f_pp
//           } else if (f[j].s == f[i].ss && f[j].ss == f[j].s) {
//             delta += f[j].me;
//           }
//         }
//         if (delta != 0) {
//           double me = f[i].me / delta;
//           if (std::abs(me) > 1e-9) {
// #pragma omp critical
//             eta_1body.push_back(OneBody(f[i].s, f[i].ss, f[i].me / delta));
//           }
//         } else {
//           std::cout << "delta_onebody = 0" << std::endl;
//         }
//       }
    }

// two-body part
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      // Gamma_pp'hh'
      if (Gamma[i].s1 > 3 && Gamma[i].s2 > 3 && Gamma[i].ss1 < 4 &&
          Gamma[i].ss2 < 4) {
        double delta = 0;
        for (const OneBody &f : f) {
          // f_pp
          if (f.s == Gamma[i].s1 && f.s == f.ss) {
            delta += f.me;
            // f_p'p'
          } else if (f.s == Gamma[i].s2 && f.s == f.ss) {
            delta += f.me;
            // f_hh
          } else if (f.s == Gamma[i].ss1 && f.s == f.ss) {
            delta += -f.me;
            // f_h'h'
          } else if (f.s == Gamma[i].ss2 && f.s == f.ss) {
            delta += -f.me;
          }
        }
        if (delta != 0) {
          double me = Gamma[i].me / delta;
          if (std::abs(me) > 1e-9) {
#pragma omp critical
            eta_2body.push_back(TwoBody(Gamma[i].s1, Gamma[i].s2,
                                        Gamma[i].ss1, Gamma[i].ss2, me));
            eta_2body.push_back(TwoBody(Gamma[i].ss2, Gamma[i].ss1,
                                        Gamma[i].s2, Gamma[i].s1, me));
          }
        } else {
          std::cout << "delta_onebody = 0" << std::endl;
        }
      }
      // Gamma_hh'pp'
//       if (Gamma[i].s1 < 4 && Gamma[i].s2 < 4 && Gamma[i].ss1 > 3 &&
//           Gamma[i].ss2 > 3) {
//         double delta = 0;
//         for (const OneBody &f : f) {
//           // f_hh
//           if (f.s == Gamma[i].s1 && f.s == f.ss) {
//             delta += -f.me;
//             // f_h'h'
//           } else if (f.s == Gamma[i].s2 && f.s == f.ss) {
//             delta += -f.me;
//             // f_pp
//           } else if (f.s == Gamma[i].ss1 && f.s == f.ss) {
//             delta += f.me;
//             // f_p'p'
//           } else if (f.s == Gamma[i].ss2 && f.s == f.ss) {
//             delta += f.me;
//           }
//         }
//         if (delta != 0) {
//           double me = Gamma[i].me / delta;
//           if (std::abs(me) > 1e-9) {
// #pragma omp critical
//             eta_2body.push_back(TwoBody(Gamma[i].s1, Gamma[i].s2,
//                                         Gamma[i].ss1, Gamma[i].ss2, me));
//           }
//         } else {
//           std::cout << "delta_onebody = 0" << std::endl;
//         }
//       }
    }
  }
}

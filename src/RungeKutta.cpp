#include "RungeKutta.hpp"
#include "TwoBody.hpp"
#include <iostream>
#include <vector>

RungeKutta::RungeKutta(Hamiltonian H, Generator eta) {
  double E = H.get_zero_body();
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();

  std::vector<OneBody> eta_1body = eta.get_one_body();
  std::vector<TwoBody> eta_2body = eta.get_two_body();

  // E flow equation
  double dEds = 0;

  int eta_1body_size = eta_1body.size();
  int f_size = f.size();
  int min_1body = std::min(eta_1body_size, f_size);
#pragma omp parallel
  {
    if (f_size == min_1body) {
#pragma omp for
      for (int i = 1; i < f_size; ++i) {
        if (f[i].s < 4 && f[i].ss > 3) {
          std::cout << "f = "<< "\t" << f[i].s << "\t" << f[i].ss << "\t" << f[i].me
                    << std::endl;
          for (const OneBody &eta : eta_1body) {
            if (eta.s == f[i].s && eta.ss == f[i].ss) {
              std::cout << "eta = " << "\t" << eta.s << "\t" << eta.ss << "\t"
                        << eta.me << std::endl;
#pragma omp critical
              dEds += 2 * eta.me * f[i].me;
            }
          }
        }
      }
    } else {
#pragma omp for
      for (int i = 1; i < eta_1body_size; ++i) {
        if (eta_1body[i].s < 4 && eta_1body[i].ss > 3) {
          std::cout << "eta = " << "\t" << eta_1body[i].s << "\t"
                    << eta_1body[i].ss << "\t" << eta_1body[i].me << std::endl;
          for (const OneBody &ff : f) {
            if (ff.s == eta_1body[i].s && ff.ss == eta_1body[i].ss) {
              std::cout << "f = " << "\t" << ff.s << "\t" << ff.ss << "\t" << ff.me
                        << std::endl;
#pragma omp critical
              dEds += 2 * eta_1body[i].me * ff.me;
            }
          }
        }
      }
    }
  }
 flow_equations.push_back(dEds);
}

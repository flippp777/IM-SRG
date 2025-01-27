#include "RungeKutta.hpp"
#include "TwoBody.hpp"
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
  int max_1body = std::max(eta_1body_size, f_size);
  if (f_size == max_1body) {
  }
}

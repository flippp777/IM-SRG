#include "Hamiltonian.hpp"
#include <iostream>

double euler_step(double x, const double dxds, const double ds) {
  return ds * dxds;
}

void Hamiltonian::euler(const double s0, const double s, const double ds) {
  int n = (int)((s - s0) / ds);
  std::cout << "n = " << n << std::endl;

  std::unordered_map<OneBodyKey, double> f_map = vec_to_map(f);
  std::unordered_map<TwoBodyKey, double> gamma_map = vec_to_map(gamma);

  for (int i = 0; i <= n; ++i) {
    set_flow_eq();
    std::cout << "E = " << E << ",\t" << "dEds = " << dEds << ";\t"
              << "s = " << s0 + i * ds << std::endl;

    E += euler_step(E, dEds, ds);

    for (const auto &[key, me] : dfds_map) {
      if (f_map.count(key) != 0) {
        f_map[key] += euler_step(f_map[key], me, ds);
      } else {
        f_map.insert({key, euler_step(0, me, ds)});
      }
    }

    for (const auto &[key, me] : dgammads_map) {
      if (gamma_map.count(key) != 0) {
        gamma_map[key] += euler_step(gamma_map[key], me, ds);
      } else {
        gamma_map.insert({key, euler_step(0, me, ds)});
      }
    }
    f = map_to_vec(f_map);
    gamma = map_to_vec(gamma_map);
  }
}

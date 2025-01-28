#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "RungeKutta.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <ostream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  // auto start = std::chrono::high_resolution_clock::now();

  Hamiltonian H = Hamiltonian(filepath);
  Generator eta = Generator(H);

  // auto stop_constructors = std::chrono::high_resolution_clock::now();
  //
  // auto duration_constructors =
  // std::chrono::duration_cast<std::chrono::seconds>(
  //     stop_constructors - start);
  // std::cout << "Constructor H and eta time: " <<
  // duration_constructors.count()
  //           << " s" << std::endl;
  // std::cout << std::endl;

  // std::vector<TwoBody> gamma = H.get_two_body();
  // std::cout << "Gamma:" << std::endl;
  // for (const TwoBody &gg : gamma) {
  //   std::cout << gg.s1 << "\t" << gg.s2 << "\t" << gg.ss1 << "\t"
  //   <<gg.ss2 << "\t" << gg.me << std::endl;
  // }
  // std::cout << std::endl;
  //
  // std::vector<TwoBody> eta_2 = eta.get_two_body();
  // std::cout << "eta_2body:" << std::endl;
  // for (const TwoBody &eta : eta_2) {
  //   std::cout << eta.s1 << "\t" << eta.s2 << "\t" << eta.ss1 << "\t" <<
  //   eta.ss2 << "\t" << eta.me << std::endl;
  // }
  // std::cout << std::endl;

  RungeKutta rk = RungeKutta(H, eta);

  std::cout << "E0 = " << H.get_zero_body() << std::endl;
  std::cout << std::endl;

  double dEds = rk.get_zero_body_flow_eq();
  std::cout << "dEds = " << dEds << std::endl;
  std::cout << std::endl;

  std::vector<OneBody> dfds = rk.get_one_body_flow_eq();
  sort(dfds);
  std::cout << "dfds:" << std::endl;
  for (const auto &element : dfds) {
    std::cout << element.s << "\t" << element.ss << "\t" << element.me << std::endl;
  }
  

  return 0;
}

#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "RungeKutta.hpp"
#include "time.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  auto start = std::chrono::high_resolution_clock::now();

  Hamiltonian H = Hamiltonian(filepath);

  auto stop = std::chrono::high_resolution_clock::now();
  print_duration("H construnction: ", start,stop);

  // H.print_zero_body();
  // H.print_one_body();
  // H.print_two_body();


  // Generator eta = Generator(H);

  // std::vector<TwoBody> eta_2 = eta.get_two_body();
  // std::cout << "eta_2body:" << std::endl;
  // for (const TwoBody &eta : eta_2) {
  //   std::cout << eta.s1 << "\t" << eta.s2 << "\t" << eta.ss1 << "\t" <<
  //   eta.ss2 << "\t" << eta.me << std::endl;
  // }
  // std::cout << std::endl;

  // RungeKutta rk = RungeKutta(H, eta);
  //

  //
  // double dEds = rk.get_zero_body_flow_eq();
  // std::cout << "dEds = " << dEds << std::endl;
  // std::cout << std::endl;
  //
  // std::vector<OneBody> dfds = rk.get_one_body_flow_eq();
  // sort(dfds);
  // std::cout << "dfds:" << std::endl;
  // for (const auto &element : dfds) {
  //   std::cout << element.s << "\t" << element.ss << "\t" << element.me << std::endl;
  // }
  

  return 0;
}

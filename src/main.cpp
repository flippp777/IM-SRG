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

  auto start_H = std::chrono::high_resolution_clock::now();

  Hamiltonian H = Hamiltonian(filepath);

  auto stop_H = std::chrono::high_resolution_clock::now();
  print_duration("H construnction: ", start_H, stop_H);

  // H.print_zero_body();
  // H.print_one_body();
  // H.print_two_body();


  auto start_eta = std::chrono::high_resolution_clock::now();
  Generator eta = Generator(H);
  auto stop_eta = std::chrono::high_resolution_clock::now();
  print_duration("eta construnction: ", start_eta, stop_eta);


  auto start_rk = std::chrono::high_resolution_clock::now();
  RungeKutta rk = RungeKutta(H, eta);
  auto stop_rk = std::chrono::high_resolution_clock::now();
  print_duration("flow eq. construnction: ", start_rk, stop_rk);
  
  double dEds = rk.get_zero_body_flow_eq();
  std::cout << "dEds = " << dEds << std::endl;
  std::cout << std::endl;

  // std::vector<OneBody> dfds = rk.get_one_body_flow_eq();
  // sort(dfds);
  // std::cout << "dfds:" << std::endl;
  // for (const auto &element : dfds) {
  //   std::cout << element.s << "\t" << element.ss << "\t" << element.me << std::endl;
  // }
  

  return 0;
}

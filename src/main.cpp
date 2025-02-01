#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "FlowEQ.hpp"
#include "TwoBody.hpp"
#include "time.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <ostream>
#include <unordered_map>
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

  std::vector<OneBody> f = H.get_one_body();
  std::cout << "f size: " << f.size() << std::endl;

  std::vector<TwoBody> g = H.get_two_body();
  std::cout << "gamma size: " << g.size() << std::endl;


  auto start_eta = std::chrono::high_resolution_clock::now();
  Generator eta = Generator(H);
  auto stop_eta = std::chrono::high_resolution_clock::now();
  print_duration("eta construnction: ", start_eta, stop_eta);


  auto start_flow = std::chrono::high_resolution_clock::now();
  FlowEQ flow = FlowEQ(H, eta);
  auto stop_flow = std::chrono::high_resolution_clock::now();
  print_duration("flow eq. construnction: ", start_flow, stop_flow);


  // rk.print_zero_body();
  // rk.print_one_body();
  
  std::vector<OneBody> dfds = flow.get_one_body_flow_eq();
  std::cout << "dfds size: " << dfds.size() << std::endl;

  std::vector<TwoBody> dgds = flow.get_two_body_flow_eq();
  std::cout << "dgds size: " << dgds.size() << std::endl;

  

  return 0;
}

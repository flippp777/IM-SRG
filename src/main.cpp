#include "FlowEQ.hpp"
#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "time.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
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
  print_duration("H construction: ", start_H, stop_H);
  
  // H.print_zero_body();
  // H.print_one_body();
  // H.print_two_body();

  auto start_eta = std::chrono::high_resolution_clock::now();
  Generator eta = Generator(H);
  auto stop_eta = std::chrono::high_resolution_clock::now();
  print_duration("eta construction: ", start_eta, stop_eta);

  FlowEQ flow = FlowEQ(H, eta);

  // rk.print_zero_body();
  // rk.print_one_body();

  auto start_flow = std::chrono::high_resolution_clock::now();
  flow.euler(0, 10, 0.1);
  auto stop_flow = std::chrono::high_resolution_clock::now();
  print_duration("flow time: ", start_flow, stop_flow);

  return 0;
}

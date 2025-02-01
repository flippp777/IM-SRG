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

  Hamiltonian H = Hamiltonian(filepath);

  // H.print_zero_body();
  // H.print_one_body();
  // H.print_two_body();

  Generator eta = Generator(H);

  FlowEQ flow = FlowEQ(H, eta);

  // rk.print_zero_body();
  // rk.print_one_body();

  auto start_flow = std::chrono::high_resolution_clock::now();
  flow.euler(0, 0.01, 0.001);
  auto stop_flow = std::chrono::high_resolution_clock::now();
  print_duration("flow time: ", start_flow, stop_flow);

  return 0;
}

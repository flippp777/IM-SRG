#include "FlowEQ.hpp"
#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "time.hpp"

#include <chrono>

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
  H.set_generators();
  auto stop_eta = std::chrono::high_resolution_clock::now();
  print_duration("eta construction: ", start_eta, stop_eta);

  // auto start_flow_eq = std::chrono::high_resolution_clock::now();
  // FlowEQ flow = FlowEQ(H);
  // auto stop_flow_eq = std::chrono::high_resolution_clock::now();
  // print_duration("flow eq. construction: ", start_flow_eq, stop_flow_eq);

  // rk.print_zero_body();
  // rk.print_one_body();

  // auto start_flow = std::chrono::high_resolution_clock::now();
  // flow.euler(0, 0.1, 0.1);
  // auto stop_flow = std::chrono::high_resolution_clock::now();
  // print_duration("flow time: ", start_flow, stop_flow);

  return 0;
}

#include "Hamiltonian.hpp"
#include "time.hpp"

#include <chrono>
#include <iostream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);

  std::vector<int> vec;
  // std::vector<OneBody> f = H.get_one_body();
  // sort(f);
  // for (const auto &ff : f){
  //   std::cout << ff.s  << "\t" << ff.ss << "\t" << ff.me <<std::endl;
  // }
  
  // H.print_zero_body();
  // H.print_one_body();
  // H.print_two_body();

  H.set_generators();

  H.set_flow_eq();

  // auto start_flow = std::chrono::high_resolution_clock::now();
  // H.euler(0, 0, 0.1);
  // auto stop_flow = std::chrono::high_resolution_clock::now();
  // print_duration("flow time: ", start_flow, stop_flow);

  return 0;
}

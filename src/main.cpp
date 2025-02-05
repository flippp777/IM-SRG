#include "Hamiltonian.hpp"

#include <chrono>
#include <iostream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);

  // std::vector<OneBody> f = H.get_one_body();
  // sort(f);
  // for (const auto &ff : f){
  //   std::cout << ff.s  << "\t" << ff.ss << "\t" << ff.me <<std::endl;
  // }
  
  // H.print_zero_body();
  // H.print_one_body();
  // H.print_two_body();

  // H.set_generators();

  // H.set_flow_eq();

  H.euler(0.001, 1, 0.003);

  return 0;
}

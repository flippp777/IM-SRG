#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "Hamiltonian.hpp"

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);

  std::cout << "E0 = " << H.get_zero_body() << std::endl;

  std::vector<OneBody> f0 = H.get_one_body();
  for (const auto &element : f0) {
    std::cout << element.s << " " << element.ss << " " << element.me << std::endl;
  }
  

  // H.Gamma0();
  // for (int i = 0; i < H.Gamma.size(); i++) {
  //   for (int j = 0; j < H.Gamma[i].size(); j++) {
  //     std::cout << H0.Gamma[i][j][0][0] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  return 0;
}

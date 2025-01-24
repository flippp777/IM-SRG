#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"

#include <iostream>
#include <ostream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);

  std::cout << "E0 = " << H.get_zero_body() << std::endl;
  std::cout << std::endl;

  std::cout << "f" << std::endl;
  std::vector<OneBody> f0 = H.get_one_body();
  for (const auto &element : f0) {
    std::cout << element.s << " " << element.ss << " " << element.me
              << std::endl;
  }
  std::cout << std::endl;

  Generator eta = Generator(H);

  std::vector<OneBody> eta_1body = eta.get_one_body();
  std::cout << "eta_1body" << std::endl;
  for (const auto &element : eta_1body) {
    std::cout << element.s << " " << element.ss << " " << element.me
              << std::endl;
  }
  std::cout << std::endl;

  std::vector<TwoBody> eta_2body = eta.get_two_body();
  std::cout << "eta_2body" << std::endl;
  for (const auto &element : eta_2body) {
    std::cout << element.s1 << " " << element.s2 << " " << element.ss1 << " "
              << element.ss2 << " " << element.me << std::endl;
  }

  return 0;
}

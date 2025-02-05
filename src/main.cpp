// #include "Hamiltonian.hpp"

#include <array>
#include <chrono>
#include <iostream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  // Hamiltonian H = Hamiltonian(filepath);
  std::array<std::array<double,140>,140> dfds;
  dfds[1][1] = 2;
    std::cout << dfds[1][1] << std::endl;
  std::vector<std::vector<std::vector<std::vector<int>>>> arr(
        140, std::vector<std::vector<std::vector<int>>>(
            140, std::vector<std::vector<int>>(
                140, std::vector<int>(140))));
  arr[1][1][1][1] = 2;
    std::cout << arr[1][1][1][1] << std::endl;

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

  // H.euler(0.0001, 0.0060, 0.0003);
  

  return 0;
}

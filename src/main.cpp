#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "RungeKutta.hpp"
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


  auto start_eta = std::chrono::high_resolution_clock::now();
  Generator eta = Generator(H);
  auto stop_eta = std::chrono::high_resolution_clock::now();
  print_duration("eta construnction: ", start_eta, stop_eta);

  auto eta2 = eta.get_two_body();
  std::cout << eta2.size() << std::endl;

  // auto start_rk = std::chrono::high_resolution_clock::now();
  // RungeKutta rk = RungeKutta(H, eta);
  // auto stop_rk = std::chrono::high_resolution_clock::now();
  // print_duration("flow eq. construnction: ", start_rk, stop_rk);
  //
  //
  // rk.print_zero_body();
  // rk.print_one_body();
  std::unordered_map<TwoBodyKey, double> test;
  std::unordered_map<TwoBodyKey, double> test_result;
  test[TwoBodyKey(1,2,3,4)] += 1;
  test[TwoBodyKey(5,6,7,8)] += 2;
  TwoBodyKey lookup_key{1,2,3,4};
  std::cout << "hello" << std::endl;
  for (auto it = test.find(lookup_key); it != test.end(); ++it)
  {
    std::cout << "hi" << std::endl;
  
    const auto& [key, val] =*it;
    std::cout << val << std::endl;
  }

  

  return 0;
}

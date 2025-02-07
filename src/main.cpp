#include "Hamiltonian.hpp"

#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "reduce_lists.hpp"
#include <array>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <ratio>
#include <vector>

int
main ()
{
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian (filepath);

  auto start = std::chrono::high_resolution_clock::now();
  H.set_full_vector();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> dur = stop-start;
  std::cout << dur.count() << std::endl;

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
  // int status = H.flow_eq_gsl (0., NULL, NULL, NULL);
  // std::cout << status << std::endl;

  // H.euler(0.0001, 0.0060, 0.0003);

  return 0;
}

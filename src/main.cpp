#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "reduce_lists.hpp"
#include <array>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <ratio>
#include <vector>
#include <print>

int
main ()
{
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian (filepath);
  
  // H.set_num_states(5);
  
  // H.set_generators();
  // std::vector<OneBody> eta1 = H.get_one_body_generator();
  // std::ofstream eta_file ("etaOneBody_s0.0.txt");
  // for (const auto &eta : eta1) {
  //   eta_file << eta.s << "\t" << eta.ss << "\t" << eta.me << "\n";
  // }

  H.euler(0., 0.1, 0.003);


  // H.runge_kutta_gsl(0, 10);



  return 0;
}

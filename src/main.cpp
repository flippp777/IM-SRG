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
  // std::vector<TwoBody> eta2 = H.get_two_body_generator();
  // std::ofstream eta_file ("etaTwoBody_s0.0.txt");
  // for (const auto &eta : eta2) {
  //   eta_file << eta.s1 << "\t" << eta.s2 << "\t" << eta.ss1 << "\t" << eta.ss2 << "\t" << eta.me << "\n";
  // }

  // H.set_flow_eq();
  // std::unordered_map<OneBodyKey, double> dfds = H.get_one_body_flow_eq_map();
  // std::ofstream dfds_file ("dfds_s0.0.txt");
  // for (const auto &[key, me] : dfds) {
  //   dfds_file << key.s << "\t" << key.ss << "\t" << me << "\n";
  // }

  H.euler(0., 0.1, 0.003);


  // H.runge_kutta_gsl(0, 10);



  return 0;
}

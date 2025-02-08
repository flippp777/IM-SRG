#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "reduce_lists.hpp"
#include <array>
#include <chrono>
#include <cstddef>
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
  
  H.set_num_states(140);

  H.euler(0., 1, 0.02);

  // H.set_full_vector();
  // H.runge_kutta_gsl(0, 0.001);



  return 0;
}

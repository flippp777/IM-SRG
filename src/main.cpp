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

int
main ()
{
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian (filepath);
  H.set_full_vector();
  
  im_srg_ode(0.0, NULL, NULL, &H);

  double dEds = H.get_zero_body_flow_eq();
  std::cout << dEds << std::endl;


  return 0;
}

#include "im_srg_ode.hpp"
#include "Hamiltonian.hpp"
#include <gsl/gsl_errno.h>
#include <iostream>
#include <vector>

int im_srg_ode(double s, const double y[], double dyds[], void *H) {
  Hamiltonian *h = static_cast<Hamiltonian *>(H);

  // Set E, f and gamma from input vector
  h->set_elements(y);

  // Calculate generators and flow equations
  h->set_generators();
  h->set_flow_eq();

  // Get flow equations into output vector
  std::vector<double> flow_eqs = h->get_full_eqs();

  std::copy(flow_eqs.begin(), flow_eqs.end(), dyds);

  std::cout << "\b\b\b\b" << (int)(16.67*h->rk_animation_couner) 
            << "\%" << std::flush;
  h->rk_animation_couner++;

  return GSL_SUCCESS;
}

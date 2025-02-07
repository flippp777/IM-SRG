#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include <gsl/gsl_errno.h>
#include <vector>

int im_srg_ode(double s, const double y[], double dyds[], void* H) {
  Hamiltonian* h = static_cast<Hamiltonian*>(H);

  // Get E, f and gamma from the full vector

  h->set_generators();
  h->set_flow_eq();

  // fill the full vector of flow equations

  
  
  return GSL_SUCCESS;
}

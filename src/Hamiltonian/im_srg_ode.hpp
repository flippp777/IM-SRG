#ifndef ODE
#define ODE
#include "Hamiltonian.hpp"

int im_srg_ode(double s, const double y[], double dyds[], void* H);
#endif // !ODE

#include "Hamiltonian.hpp"
#include <cstddef>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <stdio.h>
#include <vector>

// void
// Hamiltonian::runge_kutta_gsl (double s0, const double s, const double ds)
// {
//   gsl_odeiv2_system sys = { Hamiltonian::flow_eq_gsl_wrapper, NULL, 12, NULL };
//
//   gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (
//       &sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
//
//   size_t f_size = f.size ();
//   size_t gamma_size = gamma.size ();
//   sort (f);
//   sort (gamma);
//   std::vector<double> y;
//   y.resize (1 + f_size + gamma_size);
//   y[0] = E;
// #pragma omp parallel for
//   for (size_t i = 0; i < f_size + gamma_size; ++i)
//     {
//       if (i >= f_size)
//         {
//           y[i + 1] = gamma[i - f_size].me;
//         }
//       else
//         {
//           y[i + 1] = f[i].me;
//         }
//     }
//
//   for (int i = 1; i <= 100; ++i)
//     {
//       double ti = i * s / 100.0;
//       int status = gsl_odeiv2_driver_apply (d, &s0, ti, y.data ());
//
//       if (status != GSL_SUCCESS)
//         {
//           std::cout << "error, retun value =" << status << std::endl;
//           break;
//         }
//     }
//   std::cout << "s = " << s0 << ", E = " << y[0] << std::endl;
// }

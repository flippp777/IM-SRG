#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <stdio.h>
#include <vector>

void
Hamiltonian::runge_kutta_gsl (double s0, const double s)
{
  size_t dimension = 1 + num_states * num_states + num_states * num_states * num_states * num_states;

  gsl_odeiv2_system sys = { im_srg_ode, NULL, dimension, this };

  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (
      &sys, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

  for (int i = 1; i <= 100; ++i)
    {
      double ti = i * s / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &s0, ti, y.data ());

      if (status != GSL_SUCCESS)
        {
          std::cout << "error, retun value =" << status << std::endl;
          break;
        }

      std::cout << "s = " << s0 << ", E = " << y[0] << std::endl;
    }

gsl_odeiv2_driver_free(d);
}

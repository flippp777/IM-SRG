#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <sys/resource.h>
#include <vector>

void
Hamiltonian::runge_kutta_gsl (double s, const double s1)
{
  std::ofstream E_out ("E_runge-kutta.txt");
  size_t dim = 1 + num_states * num_states
               + num_states * num_states * num_states * num_states;

  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc (T, dim);
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new (1e-6, 0.0);
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc (dim);

  gsl_odeiv2_system sys = { im_srg_ode, NULL, dim, this };

  std::vector<double> y = full_vector ();

  double h = 1e-6;

  while (s < s1)
    {
      int status
          = gsl_odeiv2_evolve_apply (e, c, step, &sys, &s, s1, &h, y.data ());

      if (status != GSL_SUCCESS)
        {
          std::cout << "error, retun value =" << status << std::endl;
          break;
        }

      std::cout << "s = " << s << ", E = " << y[0] << std::endl;
      E_out << s << "\t" << y[0] << "\n";

    }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (step);
}

#include "Hamiltonian.hpp"
#include "im_srg_ode.hpp"
#include <chrono>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <ostream>
#include <ratio>
#include <sstream>
#include <stdio.h>
#include <vector>

void Hamiltonian::runge_kutta_gsl(double s, const double s1, double h) {
  std::cout << "### Starting rk integration ###" << std::endl;
  std::cout << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  std::ostringstream E_path;
  E_path << "E_rk_h" << h << "_s" << s << ".txt";
  std::ofstream E_out(E_path.str());

  size_t dim = 1 + f_index.size() + gamma_index.size();

  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(T, dim);
  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-6, 0.0);
  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(dim);

  gsl_odeiv2_system sys = {im_srg_ode, NULL, dim, this};

  std::vector<double> y = full_vector();
  // if (y.size() != dim) {
  //   std::cerr << "Error: Vector size mismatch" << std::endl;
  //   return;
  // }

  while (s < s1) {
    auto start_step = std::chrono::high_resolution_clock::now();
    std::cout << "0\%" << std::flush;
    rk_animation_couner = 1;

    int status =
        gsl_odeiv2_evolve_apply(e, c, step, &sys, &s, s1, &h, y.data());

    if (status != GSL_SUCCESS) {
      std::cout << "error, retun value =" << status << std::endl;
      break;
    }

    auto stop_step = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double =
        stop_step - start_step;
    
    std::cout << "\b\b\b\b" <<"s = " << s << ", E = " << y[0] << std::endl;
    std::cout << "Step time: " << ms_double.count() * 1e-3 << "s"
              << std::endl;
    std::cout << std::endl;
    E_out << s << "\t" << y[0] << std::endl;
  }

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(step);

  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = stop - start;
  std::cout << "Integration time: " << ms_double.count() * 1e-3/60 << "m"
            << std::endl;
}

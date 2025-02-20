#include "Hamiltonian/Hamiltonian.hpp"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <ostream>
#include <string>

// Usage: Execute the IM-SRG program **in the IM-SRG directory** with flags for 
//        different options to integrate the flow equations
//        from 0 to a desired point.
//
//
// -i (string): Integration method, use "rk" for Runge-Kutta or "e" for Euler
//  defaut: rk
//
// -N (int): Number of single particle states used in the Basis (0-140)
//  default: 140
//
// -s1 (double): Integrate up to a flow parameter of s1
//  defaut: 1
//
// -ds (double): Initial step size in the case of Runge-Kutta and
//               step size in the case of Euler
//  default: 0.1
//
// 
//
// Example: ./IM-SRG -i rk -N 140 -s1 5 -ds 1e-6
//    Integrates the flow equations with Runge-Kutta and 140 single
//    particle states from s=0 to s=5 with a initial step size
//    of 10^(-6)

int main(int argc, char *argv[]) {
  std::string integration_method = "rk";
  int N = 140;
  double s1 = 1;
  double ds = 0.1;
  bool output = 0;
  int opt;

  while ((opt = getopt(argc, argv, ":i:N:s:d:og:")) != -1) {
    switch (opt) {
    case 'i':
      integration_method = optarg;
      break;
    case 'N':
      N = std::stoi(optarg);
      break;
    case 's':
      s1 = std::stod(optarg);
      break;
    case 'd':
      ds = std::stod(optarg);
      break;
    case 'o':
      output = 1;
      break;
    case '?':
      std::cout << "Wrong specifications" << std::endl;
      return 0;
    }
  }

  std::cout << "Inetgration method: " << integration_method << std::endl;
  std::cout << "Number of single particle states: " << N << std::endl;
  std::cout << "Integrating to s = " << s1 << std::endl;

  const char integral = integration_method.c_str()[0];
  switch (integral) {
  case 'r':
    std::cout << "Initial step size: " << ds << std::endl;
    break;
  case 'e':
    std::cout << "Step size: " << ds << std::endl;
    break;
  case '?':
    std::cout << "Wrong integration method specification" << std::endl;
    return 0;
  }
  std::cout << std::endl;

  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);

  switch (integral) {
  case 'r':
    H.set_num_states(N);
    H.set_index_maps();
    H.runge_kutta_gsl(0, s1, ds);
    break;
  case 'e':
    H.euler_print(0, s1, ds);
    break;
  }

  return 0;
}

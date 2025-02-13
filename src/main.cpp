#include "Hamiltonian/Hamiltonian.hpp"
#include "MatrixRepresentation/MatrixRepresentation.hpp"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <ostream>
#include <string>

int main(int argc, char *argv[]) {
  std::string integration_method = "rk";
  int N = 140;
  double s1 = 1;
  double ds = 0.1;
  bool output = 0;
  int opt;

  while ((opt = getopt(argc, argv, ":i:N:s:d:o")) != -1) {
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
  std::cout << std::endl;

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

  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);

  MatrixRepresentation Matrix = MatrixRepresentation(H);
  Matrix.export_two_body("Gamma_s0.txt");

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

  Matrix.set_one_body(H.get_one_body());
  Matrix.set_two_body(H.get_two_body());

  Matrix.export_one_body("f_s0.3.txt");
  Matrix.export_two_body("Gamma_s0.3.txt");

  return 0;
}

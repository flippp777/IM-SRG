#include "Hamiltonian/Hamiltonian.hpp"
#include "matrixElementClasses/TwoBody.hpp"
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

int
main (int argc, char *argv[])
{
  std::string integration_method = "rk";
  int N = 140;
  double s1 = 3;
  double ds = 0.1;
  bool output = 0;
  int opt;

  while ((opt = getopt (argc, argv, ":i:N:s:d:o")) != -1)
    {
      switch (opt)
        {
        case 'i':
          integration_method = optarg;
          break;
        case 'N':
          N = std::stoi (optarg);
          break;
        case 's':
          s1 = std::stod (optarg);
          break;
        case 'd':
          ds = std::stod (optarg);
          break;
        case 'o':
          output = 1;
          break;
        case '?':
          std::cout << "Wrong specifications" << std::endl;
          return 0;
        }
    }

  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian (filepath);
  H.set_num_states (140);
  H.set_index_maps ();

  const char integral = integration_method.c_str ()[0];
  switch (integral)
    {
    case 'e':
      H.euler_print (0, s1, ds);
    case 'r':
      H.set_num_states (N);
      H.runge_kutta_gsl (0, s1, ds);
    case '?':
      std::cout << "Wrong integration method specification" << std::endl;
      return 0;
    }

  return 0;
}

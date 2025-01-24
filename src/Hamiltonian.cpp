#include "Hamiltonian.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <ostream>
#include <sstream>
#include <vector>

Hamiltonian::Hamiltonian(std::string filename) {
  // Open file
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "error: could not open file" << std::endl;
    return;
  }

  // save two-body matrix elemetns in Gamma
  std::string line;
  int i = 0;
  while (std::getline(file, line)) {
    std::vector<double> row;
    std::istringstream iss(line);
    double num;
    if (i > 145) {
      TwoBody gamma;
      iss >> gamma.s1 >> gamma.s2 >> gamma.ss1 >> gamma.ss2 >> gamma.me;
      Gamma.push_back(gamma);
    }
    i++;
  }
  file.close();

  // Calculate and set zero-body
  double E0 = 0.0;
#pragma omp parallel for reduction(+ : E0)
  for (int i = 0; i < Gamma.size(); ++i) {
    if (Gamma[i].s1 < 4 && Gamma[i].s2 < 4 && Gamma[i].s1 == Gamma[i].ss1 &&
        Gamma[i].s2 == Gamma[i].ss2) {
      E0 += Gamma[i].me;
    }
  }
  E = E0;

  // calculate and set one-body
  std::vector<OneBody> f0;
#pragma omp parallel for collapse(2)
  for (int i = 0; i < num_states; ++i) {
    for (int j = 0; j < num_states; ++j) {
      OneBody fij = OneBody(i, j, 0);
      for (int k = 0; k < Gamma.size(); ++k) {
        if (Gamma[k].s1 == i && Gamma[k].ss1 == j && Gamma[k].s2 < 4 &&
            Gamma[k].ss2 == Gamma[k].s2) {
          fij.me += Gamma[k].me;
        }
      }
      if (fij.me != 0) {
        f0.push_back(fij);
      }
    }
  }
  f = f0;
}

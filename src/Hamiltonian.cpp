#include "Hamiltonian.hpp"
#include "TwoBody.hpp"

#include <algorithm>
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
    std::istringstream iss(line);
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
      std::vector<int> indices = {Gamma[i].s1, Gamma[i].s2, Gamma[i].ss1,
                                  Gamma[i].ss2};
      int permutation = 0;
      do {
        E0 += Gamma[i].me * std::pow(-1, permutation);
        ++permutation;
      } while (std::next_permutation(indices.begin(), indices.end()));
      // std::cout << permutation << std::endl;
    }
  }
  E = E0;

  // calculate and set one-body
#pragma omp parallel for collapse(2)
  for (int i = 0; i < num_states; ++i) {
    for (int j = 0; j < num_states; ++j) {
      OneBody fij = OneBody(i, j, 0);
      OneBody fji = OneBody(j, i, 0);
      for (int k = 0; k < Gamma.size(); ++k) {
        // <ia|H|ja>n_a
        if (Gamma[k].s1 == i && Gamma[k].ss1 == j && Gamma[k].s2 < 4 &&
            Gamma[k].ss2 == Gamma[k].s2) {
          fij.me += Gamma[k].me;
          fji.me += Gamma[k].me;
        }
        // <ai|H|ja>n_a
        if (Gamma[k].s2 == i && Gamma[k].ss1 == j && Gamma[k].s1 < 4 &&
            Gamma[k].ss2 == Gamma[k].s1) {
          fij.me += -Gamma[k].me;
          fji.me += -Gamma[k].me;
        }
        // <ia|H|aj>n_a
        if (Gamma[k].s1 == i && Gamma[k].ss2 == j && Gamma[k].s2 < 4 &&
            Gamma[k].ss1 == Gamma[k].s2) {
          fij.me += -Gamma[k].me;
          fji.me += -Gamma[k].me;
        }
        // <ai|H|aj>n_a
        if (Gamma[k].s2 == i && Gamma[k].ss2 == j && Gamma[k].s1 < 4 &&
            Gamma[k].ss1 == Gamma[k].s1) {
          fij.me += Gamma[k].me;
          fji.me += Gamma[k].me;
        }
      }
      if (fij.me != 0 && std::abs(fij.me) > 1e-9) {
#pragma omp critical
        {
          f.push_back(fij);
          if (j != i) {
            f.push_back(fji);
          }
        }
      }
    }
  }
  // sort like Gamma
  sort(f);

  const int gamma_size = Gamma.size();
  std::vector<TwoBody> gamma_full;
  for (int i = 0;i<gamma_size;++i) {
    std::vector<int> indices = {Gamma[i].s1, Gamma[i].s2, Gamma[i].ss1,
                                Gamma[i].ss2};
    int permutation = 0;
    do {
      gamma_full.push_back(TwoBody(indices[0], indices[1], indices[2], indices[3], Gamma[i].me * std::pow(-1, permutation)));
      ++permutation;
    } while (std::next_permutation(indices.begin(), indices.end()));
  }
  Gamma = gamma_full;
}

const void Hamiltonian::print_zero_body(){
  std::cout << "E: " << E << std::endl;
}

const void Hamiltonian::print_one_body() {
  std::cout << "f: " << std::endl;
  for (const OneBody &ff : f) {
    std::cout << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
  }
  std::cout << std::endl;
}

const void Hamiltonian::print_two_body() {
  std::cout << "Gamma:" << std::endl;
  for (const TwoBody &gg : Gamma) {
    std::cout << gg.s1 << "\t" << gg.s2 << "\t" << gg.ss1 << "\t"
    <<gg.ss2 << "\t" << gg.me << std::endl;
  }
  std::cout << std::endl;
}

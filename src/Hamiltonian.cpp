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
      OneBody fij;
      fij.s = i;
      fij.ss = j;
      fij.me = 0;
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

// void Hamiltonian::f0() {
//   f.resize(sqrt_matrix_size, std::vector<double>(sqrt_matrix_size,
//   0.0));
//
// #pragma omp parallel for collapse(2)
//   for (int i = 0; i < sqrt_matrix_size; i++) {
//     for (int j = 0; j < sqrt_matrix_size; j++) {
//       double f_ij = 0.0;
//       for (const auto &element : elements) {
//         if (element.s1 == i && element.s2 < 4 && element.ss1 == j &&
//             element.ss2 == element.s2) {
//           f_ij += element.me;
//         } else if (element.s2 == i && element.s1 < 4 && element.ss2 ==
//         j &&
//                    element.ss1 == element.s1) {
//           f_ij += element.me;
//         }
//       }
//       if (abs(f_ij) > 1e-10) {
//         f[i][j] = f_ij;
//         f[j][i] = f_ij;
//       }
//     }
//   }
// }
//
// void Hamiltonian::Gamma0() {
//   Gamma.resize(sqrt_matrix_size,
//   std::vector<std::vector<std::vector<double>>>(
//                         sqrt_matrix_size,
//                         std::vector<std::vector<double>>(
//                                  sqrt_matrix_size,
//                                  std::vector<double>(sqrt_matrix_size,
//                                  0.0))));
// #pragma omp parallel for collapse(4)
//   for (int i = 0; i < sqrt_matrix_size; i++) {
//     for (int j = 0; j < sqrt_matrix_size; j++) {
//       for (int k = 0; k < sqrt_matrix_size; k++) {
//         for (int l = 0; l < sqrt_matrix_size; l++) {
//           for (const auto &element : elements) {
//             if (element.s1 == i && element.s2 == j && element.ss1 == k
//             &&
//                 element.ss2 == l) {
//               double gamma_ijkl = element.me;
//
//               std::array<int, 4> indices = {i, j, k, l};
//               int permutation = 0;
//               if (abs(gamma_ijkl) > 1e-10) {
//                 do {
//                   Gamma[indices[0]][indices[1]][indices[2]][indices[3]]
//                   =
//                       gamma_ijkl * std::pow(-1, permutation);
//                   permutation++;
//                 } while (
//                     std::next_permutation(indices.begin(),
//                     indices.end()));
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }

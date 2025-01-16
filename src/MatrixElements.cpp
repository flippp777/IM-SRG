#include "MatrixElements.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <vector>

MatrixElements::MatrixElements(std::string filename) {
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "error: could not open file" << std::endl;
    return;
  }

  std::string line;
  int i = 0;
  while (std::getline(file, line)) {
    std::vector<double> row;
    std::istringstream iss(line);
    double num;
    if (i > 2 && i < 142) {
      while (iss >> num) {
        row.push_back(num);
      }
      states.push_back(static_cast<int>(row[1]));
    } else if (i > 145) {
      MatrixElement element;
      iss >> element.s1 >> element.s2 >> element.ss1 >> element.ss2 >>
          element.me;
      elements.push_back(element);
    }
    i++;
  }
  file.close();
}

std::vector<MatrixElement>
MatrixElements::truncateElements(int Nmax) const {
  std::vector<MatrixElement> truncated_elements;
  for (const auto &element : elements) {
    if (states[element.s1] + states[element.s2] <= Nmax &&
        states[element.ss1] + states[element.ss2] <= Nmax) {
      truncated_elements.push_back(element);
    }
  }
  return truncated_elements;
}

void MatrixElements::setTruncatedElements(int Nmax) {
  elements = truncateElements(Nmax);
}

void MatrixElements::E0() {
  double E0 = 0.0;
  for (const auto &element : elements) {
    if (element.s1 < 4 && element.s2 < 4 && element.s1 == element.ss1 &&
        element.s2 == element.ss2) {
      E0 += element.me;
    }
  }
  E = E0;
}

void MatrixElements::f0() {
  f.resize(states.size(), std::vector<double>(states.size(), 0.0));
  int steps = states.size() * states.size();
  int step_count = 0;

#pragma omp parallel for collapse(2)
  for (int i = 0; i < states.size(); i++) {
    for (int j = 0; j < states.size(); j++) {
      double f_ij = 0.0;
      for (const auto &element : elements) {
        if (element.s1 == i && element.s2 < 4 && element.ss1 == j &&
            element.ss2 == element.s2) {
          f_ij += element.me;
        } else if (element.s2 == i && element.s1 < 4 && element.ss2 == j &&
                   element.ss1 == element.s1) {
          f_ij += element.me;
        }
      }
      if (abs(f_ij) > 1e-10) {
        f[i][j] = f_ij;
        f[j][i] = f_ij;
      }
    }
  }
}

void MatrixElements::Gamma0() {
  Gamma.resize(139, std::vector<std::vector<std::vector<double>>>(
                        139, std::vector<std::vector<double>>(
                                 139, std::vector<double>(139, 0.0))));
#pragma omp parallel for collapse(4)
  for (int i = 0; i < states.size(); i++) {
    for (int j = 0; j < states.size(); j++) {
      for (int k = 0; k < states.size(); k++) {
        for (int l = 0; l < states.size(); l++) {
          for (const auto &element : elements) {
            if (element.s1 == i && element.s2 == j && element.ss1 == k &&
                element.ss2 == l) {
              double gamma_ijkl = element.me;

              std::array<int, 4> indices = {i, j, k, l};
              int permutation = 0;
              if (abs(gamma_ijkl) > 1e-10) {
                do {
                  Gamma[indices[0]][indices[1]][indices[2]][indices[3]] =
                      gamma_ijkl * std::pow(-1, permutation);
                  permutation++;
                } while (
                    std::next_permutation(indices.begin(), indices.end()));
              }
            }
          }
        }
      }
    }
  }
}

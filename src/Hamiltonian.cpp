#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"

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
      TwoBody gam;
      iss >> gam.s1 >> gam.s2 >> gam.ss1 >> gam.ss2 >> gam.me;
      gamma.push_back(gam);
    }
    i++;
  }
  file.close();

  // full gamma
  const int gamma_init_size = gamma.size();
  std::vector<TwoBody> gamma_full;
  for (const auto &g : gamma) {
    if (g.s1 == g.ss1 && g.s2 == g.ss2) {
      // <ab|H|ab>
      gamma_full.push_back(TwoBody(g.s1, g.s2, g.ss1, g.ss2, g.me));
      // <ba|H|ab>
      gamma_full.push_back(TwoBody(g.s2, g.s1, g.ss1, g.ss2, -g.me));
      // <ba|H|ba>
      gamma_full.push_back(TwoBody(g.s2, g.s1, g.ss2, g.ss1, g.me));
      // <ab|H|ba>
      gamma_full.push_back(TwoBody(g.s1, g.s2, g.ss2, g.ss1, -g.me));
    } else {
      // <ab|H|cd>
      gamma_full.push_back(TwoBody(g.s1, g.s2, g.ss1, g.ss2, g.me));
      // <ba|H|cd>
      gamma_full.push_back(TwoBody(g.s2, g.s1, g.ss1, g.ss2, -g.me));
      // <ba|H|dc>
      gamma_full.push_back(TwoBody(g.s2, g.s1, g.ss2, g.ss1, g.me));
      // <ab|H|dc>
      gamma_full.push_back(TwoBody(g.s1, g.s2, g.ss2, g.ss1, -g.me));
      // <cd|H|ab>
      gamma_full.push_back(TwoBody(g.ss1, g.ss2, g.s1, g.s2, g.me));
      // <dc|H|ab>
      gamma_full.push_back(TwoBody(g.ss2, g.ss1, g.s1, g.s2, -g.me));
      // <dc|H|ba>
      gamma_full.push_back(TwoBody(g.ss2, g.ss1, g.s2, g.s1, g.me));
      // <cd|H|ba>
      gamma_full.push_back(TwoBody(g.ss1, g.ss2, g.s2, g.s1, -g.me));
    }
  }
  gamma = gamma_full;

  std::vector<OneBody> f_new;
  const int gamma_size = gamma.size();

#pragma omp parallel
  {
#pragma omp for reduction(+ : E)
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 < 4 && gamma[i].s2 < 4 && gamma[i].s1 == gamma[i].ss1 &&
          gamma[i].s2 == gamma[i].ss2) {
        E += gamma[i].me;
      }
    }

    // calculate and set one-body
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s2 < 4 && gamma[i].ss2 < 4) {
#pragma omp critical
        f = f + OneBody(gamma[i].s1, gamma[i].ss1, gamma[i].me);
      }
    }
    // remove 0
#pragma omp for
    for (int i = 0; i < f.size(); ++i) {
      if (abs(f[i].me) > 10 - 9) {
#pragma omp critical
        f_new.push_back(f[i]);
      }
    }
  }
  f = f_new;
}

const void Hamiltonian::print_zero_body() {
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
  for (const TwoBody &gg : gamma) {
    std::cout << gg.s1 << "\t" << gg.s2 << "\t" << gg.ss1 << "\t" << gg.ss2
              << "\t" << gg.me << std::endl;
  }
  std::cout << std::endl;
}

#include "Hamiltonian.hpp"
#include "../matrixElementClasses/OneBody.hpp"
#include "../matrixElementClasses/TwoBody.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <ratio>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

Hamiltonian::Hamiltonian() : E(), f(), gamma() {}

Hamiltonian::Hamiltonian(const double E, const std::vector<OneBody> f,
                         const std::vector<TwoBody> gamma)
    : E(E), f(f), gamma(gamma) {}

Hamiltonian::Hamiltonian(std::string filename) {
  std::cout << "Constructing Hamiltonian";
  auto start = std::chrono::high_resolution_clock::now();

  // Open file
  std::ifstream file(filename);
  if (!file) {
    std::cerr << "error: could not open file" << std::endl;
    return;
  }
  std::cout << "." << std::flush;

  // save two-body matrix elemetns in gamma
  gamma.reserve(7300920);
  std::string line;
  int i = 0;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    if (i > 144) {
      TwoBody gam;
      iss >> gam.s1 >> gam.s2 >> gam.ss1 >> gam.ss2 >> gam.me;
      gamma.push_back(gam);
    }
    i++;
  }
  file.close();

  std::cout << "." << std::flush;

  // full gamma
  std::vector<TwoBody> gamma_full;
  gamma_full.reserve(7500000);
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

  std::cout << "." << std::flush;

  std::vector<OneBody> f_new;
  const int gamma_size = gamma.size();

#pragma omp parallel
  {
#pragma omp for reduction(+ : E)
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 < 4 && gamma[i].s2 < 4 &&
          gamma[i].s1 == gamma[i].ss1 && gamma[i].s2 == gamma[i].ss2) {
        E += 0.5 * gamma[i].me;
      }
    }

    // calculate and set one-body
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s2 < 4 && gamma[i].s2 == gamma[i].ss2) {
#pragma omp critical
        f = f + OneBody(gamma[i].s1, gamma[i].ss1, gamma[i].me);
      }
    }
    // remove values close to 0
#pragma omp for
    for (size_t i = 0; i < f.size(); ++i) {
      if (std::abs(f[i].me) > 1e-9) {
#pragma omp critical
        f_new.push_back(f[i]);
      }
    }
  }
  f = f_new;
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> ms_double = stop - start;
  std::cout << "\t" << " done (" << ms_double.count() << "ms)"
            << std::endl;
}

const void Hamiltonian::set_index_maps() {
  int idx = 1;
  for (int i = 0; i < num_states; i++) {
    for (int j = i; j < num_states; ++j) {
      f_index.insert({{i, j}, idx});
      index_f.insert({idx, {i, j}});
      idx++;
    }
  }
  for (int i = 0; i < num_states; ++i) {
    for (int j = i + 1; j < num_states; ++j) {
      for (int k = j; k < num_states; ++k) {
        for (int l = k + 1; l < num_states; ++l) {
          gamma_index.insert({{i, j, k, l}, idx});
          index_gamma.insert({idx, {i, j, k, l}});
          idx++;
        }
      }
    }
  }
}

const std::vector<double> Hamiltonian::full_vector() {

  std::unordered_map<OneBodyKey, double> f_map;
  for (const auto &el : f) {
    f_map[{el.s, el.ss}] = el.me;
  }

  std::unordered_map<TwoBodyKey, double> gamma_map;
  for (const auto &el : gamma) {
    gamma_map[{el.s1, el.s2, el.ss1, el.ss2}] = el.me;
  }

  std::vector<double> y(1 + f_index.size() + gamma_index.size(), 0.0);
  y[0] = E;

  int f_size = f.size();
  int gamma_size = gamma.size();

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < f_size; ++i) {
      auto idx = f_index.find({f[i].s, f[i].ss});
      if (idx != f_index.end()) {
        y[idx->second] = f[i].me;
      }
    }
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      auto idx = gamma_index.find(
          {gamma[i].s1, gamma[i].s2, gamma[i].ss1, gamma[i].ss2});
      if (idx != gamma_index.end()) {
        y[idx->second] = gamma[i].me;
      }
    }
  }
  return y;
}

const void Hamiltonian::set_elements(const double y[]) {
  f.clear();
  gamma.clear();

  f.reserve(1000);
  gamma.reserve(20000000);

  int f_index_size = f_index.size();
  int gamma_index_size = gamma_index.size();

  E = y[0];

#pragma omp parallel
  {
#pragma omp for
    for (int i = 1; i < 1 + f_index_size; ++i) {
      double me = y[i];
      if (std::abs(me) > 1e-9) {
        auto idx = index_f.find(i);
        if (idx != index_f.end()) {
          const int a = std::get<0>(idx->second);
          const int b = std::get<1>(idx->second);
          if (a != b) {
#pragma omp critical
            {
              f.push_back(OneBody(a, b, me));
              f.push_back(OneBody(b, a, me));
            }
          } else {
#pragma omp critical
            f.push_back(OneBody(a, b, me));
          }
        } else {
#pragma omp critical
          std::cout << "find(" << i << ")" << std::endl;
        }
      }
    }
#pragma omp for
    for (int i = 1 + f_index_size; i < 1 + f_index_size + gamma_index_size;
         ++i) {
      double me = y[i];
      if (std::abs(me) > 1e-9) {
        auto idx = index_gamma.find(i);
        if (idx != index_gamma.end()) {
          int a = std::get<0>(idx->second);
          int b = std::get<1>(idx->second);
          int c = std::get<2>(idx->second);
          int d = std::get<3>(idx->second);
#pragma omp critical
          {
            if (a == c && b == d) {
              // <ab|H|ab>
              gamma.push_back(TwoBody(a, b, c, d, me));
              // <ba|H|ab>
              gamma.push_back(TwoBody(b, a, c, d, -me));
              // <ba|H|ba>
              gamma.push_back(TwoBody(b, a, d, c, me));
              // <ab|H|ba>
              gamma.push_back(TwoBody(a, b, d, c, -me));
            } else {
              // <ab|H|cd>
              gamma.push_back(TwoBody(a, b, c, d, me));
              // <ba|H|cd>
              gamma.push_back(TwoBody(b, a, c, d, -me));
              // <ba|H|dc>
              gamma.push_back(TwoBody(b, a, d, c, me));
              // <ab|H|dc>
              gamma.push_back(TwoBody(a, b, d, c, -me));
              // <cd|H|ab>
              gamma.push_back(TwoBody(c, d, a, b, me));
              // <dc|H|ab>
              gamma.push_back(TwoBody(d, c, a, b, -me));
              // <dc|H|ba>
              gamma.push_back(TwoBody(d, c, b, a, me));
              // <cd|H|ba>
              gamma.push_back(TwoBody(c, d, b, a, -me));
            }
          }
        } else {
#pragma omp critical
          std::cout << "find(" << i << ")" << std::endl;
        }
      }
    }
  }
  f.shrink_to_fit();
  gamma.shrink_to_fit();
}

std::vector<double> Hamiltonian::get_full_eqs() {
  std::vector<double> dyds(1 + f_index.size() + gamma_index.size(), 0.0);
  dyds[0] = dEds;

  for (const auto &[key, me] : dfds_map) {
    auto idx = f_index.find(key);
    if (idx != f_index.end())
      dyds[idx->second] = me;
  }
  for (const auto &[key, me] : dgammads_map) {
    auto idx = gamma_index.find(key);
    if (idx != gamma_index.end())
      dyds[idx->second] = me;
  }

  return dyds;
}

const void Hamiltonian::print_zero_body() {
  std::cout << "E: " << E << std::endl;
}

const void Hamiltonian::print_one_body() {
  std::cout << "f: " << std::endl;
  for (const OneBody &ff : f) {
    std::cout << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
  }
  std::cout << "f size: " << f.size() << std::endl;
  std::cout << std::endl;
}

const void Hamiltonian::print_two_body() {
  std::cout << "gamma:" << std::endl;
  for (const TwoBody &gg : gamma) {
    std::cout << gg.s1 << "\t" << gg.s2 << "\t" << gg.ss1 << "\t" << gg.ss2
              << "\t" << gg.me << std::endl;
  }
  std::cout << "gamma size: " << gamma.size() << std::endl;
  std::cout << std::endl;
}

const void Hamiltonian::print_one_body_generator() {
  std::cout << "f: " << std::endl;
  for (const OneBody &eta1 : eta_1body) {
    std::cout << eta1.s << "\t" << eta1.ss << "\t" << eta1.me << std::endl;
  }
  std::cout << std::endl;
}

const void Hamiltonian::print_two_body_generator() {
  std::cout << "Gamma:" << std::endl;
  for (const TwoBody &eta2 : eta_2body) {
    std::cout << eta2.s1 << "\t" << eta2.s2 << "\t" << eta2.ss1 << "\t"
              << eta2.ss2 << "\t" << eta2.me << std::endl;
  }
  std::cout << std::endl;
}

const void Hamiltonian::print_one_body_flow_eq() {
  std::cout << "f: " << std::endl;
  for (const OneBody &df : dfds_vec) {
    std::cout << df.s << "\t" << df.ss << "\t" << df.me << std::endl;
  }
  std::cout << std::endl;
}

const void Hamiltonian::print_two_body_flow_eq() {
  std::cout << "Gamma:" << std::endl;
  for (const TwoBody &dg : dgammads_vec) {
    std::cout << dg.s1 << "\t" << dg.s2 << "\t" << dg.ss1 << "\t" << dg.ss2
              << "\t" << dg.me << std::endl;
  }
  std::cout << std::endl;
}

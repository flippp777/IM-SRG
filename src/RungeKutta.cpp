#include "RungeKutta.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "TwoBodyQA.hpp"
#include <iostream>
#include <unordered_map>
#include <vector>

RungeKutta::RungeKutta(Hamiltonian H, Generator eta) {
  set_zero_body_flow_eq(H, eta);
  set_one_body_flow_eq(H, eta);
  set_two_body_flow_eq(H, eta);
}

void RungeKutta::set_zero_body_flow_eq(Hamiltonian H, Generator eta) {
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> gamma = H.get_two_body();

  std::vector<OneBody> eta1 = eta.get_one_body();
  std::vector<TwoBody> eta2 = eta.get_two_body();

  // E flow equation
  dEds = 0;

  const int f_size = f.size();
  const int eta1_size = eta1.size();
  const int gamma_size = gamma.size();
  const int eta2_size = eta2.size();

#pragma omp parallel
  {
#pragma omp for reduction(+ : dEds)
    for (int i = 0; i < f_size; ++i) {
      if (f[i].s > 3 && f[i].ss < 4) {
        for (const OneBody &eta : eta1) {
          if (eta.s == f[i].ss && eta.ss == f[i].s) {
            dEds += f[i].me * eta.me;
          }
        }
      } else if (f[i].s < 4 && f[i].s > 3) {
        for (const OneBody &eta : f) {
          if (eta.s == f[i].ss && eta.ss == f[i].s) {
            dEds += -f[i].me * eta.me;
          }
        }
      }
    }
#pragma omp for reduction(+ : dEds)
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 > 3 && gamma[i].s2 > 3 && gamma[i].ss1 < 4 &&
          gamma[i].ss2 < 4) {
        for (const TwoBody &eta : eta2) {
          if (eta.s1 == gamma[i].ss1 && eta.s2 == gamma[i].ss2 &&
              eta.ss1 == gamma[i].s1 && eta.ss2 == gamma[i].s2) {

            dEds += 0.5 * gamma[i].me * eta.me;
          }
        }
      }
    }
  }
}

void RungeKutta::set_one_body_flow_eq(Hamiltonian H, Generator eta) {
  int num_states = H.get_num_states();

  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> gamma = H.get_two_body();

  std::vector<OneBody> eta1 = eta.get_one_body();
  std::vector<TwoBody> eta2 = eta.get_two_body();

  const int f_size = f.size();
  const int eta1_size = eta1.size();
  const int gamma_size = gamma.size();
  const int eta2_size = eta2.size();

  std::unordered_map<OneBodyKey, double> dfds_map;

#pragma omp parallel
  {
    // first sum
    std::unordered_map<OneBodyKey, double> local_dfds;
    local_dfds.reserve(eta1_size * f_size);
#pragma omp for collapse(2)
    for (int i = 0; i < f_size; ++i) {
      for (int j = 0; j < eta1_size; ++j) {
        // first sum
        if (eta1[j].ss == f[i].s) {
          local_dfds[OneBodyKey(eta1[j].s, f[i].ss)] += eta1[j].me * f[i].me;
          local_dfds[{f[i].ss, eta1[j].s}] += eta1[j].me * f[i].me;
        }
      }
    }
// second sum
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 > 3 && gamma[i].ss1 < 4) {
        for (int j = 0; j < eta1_size; ++j) {
          if (eta1[j].s == gamma[i].ss1 && eta1[j].ss == gamma[i].s1) {
            local_dfds[{gamma[i].s2, gamma[i].ss2}] += eta1[j].me * gamma[i].me;
          }
        }
      } else if (gamma[i].s1 < 4 && gamma[i].ss1 > 3) {
        for (int j = 0; j < eta1_size; ++j) {
          if (eta1[j].s == gamma[i].ss1 && eta1[j].ss == gamma[i].s1) {
            local_dfds[{gamma[i].s2, gamma[i].ss2}] +=
                -eta1[j].me * gamma[i].me;
          }
        }
      }
    }
#pragma omp for
    for (int i = 0; i < eta2_size; ++i) {
      if (eta2[i].s1 > 3 && eta2[i].ss1 < 4) {
        for (int j = 0; j < f_size; ++j) {
          if (f[j].s == eta2[i].ss1 && f[j].ss == eta2[i].s1) {
            local_dfds[{eta2[i].s2, eta2[i].ss2}] += f[j].me * eta2[i].me;
          }
        }
      } else if (eta2[i].s1 < 4 && eta2[i].ss1 > 3) {
        for (int j = 0; j < f_size; ++j) {
          if (f[j].s == eta2[i].ss1 && f[j].ss == eta2[i].s1) {
            local_dfds[{eta2[i].s2, eta2[i].ss2}] += -f[j].me * eta2[i].me;
          }
        }
      }
    }

// third sum
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 < 4 && gamma[i].s2 < 4 && gamma[i].ss1 > 3) {
        for (int j = 0; j < eta2_size; ++j) {
          if (eta2[j].s1 == gamma[i].ss1 && eta2[j].ss1 == gamma[i].s1 &&
              eta2[j].ss2 == gamma[i].s2) {
            local_dfds[{eta2[j].s2, gamma[i].ss2}] +=
                0.5 * eta2[j].me * gamma[i].me;
            local_dfds[{gamma[i].ss2, eta2[j].s2}] +=
                0.5 * eta2[j].me * gamma[i].me;
          }
        }
      }
      if (gamma[i].s1 > 3 && gamma[i].s2 > 3 && gamma[i].ss1 < 4) {
        for (int j = 0; j < eta2_size; ++j) {
          if (eta2[j].s1 == gamma[i].ss1 && eta2[j].ss1 == gamma[i].s1 &&
              eta2[j].ss2 == gamma[i].s2) {
            local_dfds[{eta2[j].s2, gamma[i].ss2}] +=
                -0.5 * eta2[j].me * gamma[i].me;
            local_dfds[{gamma[i].ss2, eta2[j].s2}] +=
                -0.5 * eta2[j].me * gamma[i].me;
          }
        }
      }
    }
#pragma omp critical
    for (const auto &[key, me] : local_dfds) {
      dfds_map[key] += me;
    }
  }
  dfds = map_to_vec(dfds_map);
}

void RungeKutta::set_two_body_flow_eq(Hamiltonian H, Generator eta) {
  int num_states = H.get_num_states();

  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> gamma = H.get_two_body();
  TwoBodyQA gamma_qa = TwoBodyQA(gamma);

  std::vector<OneBody> eta1 = eta.get_one_body();
  std::vector<TwoBody> eta2 = eta.get_two_body();
  TwoBodyQA eta2_qa = TwoBodyQA(eta2);

  const int f_size = f.size();
  const int eta1_size = eta1.size();
  const int gamma_size = gamma.size();
  const int eta2_size = eta2.size();

  std::unordered_map<TwoBodyKey, double> global_dgammads_map;

#pragma omp parallel
  {
    std::unordered_map<TwoBodyKey, double> local_dgammads;

// first sum
// gamma terms
#pragma omp for
    for (int i = 0; i < eta1_size; ++i) {
      for (const auto &g : gamma_qa.check_s1(eta1[i].ss)) {
        local_dgammads[{eta1[i].s, g->s2, g->ss1, g->ss2}] +=
            eta1[i].me * g->me;
        local_dgammads[{g->s2, eta1[i].s, g->ss1, g->ss2}] +=
            -eta1[i].me * g->me;
      }
      for (const auto &g : gamma_qa.check_ss1(eta1[i].s)) {
        local_dgammads[{g->s1, g->s2, eta1[i].ss, g->ss2}] +=
            -eta1[i].me * g->me;
        local_dgammads[{g->s1, g->s2, g->ss2, eta1[i].ss}] +=
            eta1[i].me * g->me;
      }
    }
// eta2 terms
#pragma omp for
    for (int i = 0; i < f_size; ++i) {
      for (const auto &e : eta2_qa.check_s1(f[i].ss)) {
        local_dgammads[{f[i].s, e->s2, e->ss1, e->ss2}] += -f[i].me * e->me;
        local_dgammads[{e->s2, f[i].s, e->ss1, e->ss2}] += f[i].me * e->me;
      }
      for (const auto &e : eta2_qa.check_ss1(f[i].s)) {
        local_dgammads[{e->s1, e->s2, f[i].ss, e->ss2}] += f[i].me * e->me;
        local_dgammads[{e->s1, e->s2, e->ss2, f[i].ss}] += -f[i].me * e->me;
      }
    }
// 2. and 3. sum
#pragma omp for
    for (int i = 0; i < eta2_size; ++i) {
      // 2. sum
      if (eta2[i].ss1 < 4 && eta2[i].ss2 < 4) {
        for (const auto &g : gamma_qa.check_s1s2(eta2[i].ss1, eta2[i].ss2)) {
          local_dgammads[{eta2[i].s1, eta2[i].s2, g->ss1, g->ss2}] +=
              -0.5 * eta2[i].me * g->me;
        }
      } else if (eta2[i].s1 < 4 && eta2[i].s2 < 4) {
        for (const auto &g : gamma_qa.check_ss1ss2(eta2[i].s1, eta2[i].s2)) {
          local_dgammads[{g->s1, g->s2, eta2[i].ss1, eta2[i].ss2}] +=
              0.5 * eta2[i].me * g->me;
        }
      } else if (eta2[i].ss1 > 3 && eta2[i].ss2 > 3) {
        for (const auto &g : gamma_qa.check_s1s2(eta2[i].ss1, eta2[i].ss2)) {
          local_dgammads[{eta2[i].s1, eta2[i].s2, g->ss1, g->ss2}] +=
              0.5 * eta2[i].me * g->me;
        }
      } else if (eta2[i].s1 > 3 && eta2[i].s2 > 3) {
        for (const auto &g : gamma_qa.check_ss1ss2(eta2[i].s1, eta2[i].s2)) {
          local_dgammads[{g->s1, g->s2, eta2[i].ss1, eta2[i].ss2}] +=
              -0.5 * eta2[i].me * g->me;
        }
        // 3. sum
      } else if (eta2[i].ss1 < 4 && eta2[i].s1 > 3) {
        for (const auto &g : gamma_qa.check_s1ss1(eta2[i].ss1, eta2[i].s1)) {
          // 1
          local_dgammads[{g->s2, eta2[i].s2, g->ss2, eta2[i].ss2}] +=
              -eta2[i].me * g->me;
          // -P_12
          local_dgammads[{eta2[i].s2, g->s2, g->ss2, eta2[i].ss2}] +=
              eta2[i].me * g->me;
          // -P_34
          local_dgammads[{g->s2, eta2[i].s2, eta2[i].ss2, g->ss2}] +=
              eta2[i].me * g->me;
          // P_12 P_34
          local_dgammads[{eta2[i].s2, g->s2, eta2[i].ss2, g->ss2}] +=
              -eta2[i].me * g->me;
        }
      } else if (eta2[i].ss1 > 3 && eta2[i].s1 < 4) {
        for (const auto &g : gamma_qa.check_s1ss1(eta2[i].ss1, eta2[i].s1)) {
          // 1
          local_dgammads[{g->s2, eta2[i].s2, g->ss2, eta2[i].ss2}] +=
              eta2[i].me * g->me;
          // -P_12
          local_dgammads[{eta2[i].s2, g->s2, g->ss2, eta2[i].ss2}] +=
              -eta2[i].me * g->me;
          // -P_34
          local_dgammads[{g->s2, eta2[i].s2, eta2[i].ss2, g->ss2}] +=
              -eta2[i].me * g->me;
          // P_12 P_34
          local_dgammads[{eta2[i].s2, g->s2, eta2[i].ss2, g->ss2}] +=
              eta2[i].me * g->me;
        }
      }
    }
    #pragma omp critical
    for (const auto & [key, me] : local_dgammads) {
      global_dgammads_map[key] += me;
    }
  }
  dGammads = map_to_vec(global_dgammads_map);
}

const void RungeKutta::print_zero_body() {
  std::cout << "dEds: " << dEds << std::endl;
}

const void RungeKutta::print_one_body() {
  std::cout << "dfds: " << std::endl;
  for (const OneBody &df : dfds) {
    std::cout << df.s << "\t" << df.ss << "\t" << df.me << std::endl;
  }
  std::cout << std::endl;
}

const void RungeKutta::print_two_body() {
  std::cout << "Gamma:" << std::endl;
  for (const TwoBody &dg : dGammads) {
    std::cout << dg.s1 << "\t" << dg.s2 << "\t" << dg.ss1 << "\t" << dg.ss2
              << "\t" << dg.me << std::endl;
  }
  std::cout << std::endl;
}

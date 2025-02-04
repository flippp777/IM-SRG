#include "FlowEQ.hpp"
#include "Integrate.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "TwoBodyQA.hpp"
#include "omp.h"
#include <iostream>
#include <unordered_map>
#include <vector>

FlowEQ::FlowEQ(const Hamiltonian H_in, const Generator eta_in) {
  H = H_in;
  eta = eta_in;
}

void FlowEQ::set_flow_eq() {
  const std::vector<OneBody> f = H.get_one_body();
  const std::vector<TwoBody> gamma = H.get_two_body();
  TwoBodyQA gamma_qa = TwoBodyQA(gamma);

  const std::vector<OneBody> eta_1body = eta.get_one_body();
  const std::vector<TwoBody> eta_2body = eta.get_two_body();
  TwoBodyQA eta2_qa = TwoBodyQA(eta_2body);

  const int f_size = f.size();
  const int eta1_size = eta_1body.size();
  const int gamma_size = gamma.size();
  const int eta2_size = eta_2body.size();

  std::unordered_map<OneBodyKey, double> dfds_map;
  std::unordered_map<TwoBodyKey, double> dgammads_map;
  std::vector<std::unordered_map<TwoBodyKey, double>> thread_dgammads(
      omp_get_max_threads());

#pragma omp parallel
  {
//////////////////////////////// dEds
#pragma omp for reduction(+ : dEds)
    for (int i = 0; i < f_size; ++i) {
      if (f[i].s > 3 && f[i].ss < 4) {
        for (const OneBody &eta : eta_1body) {
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
        for (const TwoBody &eta : eta_2body) {
          if (eta.s1 == gamma[i].ss1 && eta.s2 == gamma[i].ss2 &&
              eta.ss1 == gamma[i].s1 && eta.ss2 == gamma[i].s2) {

            dEds += 0.5 * gamma[i].me * eta.me;
          }
        }
      }
    }
    ///////////////////////////////////////////////// dfds
    // first sum
    std::unordered_map<OneBodyKey, double> local_dfds;
    local_dfds.reserve(eta1_size * f_size);
#pragma omp for collapse(2)
    for (int i = 0; i < f_size; ++i) {
      for (int j = 0; j < eta1_size; ++j) {
        // first sum
        if (eta_1body[j].ss == f[i].s) {
          local_dfds[OneBodyKey(eta_1body[j].s, f[i].ss)] +=
              eta_1body[j].me * f[i].me;
          local_dfds[{f[i].ss, eta_1body[j].s}] += eta_1body[j].me * f[i].me;
        }
      }
    }
// second sum
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 > 3 && gamma[i].ss1 < 4) {
        for (int j = 0; j < eta1_size; ++j) {
          if (eta_1body[j].s == gamma[i].ss1 &&
              eta_1body[j].ss == gamma[i].s1) {
            local_dfds[{gamma[i].s2, gamma[i].ss2}] +=
                eta_1body[j].me * gamma[i].me;
          }
        }
      } else if (gamma[i].s1 < 4 && gamma[i].ss1 > 3) {
        for (int j = 0; j < eta1_size; ++j) {
          if (eta_1body[j].s == gamma[i].ss1 &&
              eta_1body[j].ss == gamma[i].s1) {
            local_dfds[{gamma[i].s2, gamma[i].ss2}] +=
                -eta_1body[j].me * gamma[i].me;
          }
        }
      }
    }
#pragma omp for
    for (int i = 0; i < eta2_size; ++i) {
      if (eta_2body[i].s1 > 3 && eta_2body[i].ss1 < 4) {
        for (int j = 0; j < f_size; ++j) {
          if (f[j].s == eta_2body[i].ss1 && f[j].ss == eta_2body[i].s1) {
            local_dfds[{eta_2body[i].s2, eta_2body[i].ss2}] +=
                f[j].me * eta_2body[i].me;
          }
        }
      } else if (eta_2body[i].s1 < 4 && eta_2body[i].ss1 > 3) {
        for (int j = 0; j < f_size; ++j) {
          if (f[j].s == eta_2body[i].ss1 && f[j].ss == eta_2body[i].s1) {
            local_dfds[{eta_2body[i].s2, eta_2body[i].ss2}] +=
                -f[j].me * eta_2body[i].me;
          }
        }
      }
    }

// third sum
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 < 4 && gamma[i].s2 < 4 && gamma[i].ss1 > 3) {
        for (int j = 0; j < eta2_size; ++j) {
          if (eta_2body[j].s1 == gamma[i].ss1 &&
              eta_2body[j].ss1 == gamma[i].s1 &&
              eta_2body[j].ss2 == gamma[i].s2) {
            local_dfds[{eta_2body[j].s2, gamma[i].ss2}] +=
                0.5 * eta_2body[j].me * gamma[i].me;
            local_dfds[{gamma[i].ss2, eta_2body[j].s2}] +=
                0.5 * eta_2body[j].me * gamma[i].me;
          }
        }
      }
      if (gamma[i].s1 > 3 && gamma[i].s2 > 3 && gamma[i].ss1 < 4) {
        for (int j = 0; j < eta2_size; ++j) {
          if (eta_2body[j].s1 == gamma[i].ss1 &&
              eta_2body[j].ss1 == gamma[i].s1 &&
              eta_2body[j].ss2 == gamma[i].s2) {
            local_dfds[{eta_2body[j].s2, gamma[i].ss2}] +=
                -0.5 * eta_2body[j].me * gamma[i].me;
            local_dfds[{gamma[i].ss2, eta_2body[j].s2}] +=
                -0.5 * eta_2body[j].me * gamma[i].me;
          }
        }
      }
    }
    for (const auto &[key, me] : local_dfds) {
#pragma omp critical
      dfds_map[key] += me;
    }
    ///////////////////////////// dGammads
    int tid = omp_get_thread_num();
    std::unordered_map<TwoBodyKey, double> &local_dgammads =
        thread_dgammads[tid];

// first sum
// gamma terms
#pragma omp for
    for (int i = 0; i < eta1_size; ++i) {
      for (const auto &g : gamma_qa.check_s1(eta_1body[i].ss)) {
        local_dgammads[{eta_1body[i].s, g->s2, g->ss1, g->ss2}] +=
            eta_1body[i].me * g->me;
        local_dgammads[{g->s2, eta_1body[i].s, g->ss1, g->ss2}] +=
            -eta_1body[i].me * g->me;
      }
      for (const auto &g : gamma_qa.check_ss1(eta_1body[i].s)) {
        local_dgammads[{g->s1, g->s2, eta_1body[i].ss, g->ss2}] +=
            -eta_1body[i].me * g->me;
        local_dgammads[{g->s1, g->s2, g->ss2, eta_1body[i].ss}] +=
            eta_1body[i].me * g->me;
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
      if (eta_2body[i].ss1 < 4 && eta_2body[i].ss2 < 4) {
        for (const auto &g :
             gamma_qa.check_s1s2(eta_2body[i].ss1, eta_2body[i].ss2)) {
          local_dgammads[{eta_2body[i].s1, eta_2body[i].s2, g->ss1, g->ss2}] +=
              -0.5 * eta_2body[i].me * g->me;
        }
      } else if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4) {
        for (const auto &g :
             gamma_qa.check_ss1ss2(eta_2body[i].s1, eta_2body[i].s2)) {
          local_dgammads[{g->s1, g->s2, eta_2body[i].ss1, eta_2body[i].ss2}] +=
              0.5 * eta_2body[i].me * g->me;
        }
      } else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3) {
        for (const auto &g :
             gamma_qa.check_s1s2(eta_2body[i].ss1, eta_2body[i].ss2)) {
          local_dgammads[{eta_2body[i].s1, eta_2body[i].s2, g->ss1, g->ss2}] +=
              0.5 * eta_2body[i].me * g->me;
        }
      } else if (eta_2body[i].s1 > 3 && eta_2body[i].s2 > 3) {
        for (const auto &g :
             gamma_qa.check_ss1ss2(eta_2body[i].s1, eta_2body[i].s2)) {
          local_dgammads[{g->s1, g->s2, eta_2body[i].ss1, eta_2body[i].ss2}] +=
              -0.5 * eta_2body[i].me * g->me;
        }
        // 3. sum
      } else if (eta_2body[i].ss1 < 4 && eta_2body[i].s1 > 3) {
        for (const auto &g :
             gamma_qa.check_s1ss1(eta_2body[i].ss1, eta_2body[i].s1)) {
          // 1
          local_dgammads[{g->s2, eta_2body[i].s2, g->ss2, eta_2body[i].ss2}] +=
              -eta_2body[i].me * g->me;
          // -P_12
          local_dgammads[{eta_2body[i].s2, g->s2, g->ss2, eta_2body[i].ss2}] +=
              eta_2body[i].me * g->me;
          // -P_34
          local_dgammads[{g->s2, eta_2body[i].s2, eta_2body[i].ss2, g->ss2}] +=
              eta_2body[i].me * g->me;
          // P_12 P_34
          local_dgammads[{eta_2body[i].s2, g->s2, eta_2body[i].ss2, g->ss2}] +=
              -eta_2body[i].me * g->me;
        }
      } else if (eta_2body[i].ss1 > 3 && eta_2body[i].s1 < 4) {
        for (const auto &g :
             gamma_qa.check_s1ss1(eta_2body[i].ss1, eta_2body[i].s1)) {
          // 1
          local_dgammads[{g->s2, eta_2body[i].s2, g->ss2, eta_2body[i].ss2}] +=
              eta_2body[i].me * g->me;
          // -P_12
          local_dgammads[{eta_2body[i].s2, g->s2, g->ss2, eta_2body[i].ss2}] +=
              -eta_2body[i].me * g->me;
          // -P_34
          local_dgammads[{g->s2, eta_2body[i].s2, eta_2body[i].ss2, g->ss2}] +=
              -eta_2body[i].me * g->me;
          // P_12 P_34
          local_dgammads[{eta_2body[i].s2, g->s2, eta_2body[i].ss2, g->ss2}] +=
              eta_2body[i].me * g->me;
        }
      }
    }
    //     for (const auto &[key, me] : local_dgammads) {
    // #pragma omp critical
    //       global_dgammads_map[key] += me;
    //     }
  }
  for (const auto &local_map : thread_dgammads) {
    for (const auto &[key, me] : local_map) {
      dgammads_map[key] += me;
    }
  }
  dfds = map_to_vec(dfds_map);
  dgammads = map_to_vec(dgammads_map);
}

/////////////////////////// Evolution
void FlowEQ::euler(const double s0, const double s, const double ds) {
  int n = (int)((s - s0) / ds);
  std::cout << "n = " << n << std::endl;

  std::unordered_map<OneBodyKey, double> f, dfds_map;
  std::unordered_map<TwoBodyKey, double> gamma, dgds_map;

  for (int i = 0; i <= n; ++i) {
    set_flow_eq();
    double E = H.get_zero_body();
    std::cout << "E = " << E << ",\t" << "dEds = " << dEds << ";\t" << "s = " << s0 + i * ds << std::endl;
#pragma omp parallel 
    {
#pragma omp sections
      {
#pragma omp section
        f = vec_to_map(H.get_one_body());
#pragma omp section
        gamma = vec_to_map(H.get_two_body());
#pragma omp section
        dfds_map = vec_to_map(dfds);
#pragma omp section
        dgds_map = vec_to_map(dgammads);
      }

#pragma omp barrier

#pragma omp single
      {
#pragma omp task
        E += euler_step(E, dEds, ds);
        for (const auto &[key, me] : dfds_map) {
#pragma omp task
          f[key] += euler_step(f[key], me, ds);
        }
        for (const auto &[key, me] : dgds_map) {
#pragma omp task
          gamma[key] += euler_step(gamma[key], me, ds);
        }
      }
#pragma omp barrier

#pragma omp single
      H.set_zero_body(E);
#pragma omp sections
      {
#pragma omp section
        H.set_one_body(map_to_vec(f));
#pragma omp section
        H.set_two_body(map_to_vec(gamma));
      }
#pragma omp single
      eta = Generator(H);
    }
    dEds = 0;
    dfds.clear();
    dgammads.clear();
  }
}

const void FlowEQ::print_zero_body() {
  std::cout << "dEds: " << dEds << std::endl;
}

const void FlowEQ::print_one_body() {
  std::cout << "dfds: " << std::endl;
  for (const OneBody &df : dfds) {
    std::cout << df.s << "\t" << df.ss << "\t" << df.me << std::endl;
  }
  std::cout << std::endl;
}

const void FlowEQ::print_two_body() {
  std::cout << "dGammads:" << std::endl;
  for (const TwoBody &dg : dgammads) {
    std::cout << dg.s1 << "\t" << dg.s2 << "\t" << dg.ss1 << "\t" << dg.ss2
              << "\t" << dg.me << std::endl;
  }
  std::cout << std::endl;
}

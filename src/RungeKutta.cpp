#include "RungeKutta.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include <iostream>
#include <vector>

RungeKutta::RungeKutta(Hamiltonian H, Generator eta) {
  set_zero_body_flow_eq(H, eta);
  set_one_body_flow_eq(H, eta);
  set_two_body_flow_eq(H, eta);
}

void RungeKutta::set_zero_body_flow_eq(Hamiltonian H, Generator eta) {
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();

  std::vector<OneBody> eta_1body = eta.get_one_body();
  std::vector<TwoBody> eta_2body = eta.get_two_body();

  // E flow equation
  dEds = 0;

  const int f_size = f.size();
  const int eta_1body_size = eta_1body.size();
  const int gamma_size = Gamma.size();
  const int eta_2body_size = eta_2body.size();

#pragma omp parallel
  {
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
      if (Gamma[i].s1 > 3 && Gamma[i].s2 > 3 && Gamma[i].ss1 < 4 &&
          Gamma[i].ss2 < 4) {
        for (const TwoBody &eta : eta_2body) {
          if (eta.s1 == Gamma[i].ss1 && eta.s2 == Gamma[i].ss2 &&
              eta.ss1 == Gamma[i].s1 && eta.ss2 == Gamma[i].s2) {

            dEds += 0.5 * Gamma[i].me * eta.me;
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

#pragma omp parallel
  {
    // first sum
    std::vector<OneBody> local_dfds;
    local_dfds.reserve(eta1_size * f_size);
#pragma omp for collapse(2)
    for (int i = 0; i < f_size; ++i) {
      for (int j = 0; j < eta1_size; ++j) {
        // first sum
        if (eta1[j].ss == f[i].s) {
          local_dfds =
              local_dfds + OneBody(eta1[j].s, f[i].ss, eta1[j].me * f[i].me);
          local_dfds =
              local_dfds + OneBody(f[i].ss, eta1[j].s, eta1[j].me * f[i].me);
        }
      }
    }

// second sum
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      if (gamma[i].s1 > 3 && gamma[i].ss1 < 4) {
        for (int j = 0; j < eta1_size; ++j) {
          if (eta1[j].s == gamma[i].ss1 && eta1[j].ss == gamma[i].s1) {
            local_dfds = local_dfds + OneBody(gamma[i].s2, gamma[i].ss2,
                                              eta1[j].me * gamma[i].me);
          }
        }
      } else if (gamma[i].s1 < 4 && gamma[i].ss1 > 3) {
        for (int j = 0; j < eta1_size; ++j) {
          if (eta1[j].s == gamma[i].ss1 && eta1[j].ss == gamma[i].s1) {
            local_dfds = local_dfds + OneBody(gamma[i].s2, gamma[i].ss2,
                                              -eta1[j].me * gamma[i].me);
          }
        }
      }
    }
#pragma omp for
    for (int i = 0; i < eta2_size; ++i) {
      if (eta2[i].s1 > 3 && eta2[i].ss1 < 4) {
        for (int j = 0; j < f_size; ++j) {
          if (f[j].s == eta2[i].ss1 && f[j].ss == eta2[i].s1) {
            local_dfds = local_dfds +
                         OneBody(eta2[i].s2, eta2[i].ss2, f[j].me * eta2[i].me);
          }
        }
      } else if (eta2[i].s1 < 4 && eta2[i].ss1 > 3) {
        for (int j = 0; j < f_size; ++j) {
          if (f[j].s == eta2[i].ss1 && f[j].ss == eta2[i].s1) {
            local_dfds = local_dfds + OneBody(eta2[i].s2, eta2[i].ss2,
                                              -f[j].me * eta2[i].me);
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
            local_dfds = local_dfds + OneBody(eta2[j].s2, gamma[i].ss2,
                                              0.5 * eta2[j].me * gamma[i].me);
            local_dfds = local_dfds + OneBody(gamma[i].ss2, eta2[j].s2,
                                              0.5 * eta2[j].me * gamma[i].me);
          }
        }
      }
      if (gamma[i].s1 > 3 && gamma[i].s2 > 3 && gamma[i].ss1 < 4) {
        for (int j = 0; j < eta2_size; ++j) {
          if (eta2[j].s1 == gamma[i].ss1 && eta2[j].ss1 == gamma[i].s1 &&
              eta2[j].ss2 == gamma[i].s2) {
            local_dfds = local_dfds + OneBody(eta2[j].s2, gamma[i].ss2,
                                              -0.5 * eta2[j].me * gamma[i].me);
            local_dfds = local_dfds + OneBody(gamma[i].ss2, eta2[j].s2,
                                              -0.5 * eta2[j].me * gamma[i].me);
          }
        }
      }
    }
#pragma omp critical
    dfds.insert(dfds.end(), local_dfds.begin(), local_dfds.end());
  }
}

void RungeKutta::set_two_body_flow_eq(Hamiltonian H, Generator eta) {
  int num_states = H.get_num_states();

  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> gamma = H.get_two_body();

  std::vector<OneBody> eta1 = eta.get_one_body();
  std::vector<TwoBody> eta2 = eta.get_two_body();

  const int f_size = f.size();
  const int eta1_size = eta1.size();
  const int gamma_size = gamma.size();
  const int eta2_size = eta2.size();

#pragma omp parallel
  {
    std::vector<TwoBody> local_dgammads;
    local_dgammads.reserve(gamma_size);
    // sum with Gamma_abcd terms
#pragma omp for
    for (int i = 0; i < gamma_size; ++i) {
      // terms in first sum with eta * Gamma
      for (int j = 0; j < eta1_size; ++j) {
        if (eta1[j].ss == gamma[i].s1) {
          local_dgammads =
              local_dgammads + TwoBody(eta1[j].s, gamma[i].s2, gamma[i].ss1,
                                       gamma[i].ss2, eta1[j].me * gamma[i].me);
          local_dgammads =
              local_dgammads + TwoBody(gamma[i].s2, eta1[j].s, gamma[i].ss1,
                                       gamma[i].ss2, -eta1[j].me * gamma[i].me);
        } else if (eta1[j].s == gamma[i].ss1) {
          local_dgammads =
              local_dgammads + TwoBody(gamma[i].s1, gamma[i].s2, eta1[j].ss,
                                       gamma[i].ss2, -eta1[j].me * gamma[i].me);
          local_dgammads =
              local_dgammads + TwoBody(gamma[i].s1, gamma[i].s2, gamma[i].ss2,
                                       eta1[j].ss, eta1[j].me * gamma[i].me);
        }
      }
    }
    #pragma omp critical
    dGammads.insert(dGammads.end(), local_dgammads.begin(), local_dgammads.end());
  }
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

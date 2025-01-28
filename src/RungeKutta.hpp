#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"

#include <vector>

#ifndef RUNGEKUTTA
#define RUNGEKUTTA

class RungeKutta {
private:
  double dEds;
  std::vector<OneBody> dfds;
  std::vector<TwoBody> dGammads;

public:
  RungeKutta(Hamiltonian H, Generator eta);

  void set_zero_body_flow_eq(Hamiltonian H, Generator eta);
  void set_one_body_flow_eq(Hamiltonian H, Generator eta);
  void set_two_body_flow_eq(Hamiltonian H, Generator eta);
  
  void evolve(double s0, double s_max, double ds);

  const double &get_zero_body_flow_eq() const {
    return dEds;
  }
  const std::vector<OneBody> &get_one_body_flow_eq() const {
    return dfds;
  }
  const std::vector<TwoBody> &get_two_body_flow_eq() const {
    return dGammads;
  }
};

#endif // !RUNGEKUTTA

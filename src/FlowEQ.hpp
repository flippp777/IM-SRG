#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"

#include <vector>

#ifndef FLOWEQ
#define FLOWEQ

class FlowEQ {
private:
  Hamiltonian H;
  Generator eta;

  double dEds;
  std::vector<OneBody> dfds;
  std::vector<TwoBody> dgammads;

public:
  FlowEQ(Hamiltonian H, Generator eta);

  void set_flow_eq();
  
  void euler(const double s0, const double s, const double ds);

  void runge_kutta(double s0, double s, double ds);

  const Hamiltonian &get_hamiltonian() const {
    return H;
  }

  const double &get_zero_body_flow_eq() const {
    return dEds;
  }
  const std::vector<OneBody> &get_one_body_flow_eq() const {
    return dfds;
  }
  const std::vector<TwoBody> &get_two_body_flow_eq() const {
    return dgammads;
  }

  const void print_zero_body();
  const void print_one_body();
  const void print_two_body();
};

#endif // !FLOWEQ

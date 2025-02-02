#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"

#include <unordered_map>
#include <vector>

#ifndef FLOWEQ
#define FLOWEQ

class FlowEQ {
private:
  Hamiltonian H;
  Generator eta;

  double dEds;
  std::vector<OneBody> dfds_vec;
  std::vector<TwoBody> dgammads_vec;

  std::unordered_map<OneBodyKey, double> dfds_map;
  std::unordered_map<TwoBodyKey, double> dgammads_map;

public:
  FlowEQ(Hamiltonian H);

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
    return dfds_vec;
  }
  const std::vector<TwoBody> &get_two_body_flow_eq() const {
    return dgammads_vec;
  }

  const void print_zero_body();
  const void print_one_body();
  const void print_two_body();
};

#endif // !FLOWEQ

#include "OneBody.hpp"
#include "TwoBody.hpp"

#include <string>
#include <vector>

#ifndef HAMILTONIAN
#define HAMILTONIAN

class Hamiltonian {
private:
  int num_states = 140;

  double E;
  std::vector<OneBody> f;
  std::vector<TwoBody> gamma;

  std::vector<OneBody> eta_1body;
  std::vector<TwoBody> eta_2body;

  double dEds;
  std::vector<OneBody> dfds_vec;
  std::vector<TwoBody> dgammads_vec;
  std::unordered_map<OneBodyKey, double> dfds_map;
  std::unordered_map<TwoBodyKey, double> dgammads_map;

public:
  Hamiltonian();
  Hamiltonian(const double E, const std::vector<OneBody> f, const std::vector<TwoBody> gamma);
  Hamiltonian(std::string filepath);

  void set_generators();
  void set_flow_eq();

  void euler(const double s0, const double s, const double ds);

  const std::vector<TwoBody> &get_two_body() const {
    return gamma;
  }
  const void set_two_body(std::vector<TwoBody> new_Gamma) {
    gamma = new_Gamma;
  }
  const std::vector<OneBody> &get_one_body() const {
    return f;
  }
  const void set_one_body(std::vector<OneBody> new_f) {
    f = new_f;
  }
  const double &get_zero_body() const {
    return E;
  }
  const void set_zero_body(double new_E) {
    E = new_E;
  }
  const int &get_num_states() const {
    return num_states;
  }

  const void print_zero_body();
  const void print_one_body();
  const void print_two_body();

};
#endif // !HAMILTONIAN

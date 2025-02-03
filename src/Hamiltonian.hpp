#include "OneBody.hpp"
#include "TwoBody.hpp"

#include <string>
#include <unordered_map>
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

  const void set_zero_body(double new_E) {
    E = new_E;
  }
  const double &get_zero_body() const {
    return E;
  }
  const void set_one_body(std::vector<OneBody> new_f) {
    f = new_f;
  }
  const std::vector<OneBody> &get_one_body() const {
    return f;
  }
  const void set_two_body(std::vector<TwoBody> new_Gamma) {
    gamma = new_Gamma;
  }
  const std::vector<TwoBody> &get_two_body() const {
    return gamma;
  }
  const int &get_num_states() const {
    return num_states;
  }

  const void set_one_body_generator(std::vector<OneBody> new_eta) {
    eta_1body = new_eta;
  }
  const std::vector<OneBody> &get_one_body_generator() const {
    return eta_1body;
  }
  const void set_two_body_generator(std::vector<TwoBody> new_eta) {
    eta_2body = new_eta;
  }
  const std::vector<TwoBody> &get_two_body_generator() const {
    return eta_2body;
  }

  const void set_zero_body_flow_eq(double new_dEds) {
    dEds = new_dEds;
  }
  const double &get_zero_body_flow_eq() const {
    return dEds;
  }
  const void set_one_body_flow_eq_vec(std::vector<OneBody> new_dfds_vec) {
    dfds_vec = new_dfds_vec;
  }
  const std::vector<OneBody> &get_one_body_flow_eq_vec() const {
    return dfds_vec;
  }
  const void set_one_body_flow_eq_map(std::unordered_map<OneBodyKey, double> new_dfds_map) {
    dfds_map = new_dfds_map;
  }
  const std::unordered_map<OneBodyKey, double> &get_one_body_flow_eq_map() const {
    return dfds_map;
  }
  const void set_two_body_flow_eq_vec(std::vector<TwoBody> new_dgammads_vec) {
    dgammads_vec = new_dgammads_vec;
  }
  const std::vector<TwoBody> &get_two_body_flow_eq_vec() const {
    return dgammads_vec;
  }
  const void set_two_body_flow_eq_map(std::unordered_map<TwoBodyKey, double> new_dgammads_map) {
    dgammads_map = new_dgammads_map;
  }
  const std::unordered_map<TwoBodyKey, double> &get_two_body_flow_eq_map() const {
    return dgammads_map;
  }

  const void print_zero_body();
  const void print_one_body();
  const void print_two_body();

};
#endif // !HAMILTONIAN

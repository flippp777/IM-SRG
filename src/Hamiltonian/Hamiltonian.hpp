#include "../matrixElementClasses/OneBody.hpp"
#include "../matrixElementClasses/TwoBody.hpp"

#include <map>
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

  std::unordered_map<OneBodyKey, int> f_index;
  std::unordered_map<TwoBodyKey, int> gamma_index;
  std::unordered_map<int, std::tuple<int, int>> index_f;
  std::unordered_map<int, std::tuple<int, int, int, int>> index_gamma;

  std::vector<OneBody> eta_1body;
  std::vector<TwoBody> eta_2body;

  double dEds;
  std::vector<OneBody> dfds_vec;
  std::vector<TwoBody> dgammads_vec;
  std::unordered_map<OneBodyKey, double> dfds_map;
  std::unordered_map<TwoBodyKey, double> dgammads_map;

public:
  int rk_animation_couner;

  Hamiltonian();
  Hamiltonian(const double E, const std::vector<OneBody> f,
              const std::vector<TwoBody> gamma);
  Hamiltonian(std::string filepath);

  const void set_num_states(int new_num) { num_states = new_num; }

  const void set_index_maps();

  const std::vector<double> get_full_vector();
  const void set_generators();
void const set_it_generators();
  const void set_flow_eq();
  const void set_flow_eq_print();
  int alt_flow_eq();

  void euler(const double s0, const double s, const double ds);
  void euler_print(const double s0, const double s, const double ds);

  const void set_elements(const double y[]);

  std::vector<double> get_full_eqs();

  void runge_kutta_gsl(double s, const double s1, const double h);

  const void set_zero_body(double new_E) { E = new_E; }
  const double &get_zero_body() { return E; }

  const void set_one_body(std::vector<OneBody> new_f) { f = new_f; }
  const std::vector<OneBody> &get_one_body() { return f; }

  const void set_two_body(std::vector<TwoBody> new_Gamma) {
    gamma = new_Gamma;
  }
  const std::vector<TwoBody> &get_two_body() { return gamma; }

  const int &get_num_states() { return num_states; }

  const void set_one_body_generator(std::vector<OneBody> new_eta) {
    eta_1body = new_eta;
  }
  const std::vector<OneBody> &get_one_body_generator() {
    return eta_1body;
  }

  const void set_two_body_generator(std::vector<TwoBody> new_eta) {
    eta_2body = new_eta;
  }
  const std::vector<TwoBody> &get_two_body_generator() {
    return eta_2body;
  }

  const void set_zero_body_flow_eq(double new_dEds) { dEds = new_dEds; }
  const double &get_zero_body_flow_eq() { return dEds; }

  const void set_one_body_flow_eq_vec(std::vector<OneBody> new_dfds_vec) {
    dfds_vec = new_dfds_vec;
  }

  const std::vector<OneBody> &get_one_body_flow_eq_vec() {
    return dfds_vec;
  }

  const void set_one_body_flow_eq_map(
      std::unordered_map<OneBodyKey, double> new_dfds_map) {
    dfds_map = new_dfds_map;
  }

  const std::unordered_map<OneBodyKey, double> &
  get_one_body_flow_eq_map() {
    return dfds_map;
  }

  const void
  set_two_body_flow_eq_vec(std::vector<TwoBody> new_dgammads_vec) {
    dgammads_vec = new_dgammads_vec;
  }

  const std::vector<TwoBody> &get_two_body_flow_eq_vec() const {
    return dgammads_vec;
  }

  const void set_two_body_flow_eq_map(
      std::unordered_map<TwoBodyKey, double> new_dgammads_map) {
    dgammads_map = new_dgammads_map;
  }

  const std::unordered_map<TwoBodyKey, double> &
  get_two_body_flow_eq_map() const {
    return dgammads_map;
  }

  const std::unordered_map<OneBodyKey, int> &get_one_body_index() const {
    return f_index;
  }

  const void print_zero_body();
  const void print_one_body();
  const void print_two_body();

  const void print_one_body_generator();
  const void print_two_body_generator();

  const void print_one_body_flow_eq();
  const void print_two_body_flow_eq();
};
#endif // !HAMILTONIAN

#include "Hamiltonian.hpp"

#ifndef GENERATOR
#define GENERATOR

class Generator {
private:
  std::vector<OneBody> eta_1body;
  std::vector<TwoBody> eta_2body;

public:
  Generator();
  Generator(Hamiltonian H);
  const std::vector<OneBody> &get_one_body() const { return eta_1body; }
  const void set_one_body(std::vector<OneBody> new_eta_1body) {
    eta_1body = new_eta_1body;
  }
  const std::vector<TwoBody> &get_two_body() const { return eta_2body; }
  const void set_two_body(std::vector<TwoBody> new_eta_2body) {
    eta_2body = new_eta_2body;
  }
};

#endif // !GENERATOR

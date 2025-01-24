#include "Hamiltonian.hpp"

#ifndef GENERATOR
#define GENERATOR

class Generator {
private:
  std::vector<OneBody> eta_1body;
  std::vector<TwoBody> eta_2body;

public:
  Generator(Hamiltonian H);
  const std::vector<OneBody> &get_one_body() const { return eta_1body; }
  const std::vector<TwoBody> &get_two_body() const { return eta_2body; }
};

#endif // !GENERATOR

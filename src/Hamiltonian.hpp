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
  std::vector<TwoBody> Gamma;

public:
  Hamiltonian(std::string filepath);

  const std::vector<TwoBody> &get_two_body() const {
    return Gamma;
  }
  const std::vector<OneBody> &get_one_body() const {
    return f;
  }
  const double &get_zero_body() const {
    return E;
  }
  const int &get_num_states() const {
    return num_states;
  }

  void f0();
  void Gamma0();
};
#endif // !HAMILTONIAN

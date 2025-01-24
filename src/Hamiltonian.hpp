#include <string>
#include <vector>

#include "OneBody.hpp"
#include "TwoBody.hpp"



class Hamiltonian {
private:
  int num_states = 139;

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

  void f0();
  void Gamma0();
};

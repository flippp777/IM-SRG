#include "Hamiltonian.hpp"

class Generator {
private:
  std::vector<OneBody> eta_1body;
  std::vector<TwoBody> eta_2body;

  public:
  Generator(Hamiltonian H);
};

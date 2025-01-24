#include "Generator.hpp"

Generator::Generator(Hamiltonian H) {
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();

  
}

#include "Generator.hpp"
#include "Hamiltonian.hpp"

#include <vector>

#ifndef RUNGEKUTTA
#define RUNGEKUTTA

class RungeKutta {
private:
  std::vector<double> flow_equations;

public:
  RungeKutta(Hamiltonian H, Generator eta);
  
  void evolve(double s0, double s_max, double ds);

  const std::vector<double> &get_flow_equations() const {
    return flow_equations;
  }
};

#endif // !RUNGEKUTTA

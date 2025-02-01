#include "Integrate.hpp"

double euler_step(double x, const double dxds, const double ds) {
  return ds * dxds;
}

double rk4_step(const double x0, const double dxds, const double ds);


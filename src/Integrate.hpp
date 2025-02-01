#ifndef INTEGRATE
#define INTEGRATE

double euler_step(const double x0, const double dxds, const double ds);

double rk4_step(double x0, double dxds, double ds);

#endif // !INTEGRATE

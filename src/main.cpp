#include "Generator.hpp"
#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "RungeKutta.hpp"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <ostream>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  auto start = std::chrono::high_resolution_clock::now();

  Hamiltonian H = Hamiltonian(filepath);
  Generator eta = Generator(H);

  auto stop_constructors = std::chrono::high_resolution_clock::now();

  auto duration_constructors = std::chrono::duration_cast<std::chrono::seconds>(
      stop_constructors - start);
  std::cout << "Constructor H and eta time: " << duration_constructors.count()
            << " s" << std::endl;
  std::cout << std::endl;

  std::vector<OneBody> f = H.get_one_body();
  std::cout << "f:" << std::endl;
  for (const OneBody &ff : f) {
    std::cout << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
  }
  std::cout << std::endl;

  std::vector<OneBody> eta_1 = eta.get_one_body();
  std::cout << "eta_1body:" << std::endl;
  for (const OneBody &eta : eta_1) {
    std::cout << eta.s << "\t" << eta.ss << "\t" << eta.me << std::endl;
  }
  std::cout << std::endl;
  RungeKutta rk = RungeKutta(H, eta);

  std::cout << "E0 = " << H.get_zero_body() << std::endl;
  std::cout << std::endl;

  std::vector<double> flow_eq = rk.get_flow_equations();
  std::cout << "flow equations:" << std::endl;
  for (const double &eq : flow_eq) {
    std::cout << eq << std::endl;
  }
  std::cout << std::endl;

  return 0;
}

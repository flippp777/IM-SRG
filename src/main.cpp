#include "Hamiltonian.hpp"
#include "TwoBody.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>

int main() {
  std::string filepath = "matrix_elements/"
                         "hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_"
                         "hwHO020.me2m";

  Hamiltonian H = Hamiltonian(filepath);
  H.set_num_states(140);
  H.set_index_maps();

  // H.runge_kutta_gsl(0, 5, 1e-5);

  std::vector<TwoBody> g = H.get_two_body();
  const std::unordered_map<OneBodyKey, int> g_mat_index = H.get_one_body_index();

  

  // std::ofstream g_file ("Gamma_s0.0.txt");
  // for (const auto& gg : g)
  // {
  //   if (g_mat_index.find({gg.s1, gg.s2}) != g_mat_index.end() && g_mat_index.find({gg.ss1, gg.ss2}) != g_mat_index.end() )
  //   g_file << g_mat_index.at({gg.s1, gg.s2}) << "\t" << g_mat_index.at({gg.ss1, gg.ss2}) << "\t" << gg.me << "\n";
  // }
  


// std::vector<double> y = H.full_vector();
// H.set_elements(y.data());

// std::vector<OneBody> f = H.get_one_body();
// std::ofstream f_file("f_after_full_s0.0.txt");
// for (const auto& eta : f) {
//   f_file << eta.s << "\t" << eta.ss << "\t" << eta.me << "\n";
// }
// for (const auto &eta : f) {
//   std::cout << eta.s << "\t" << eta.ss << "\t" << eta.me << std::endl;
// }

// std::vector<TwoBody> gamma = H.get_two_body();
// for (const auto &g : gamma) {
//   std::cout << g.s1 << " " << g.s2 << " "
//             << g.ss1 << " " << g.ss2 << " " << g.me << std::endl;
// }

// H.set_flow_eq();
// H.get_full_eqs();

// H.set_generators();
// std::vector<TwoBody> eta2 = H.get_two_body_generator();
// std::ofstream eta_file ("etaTwoBody_s0.0.txt");
// for (const auto &eta : eta2) {
//   eta_file << eta.s1 << "\t" << eta.s2 << "\t" << eta.ss1 << "\t" <<
//   eta.ss2 << "\t" << eta.me << "\n";
// }

// H.set_flow_eq();
// std::unordered_map<OneBodyKey, double> dfds =
// H.get_one_body_flow_eq_map(); std::ofstream dfds_file
// ("dfds_s0.0.txt"); for (const auto &[key, me] : dfds) {
//   dfds_file << key.s << "\t" << key.ss << "\t" << me << "\n";
// }

// H.euler(0., 3, 0.001);

return 0;
}

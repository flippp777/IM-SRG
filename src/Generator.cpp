#include "Generator.hpp"
#include <iostream>

Generator::Generator(Hamiltonian H) {
  std::vector<OneBody> f = H.get_one_body();
  std::vector<TwoBody> Gamma = H.get_two_body();
  int num_states = H.get_num_states();

#pragma omp parallel for collapse(2)
  for (int i = 0; i < num_states; ++i) {
    for (int j = 0; j < num_states; ++j) {
      double f_ph = 0;
      for (int k = 0; k < f.size(); ++k) {
        if (f[k].s == i && f[k].s > 3 && f[k].ss == j && f[k].ss < 4) {
          f_ph += f[k].me;
        }
      }
      double delta = 0;
      for (int k = 0; k < f.size(); ++k) {
        if (f[k].s == f[k].ss && f[k].s > 3) {
          delta = f[k].me;
        } else if (f[k].s == f[k].ss && f[k].s < 3) {
          delta += -f[k].me;
        }
      }
      if (f_ph != 0 && delta != 0) {
        eta_1body.push_back(OneBody(i, j, f_ph / delta));
      } else if (delta == 0){
        std::cout << "delta = 0" << std::endl;
      }
    }
  }
}

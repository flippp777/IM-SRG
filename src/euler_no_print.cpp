#include "Hamiltonian.hpp"
#include <fstream>
#include <iostream>
#include <sstream>


void
Hamiltonian::euler_no_print (const double s0, const double s, const double ds)
{

  std::ostringstream E_filename;
  E_filename << "im_srg_euler_E_s0" << s0 << "_s" << s << "_ds" << ds
             << ".txt";
  std::ofstream E_file (E_filename.str ());

  int n = (int)((s - s0) / ds);

  std::unordered_map<OneBodyKey, double> f_map = vec_to_map (f);
  std::unordered_map<TwoBodyKey, double> gamma_map = vec_to_map (gamma);

  std::ostringstream name1;
  name1 << "im_srg_euler_f_s" << s0 << ".txt";
  std::ofstream f_file1 (name1.str ());
  for (const auto &ff : f)
    {
      f_file1 << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
    }

  set_generators ();
  set_flow_eq_print ();

  for (int i = 0; i < n; ++i)
    {
      set_generators ();
      set_flow_eq_print ();

      E_file << i*ds << "\t" << E << "\t" << dEds << std::endl;

      E += euler_step (dEds, ds);

      for (const auto &[key, me] : dfds_map)
        {
          if (f_map.count (key) != 0)
            {
              f_map[key] += euler_step (me, ds);
            }
          else
            {
              f_map.insert ({ key, euler_step (me, ds) });
            }
        }

      std::cout << "\t" << "Euler step for G";

      for (const auto &[key, me] : dgammads_map)
        {
          if (gamma_map.count (key) != 0)
            {
              gamma_map[key] += euler_step (me, ds);
            }
          else
            {
              gamma_map.insert ({ key, euler_step (me, ds) });
            }
        }

      f = map_to_vec (f_map);
      gamma = map_to_vec (gamma_map);

      if (i == (int)(n / 2))
        {
          std::ostringstream name2;
          name2 << "im_srg_euler_f_s" << i * ds << ".txt";
          std::ofstream f_file (name2.str ());
          for (const auto &ff : f)
            {
              f_file << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
            }
        }
    }

  std::ostringstream name3;
  name3 << "im_srg_euler_f_s" << s << ".txt";
  std::ofstream ff_file (name3.str ());
  for (const auto &ff : f)
    {
      ff_file << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
    }
}

#include "Hamiltonian.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
#include <ratio>
#include <sstream>
#include <string>

double
euler_step (const double dxds, const double ds)
{
  return ds * dxds;
}

void
Hamiltonian::euler (const double s0, const double s, const double ds)
{
  auto start = std::chrono::high_resolution_clock::now ();

  std::ostringstream E_filename;
  E_filename << "im_srg_euler_E_s0" << s0 << "_s" << s << "_ds" << ds
             << ".txt";
  std::ofstream E_file (E_filename.str ());

  int n = (int)((s - s0) / ds);
  std::cout << "##########   " << "s0 = " << s0 << "  <->  " << "s = " << s
            << "  <->  " << "ds = " << ds << "  <->  " << n << " steps"
            << "   ##########" << std::endl;
  std::cout << std::endl;

  std::cout << "Creating f and gamma maps" << std::flush;
  std::cout << "." << std::flush;

  std::unordered_map<OneBodyKey, double> f_map = vec_to_map (f);

  std::cout << "." << std::flush;

  std::unordered_map<TwoBodyKey, double> gamma_map = vec_to_map (gamma);

  std::cout << "." << std::flush;
  auto stop_map = std::chrono::high_resolution_clock::now ();
  std::chrono::duration<double, std::milli> ms_double = stop_map - start;
  std::cout << " done (" << ms_double.count () << "ms)" << std::endl;

  std::ostringstream name1;
  name1 << "im_srg_euler_f_s" << s0 << ".txt";
  std::ofstream f_file1 (name1.str ());
  for (const auto &ff : f)
    {
      f_file1 << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
    }

  set_generators ();
  set_flow_eq_print ();

      std::cout << std::endl;
  std::cout << "E = " << E << ",\t" << "dEds = " << dEds << ";\t"
            << "s = " << s0 << std::endl;
  std::cout << std::endl;
      std::cout << "f size: " << f.size () << std::endl;
      std::cout << "dfds size: " << dfds_map.size () << std::endl;
      std::cout << std::endl;
      std::cout << "G size: " << gamma.size () << std::endl;
      std::cout << "dGds size: " << dgammads_map.size () << std::endl;
      std::cout << std::endl;

  for (int i = 0; i <= n; ++i)
    {
      auto start_step = std::chrono::high_resolution_clock::now ();

      std::cout << "----------------   " << "Step " << i
                << "   ---------------" << std::endl;

      set_generators ();
      set_flow_eq_print ();

      E_file << s0 << "\t" << E << "\t" << dEds << std::endl;

      std::cout << "Calculating Euler steps" << std::endl;

      std::cout << "\t" << "Euler step for E";
      auto start_eE = std::chrono::high_resolution_clock::now ();

      E += euler_step (dEds, ds);

      auto stop_eE = std::chrono::high_resolution_clock::now ();
      std::chrono::duration<double, std::milli> ms_double = stop_eE - start_eE;
      std::cout << "..." << "\t" << "\t" << " done (" << ms_double.count ()
                << "ms)" << std::endl;

      std::cout << "\t" << "Euler step for f";

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

      auto stop_ef = std::chrono::high_resolution_clock::now ();
      ms_double = stop_ef - start_eE;
      std::cout << "..." << "\t"
                << "\t"
                   " done ("
                << ms_double.count () << "ms)" << std::endl;

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

      auto stop_eG = std::chrono::high_resolution_clock::now ();
      ms_double = stop_eG - stop_ef;
      std::cout << "..." << "\t" << "\t" << " done (" << ms_double.count ()
                << "ms)" << std::endl;

      std::cout << "\t" << "f in vector";

      f = map_to_vec (f_map);

      auto stop_f = std::chrono::high_resolution_clock::now ();
      ms_double = stop_f - stop_eG;
      std::cout << "..." << "\t" << "\t" << "\t" << " done ("
                << ms_double.count () << "ms)" << std::endl;

      std::cout << "\t" << "G in vector";

      gamma = map_to_vec (gamma_map);

      auto stop_G = std::chrono::high_resolution_clock::now ();
      ms_double = stop_G - stop_f;
      std::cout << "..." << "\t" << "\t" << "\t" << " done ("
                << ms_double.count () << "ms)" << std::endl;

      ms_double = stop_G - start_eE;
      std::cout << "\t" << "|-> done(" << ms_double.count () << "ms)"
                << std::endl;

      std::cout << std::endl;
      std::cout << "E = " << E << ",\t" << "dEds = " << dEds << ";\t"
                << "s = " << (i + 1) * ds << std::endl;
      std::cout << std::endl;
      std::cout << "f size: " << f.size () << std::endl;
      std::cout << "dfds size: " << dfds_map.size () << std::endl;
      std::cout << std::endl;
      std::cout << "G size: " << gamma.size () << std::endl;
      std::cout << "dGds size: " << dgammads_map.size () << std::endl;
      std::cout << std::endl;

      ms_double = stop_G - start_step;
      std::cout << "Step time: " << ms_double.count () * 1 / 1000 << "s"
                << std::endl;
      std::cout << std::endl;

      if (i == (int)n / 2)
        {
          std::ostringstream name2;
          name2 << "im_srg_euler_f_s" << i * ds << ".txt";
          std::ofstream f_file (name1.str ());
          for (const auto &ff : f)
            {
              f_file << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
            }
        }
    }
  std::cout << std::endl;
  auto stop = std::chrono::high_resolution_clock::now ();
  ms_double = stop - start;
  std::cout << "Flow time: " << ms_double.count () * 1 / 1000 << "s"
            << std::endl;

  std::ostringstream name3;
  name3 << "im_srg_euler_f_s" << s << ".txt";
  std::ostringstream ff_file (name3.str ());
  for (const auto &ff : f)
    {
      ff_file << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
    }
}

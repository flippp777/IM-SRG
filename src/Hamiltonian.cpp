#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <ratio>
#include <sstream>
#include <tuple>
#include <vector>

Hamiltonian::Hamiltonian () : E (), f (), gamma () {}

Hamiltonian::Hamiltonian (const double E, const std::vector<OneBody> f,
                          const std::vector<TwoBody> gamma)
    : E (E), f (f), gamma (gamma)
{
}

Hamiltonian::Hamiltonian (std::string filename)
{
  std::cout << "Constructing Hamiltonian";
  auto start = std::chrono::high_resolution_clock::now ();

  // Open file
  std::ifstream file (filename);
  if (!file)
    {
      std::cerr << "error: could not open file" << std::endl;
      return;
    }
  std::cout << "." << std::flush;

  // save two-body matrix elemetns in gamma
  gamma.reserve (7300920);
  std::string line;
  int i = 0;
  while (std::getline (file, line))
    {
      std::istringstream iss (line);
      if (i > 144)
        {
          TwoBody gam;
          iss >> gam.s1 >> gam.s2 >> gam.ss1 >> gam.ss2 >> gam.me;
          gamma.push_back (gam);
        }
      i++;
    }
  file.close ();

  std::cout << "." << std::flush;

  // full gamma
  const int gamma_init_size = gamma.size ();
  std::vector<TwoBody> gamma_full;
  gamma_full.reserve (7500000);
  for (const auto &g : gamma)
    {
      if (g.s1 == g.ss1 && g.s2 == g.ss2)
        {
          // <ab|H|ab>
          gamma_full.push_back (TwoBody (g.s1, g.s2, g.ss1, g.ss2, g.me));
          // <ba|H|ab>
          gamma_full.push_back (TwoBody (g.s2, g.s1, g.ss1, g.ss2, -g.me));
          // <ba|H|ba>
          gamma_full.push_back (TwoBody (g.s2, g.s1, g.ss2, g.ss1, g.me));
          // <ab|H|ba>
          gamma_full.push_back (TwoBody (g.s1, g.s2, g.ss2, g.ss1, -g.me));
        }
      else
        {
          // <ab|H|cd>
          gamma_full.push_back (TwoBody (g.s1, g.s2, g.ss1, g.ss2, g.me));
          // <ba|H|cd>
          gamma_full.push_back (TwoBody (g.s2, g.s1, g.ss1, g.ss2, -g.me));
          // <ba|H|dc>
          gamma_full.push_back (TwoBody (g.s2, g.s1, g.ss2, g.ss1, g.me));
          // <ab|H|dc>
          gamma_full.push_back (TwoBody (g.s1, g.s2, g.ss2, g.ss1, -g.me));
          // <cd|H|ab>
          gamma_full.push_back (TwoBody (g.ss1, g.ss2, g.s1, g.s2, g.me));
          // <dc|H|ab>
          gamma_full.push_back (TwoBody (g.ss2, g.ss1, g.s1, g.s2, -g.me));
          // <dc|H|ba>
          gamma_full.push_back (TwoBody (g.ss2, g.ss1, g.s2, g.s1, g.me));
          // <cd|H|ba>
          gamma_full.push_back (TwoBody (g.ss1, g.ss2, g.s2, g.s1, -g.me));
        }
    }
  gamma = gamma_full;

  std::cout << "." << std::flush;

  std::vector<OneBody> f_new;
  const int gamma_size = gamma.size ();

#pragma omp parallel
  {
#pragma omp for reduction(+ : E)
    for (int i = 0; i < gamma_size; ++i)
      {
        if (gamma[i].s1 < 4 && gamma[i].s2 < 4 && gamma[i].s1 == gamma[i].ss1
            && gamma[i].s2 == gamma[i].ss2)
          {
            E += 0.5 * gamma[i].me;
          }
      }

      // calculate and set one-body
#pragma omp for
    for (int i = 0; i < gamma_size; ++i)
      {
        if (gamma[i].s2 < 4 && gamma[i].s2 == gamma[i].ss2)
          {
#pragma omp critical
            f = f + OneBody (gamma[i].s1, gamma[i].ss1, gamma[i].me);
          }
      }
      // remove values close to 0
#pragma omp for
    for (int i = 0; i < f.size (); ++i)
      {
        if (std::abs (f[i].me) > 1e-9)
          {
#pragma omp critical
            f_new.push_back (f[i]);
          }
      }
  }
  f = f_new;
  auto stop = std::chrono::high_resolution_clock::now ();
  std::chrono::duration<double, std::milli> ms_double = stop - start;
  std::cout << "\t" << " done (" << ms_double.count () << "ms)" << std::endl;
  ;
}

const void
Hamiltonian::set_full_vector ()
{
  std::cout << "Creating the full vector" << std::flush;

  std::unordered_map<OneBodyKey, double> f_map;
  for (const auto &el : f)
    {
      f_map[{ el.s, el.ss }] = el.me;
    }
  std::unordered_map<TwoBodyKey, double> gamma_map;
  for (const auto &el : gamma)
    {
      gamma_map[{ el.s1, el.s2, el.ss1, el.ss2 }] = el.me;
    }

  std::vector<double> f_full (num_states * num_states, 0.0);
  std::vector<double> gamma_full (num_states * num_states * num_states * num_states, 0.0);

#pragma omp parallel
  {
#pragma omp for collapse(2) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        f_full[i * num_states + j] = f_map.count ({ i, j }) ? f_map[{ i, j }] : 0;

#pragma omp for collapse(4) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        for (int k = 0; k < num_states; ++k)
          for (int l = 0; l < num_states; ++l)
            gamma_full[i * num_states * num_states * num_states + j * num_states * num_states + k * num_states + l]
                = gamma_map.count ({ i, j, k, l }) ? gamma_map[{ i, j, k, l }]
                                                   : 0;

#pragma omp barrier
#pragma omp single
    {
      y.resize (1 + num_states * num_states + num_states * num_states * num_states * num_states);
      y[0] = E;
    }

#pragma omp for collapse(2) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        y[1 + i * num_states + j] = f_full[i * num_states + j];

#pragma omp for collapse(4) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        for (int k = 0; k < num_states; ++k)
          for (int l = 0; l < num_states; ++l)
            y[1 + num_states * num_states
              + (i * num_states * num_states * num_states + j * num_states * num_states + k * num_states + l)]
                = gamma_full[i * num_states * num_states * num_states + j * num_states * num_states + k * num_states
                             + l];
  }
  printf ("... done\n");
}

const void
Hamiltonian::set_elements (const double y[])
{
  f.clear ();
  f.shrink_to_fit ();
  gamma.clear ();

  f.reserve (num_states * num_states);
  gamma.reserve (20000000);

  E = y[0];

#pragma omp parallel
  {
#pragma omp for collapse(2) nowait
    for (int i = 0; i < num_states; ++i)
      {
        for (int j = 0; j < num_states; ++j)
          {
            if (y[1 + i * num_states + j] != 0)
              {
#pragma omp critical
                f.push_back (OneBody (i, j, y[1 + i * num_states + j]));
              }
          }
      }
#pragma omp for collapse(4) nowait
    for (int i = 0; i < num_states; ++i)
      {
        for (int j = 0; j < num_states; ++j)
          {
            for (int k = 0; k < num_states; ++k)
              {
                for (int l = 0; l < num_states; ++l)
                  {
                    if (y[1 + num_states * num_states
                          + (i * num_states * num_states * num_states + j * num_states * num_states + k * num_states
                             + l)]
                        != 0)
                      {
#pragma omp critical
                        gamma.push_back (
                            TwoBody (i, j, k, l,
                                     y[1 + num_states * num_states
                                       + (i * num_states * num_states * num_states + j * num_states * num_states
                                          + k * num_states + l)]));
                      }
                  }
              }
          }
      }
  }
}

std::vector<double>
Hamiltonian::get_full_eqs ()
{
  std::vector<double> dyds;
  std::vector<double> dfds_full;
  std::vector<double> dgammads_full;

  dfds_full.resize (num_states * num_states, 0.0);
  dgammads_full.resize (num_states * num_states * num_states * num_states, 0.0);

#pragma omp parallel
  {
#pragma omp for collapse(2) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        dfds_full[i * num_states + j]
            = dfds_map.count ({ i, j }) ? dfds_map[{ i, j }] : 0;

#pragma omp for collapse(4) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        for (int k = 0; k < num_states; ++k)
          for (int l = 0; l < num_states; ++l)
            dgammads_full[i * num_states * num_states * num_states + j * num_states * num_states + k * num_states + l]
                = dgammads_map.count ({ i, j, k, l })
                      ? dgammads_map[{ i, j, k, l }]
                      : 0;

#pragma omp singe
    {
      dyds.resize (1 + num_states * num_states + num_states * num_states * num_states * num_states);
      dyds[0] = E;
    }
#pragma omp barrier

#pragma omp for collapse(2) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        dyds[1 + i * num_states + j] = dfds_full[i * num_states + j];

#pragma omp for collapse(4) nowait
    for (int i = 0; i < num_states; ++i)
      for (int j = 0; j < num_states; ++j)
        for (int k = 0; k < num_states; ++k)
          for (int l = 0; l < num_states; ++l)
            dyds[1 + num_states * num_states
                 + (i * num_states * num_states * num_states + j * num_states * num_states + k * num_states + l)]
                = dgammads_full[i * num_states * num_states * num_states + j * num_states * num_states + k * num_states
                                + l];
  }
  return dyds;
}

const void
Hamiltonian::print_zero_body ()
{
  std::cout << "E: " << E << std::endl;
}

const void
Hamiltonian::print_one_body ()
{
  std::cout << "f: " << std::endl;
  for (const OneBody &ff : f)
    {
      std::cout << ff.s << "\t" << ff.ss << "\t" << ff.me << std::endl;
    }
  std::cout << "f size: " << f.size () << std::endl;
  std::cout << std::endl;
}

const void
Hamiltonian::print_two_body ()
{
  std::cout << "gamma:" << std::endl;
  for (const TwoBody &gg : gamma)
    {
      std::cout << gg.s1 << "\t" << gg.s2 << "\t" << gg.ss1 << "\t" << gg.ss2
                << "\t" << gg.me << std::endl;
    }
  std::cout << "gamma size: " << gamma.size () << std::endl;
  std::cout << std::endl;
}

const void
Hamiltonian::print_one_body_generator ()
{
  std::cout << "f: " << std::endl;
  for (const OneBody &eta1 : eta_1body)
    {
      std::cout << eta1.s << "\t" << eta1.ss << "\t" << eta1.me << std::endl;
    }
  std::cout << std::endl;
}

const void
Hamiltonian::print_two_body_generator ()
{
  std::cout << "Gamma:" << std::endl;
  for (const TwoBody &eta2 : eta_2body)
    {
      std::cout << eta2.s1 << "\t" << eta2.s2 << "\t" << eta2.ss1 << "\t"
                << eta2.ss2 << "\t" << eta2.me << std::endl;
    }
  std::cout << std::endl;
}

const void
Hamiltonian::print_one_body_flow_eq ()
{
  std::cout << "f: " << std::endl;
  for (const OneBody &df : dfds_vec)
    {
      std::cout << df.s << "\t" << df.ss << "\t" << df.me << std::endl;
    }
  std::cout << std::endl;
}

const void
Hamiltonian::print_two_body_flow_eq ()
{
  std::cout << "Gamma:" << std::endl;
  for (const TwoBody &dg : dgammads_vec)
    {
      std::cout << dg.s1 << "\t" << dg.s2 << "\t" << dg.ss1 << "\t" << dg.ss2
                << "\t" << dg.me << std::endl;
    }
  std::cout << std::endl;
}

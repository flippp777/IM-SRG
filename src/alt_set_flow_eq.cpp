#include "Hamiltonian.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include "TwoBodyQA.hpp"
#include "reduce_lists.hpp"
#include <chrono>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <vector>

int
Hamiltonian::flow_eq_gsl ()
{
  dEds = 0;
  dfds_vec.clear ();
  dfds_map.clear ();
  dgammads_vec.clear ();
  dgammads_map.clear ();

  std::cout << "Constructing flow equations" << std::endl;
  auto start = std::chrono::high_resolution_clock::now ();

  std::cout << "\t" << "Creating maps" << std::flush;

  TwoBodyQA gamma_qa = TwoBodyQA (gamma, 1);
  std::cout << "." << std::flush;
  TwoBodyQA eta2_qa = TwoBodyQA (eta_2body, 1);

  std::cout << "." << std::flush;

  const int f_size = f.size ();
  const int eta1_size = eta_1body.size ();
  const int gamma_size = gamma.size ();
  const int eta2_size = eta_2body.size ();

  std::cout << "." << std::flush;

  int max_threads = omp_get_max_threads ();
  std::vector<std::vector<OneBody> > thread_dfds (max_threads);
  std::vector<std::vector<TwoBody> > thread_dgammads (max_threads);

  auto stop_maps = std::chrono::high_resolution_clock::now ();
  std::chrono::duration<double, std::milli> ms_double = stop_maps - start;
  std::cout << "\t\t" << " done (" << ms_double.count () << "ms)" << std::endl;

  std::cout << "\t" << "Constructing new flow eqs" << std::flush;

#pragma omp parallel
  {
//////////////////////////////// dEds
#pragma omp for reduction(+ : dEds)
    for (int i = 0; i < f_size; ++i)
      {
        for (int j = 0; j < eta1_size; ++j)
          {
            if (f[i].s == eta_1body[j].ss && f[i].ss == eta_1body[j].s
                && f[i].s > 3 && f[i].ss < 4)
              {
                dEds += f[i].me * eta_1body[j].me;
              }
            else if (f[i].s == eta_1body[j].ss && f[i].ss == eta_1body[j].s
                     && f[i].s < 4 && f[i].s > 3)
              {
                dEds += -f[i].me * eta_1body[j].me;
              }
          }
      }

#pragma omp for reduction(+ : dEds) nowait
    for (int i = 0; i < eta2_size; ++i)
      {
        if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4 && eta_2body[i].ss1 > 3
            && eta_2body[i].ss2 > 3)
          {
            dEds += 0.5 * eta_2body[i].me
                    * gamma_qa
                          .map[{ eta_2body[i].ss1, eta_2body[i].ss2,
                                 eta_2body[i].s1, eta_2body[i].s2 }]
                          ->me;
          }
      }

#pragma omp single
    std::cout << "." << std::flush;
    ///////////////////////////////////////////////// dfds
    // first sum
    int tid = omp_get_thread_num ();
    std::vector<OneBody> &local_dfds = thread_dfds[tid];
    local_dfds.reserve (f_size);

#pragma omp for collapse(2)
    for (int i = 0; i < f_size; ++i)
      {
        for (int j = 0; j < eta1_size; ++j)
          {
            // first sum
            if (eta_1body[j].ss == f[i].s)
              {
                local_dfds.push_back (OneBody (eta_1body[j].s, f[i].ss,
                                               eta_1body[j].me * f[i].me));
                local_dfds.push_back (OneBody (f[i].ss, eta_1body[j].s,
                                               eta_1body[j].me * f[i].me));
              }
          }
      }

// second sum
#pragma omp for
    for (int i = 0; i < eta1_size; ++i)
      {
        if (eta_1body[i].s < 4 && eta_1body[i].ss > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_1body[i].ss, eta_1body[i].s))
              {
                local_dfds.push_back (
                    OneBody (g->s2, g->ss2, eta_1body[i].me * g->me));
              }
          }
        else if (eta_1body[i].s > 3 && eta_1body[i].ss < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_1body[i].ss, eta_1body[i].s))
              {
                local_dfds.push_back (
                    OneBody (g->s2, g->ss2, -eta_1body[i].me * g->me));
              }
          }
      }

#pragma omp for
    for (int i = 0; i < f_size; ++i)
      {
        if (f[i].s < 4 && f[i].ss > 3)
          {
            for (const auto &e : eta2_qa.check_s1ss1 (f[i].ss, f[i].s))
              {
                local_dfds.push_back (
                    OneBody (e->s2, e->ss2, -f[i].me * e->me));
              }
          }
        else if (f[i].s > 3 && f[i].ss < 4)
          {
            for (const auto &e : gamma_qa.check_s1ss1 (f[i].ss, f[i].s))
              {
                local_dfds.push_back (
                    OneBody (e->s2, e->ss2, f[i].me * e->me));
              }
          }
      }

#pragma omp single
    std::cout << "." << std::flush;

// third sum
#pragma omp for
    for (int i = 0; i < eta2_size; ++i)
      {
        if (eta_2body[i].ss1 < 4 && eta_2body[i].ss2 < 4
            && eta_2body[i].s1 > 3)
          {
            for (const auto &g : gamma_qa.check_s1s2ss1 (
                     eta_2body[i].ss1, eta_2body[i].ss2, eta_2body[i].s1))
              {
                local_dfds.push_back (OneBody (eta_2body[i].s2, g->ss2,
                                               0.5 * eta_2body[i].me * g->me));
                local_dfds.push_back (OneBody (g->ss2, eta_2body[i].s2,
                                               0.5 * eta_2body[i].me * g->me));
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3
                 && eta_2body[i].s1 < 4)
          {
            for (const auto &g : gamma_qa.check_s1s2ss1 (
                     eta_2body[i].ss1, eta_2body[i].ss2, eta_2body[i].s1))
              {
                local_dfds.push_back (OneBody (
                    eta_2body[i].s2, g->ss2, -0.5 * eta_2body[i].me * g->me));
                local_dfds.push_back (OneBody (
                    g->ss2, eta_2body[i].s2, -0.5 * eta_2body[i].me * g->me));
              }
          }
      }

    ///////////////////////////// dGammads
    std::vector<TwoBody> &local_dgammads = thread_dgammads[tid];
    local_dgammads.reserve (gamma_size);

// first sum
// gamma terms
#pragma omp for
    for (int i = 0; i < eta1_size; ++i)
      {
        for (const auto &g : gamma_qa.check_s1 (eta_1body[i].ss))
          {
            local_dgammads.push_back (TwoBody (eta_1body[i].s, g->s2, g->ss1,
                                               g->ss2,
                                               eta_1body[i].me * g->me));
            local_dgammads.push_back (TwoBody (g->s2, eta_1body[i].s, g->ss1,
                                               g->ss2,
                                               -eta_1body[i].me * g->me));
          }
        for (const auto &g : gamma_qa.check_ss1 (eta_1body[i].s))
          {
            local_dgammads.push_back (TwoBody (g->s1, g->s2, eta_1body[i].ss,
                                               g->ss2,
                                               -eta_1body[i].me * g->me));
            local_dgammads.push_back (TwoBody (g->s1, g->s2, g->ss2,
                                               eta_1body[i].ss,
                                               eta_1body[i].me * g->me));
          }
      }

#pragma omp single
    std::cout << "." << std::flush;

// eta2 terms
#pragma omp for
    for (int i = 0; i < f_size; ++i)
      {
        for (const auto &e : eta2_qa.check_s1 (f[i].ss))
          {
            local_dgammads.push_back (
                TwoBody (f[i].s, e->s2, e->ss1, e->ss2, -f[i].me * e->me));
            local_dgammads.push_back (
                TwoBody (e->s2, f[i].s, e->ss1, e->ss2, f[i].me * e->me));
          }
        for (const auto &e : eta2_qa.check_ss1 (f[i].s))
          {
            local_dgammads.push_back (
                TwoBody (e->s1, e->s2, f[i].ss, e->ss2, f[i].me * e->me));
            local_dgammads.push_back (
                TwoBody (e->s1, e->s2, e->ss2, f[i].ss, -f[i].me * e->me));
          }
      }

      // 2. and 3. sum

#pragma omp for
    for (int i = 0; i < eta2_size; ++i)
      {
        // 2. sum
        if (eta_2body[i].ss1 < 4 && eta_2body[i].ss2 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1s2 (eta_2body[i].ss1, eta_2body[i].ss2))
              {
                local_dgammads.push_back (
                    TwoBody (eta_2body[i].s1, eta_2body[i].s2, g->ss1, g->ss2,
                             -0.5 * eta_2body[i].me * g->me));
              }
          }
        else if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_ss1ss2 (eta_2body[i].s1, eta_2body[i].s2))
              {
                local_dgammads.push_back (
                    TwoBody (g->s1, g->s2, eta_2body[i].ss1, eta_2body[i].ss2,
                             0.5 * eta_2body[i].me * g->me));
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1s2 (eta_2body[i].ss1, eta_2body[i].ss2))
              {
                local_dgammads.push_back (
                    TwoBody (eta_2body[i].s1, eta_2body[i].s2, g->ss1, g->ss2,
                             0.5 * eta_2body[i].me * g->me));
              }
          }
        else if (eta_2body[i].s1 > 3 && eta_2body[i].s2 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_ss1ss2 (eta_2body[i].s1, eta_2body[i].s2))
              {
                local_dgammads.push_back (
                    TwoBody (g->s1, g->s2, eta_2body[i].ss1, eta_2body[i].ss2,
                             -0.5 * eta_2body[i].me * g->me));
              }

            // 3. sum
          }
        else if (eta_2body[i].ss1 < 4 && eta_2body[i].s1 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_2body[i].ss1, eta_2body[i].s1))
              {
                // 1
                local_dgammads.push_back (TwoBody (g->s2, eta_2body[i].s2,
                                                   g->ss2, eta_2body[i].ss2,
                                                   -eta_2body[i].me * g->me));
                // -P_12
                local_dgammads.push_back (TwoBody (eta_2body[i].s2, g->s2,
                                                   g->ss2, eta_2body[i].ss2,
                                                   eta_2body[i].me * g->me));
                // -P_34
                local_dgammads.push_back (TwoBody (g->s2, eta_2body[i].s2,
                                                   eta_2body[i].ss2, g->ss2,
                                                   eta_2body[i].me * g->me));
                // P_12 P_34
                local_dgammads.push_back (TwoBody (eta_2body[i].s2, g->s2,
                                                   eta_2body[i].ss2, g->ss2,
                                                   -eta_2body[i].me * g->me));
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].s1 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_2body[i].ss1, eta_2body[i].s1))
              {
                // 1
                local_dgammads.push_back (TwoBody (g->s2, eta_2body[i].s2,
                                                   g->ss2, eta_2body[i].ss2,
                                                   eta_2body[i].me * g->me));
                // -P_12
                local_dgammads.push_back (TwoBody (eta_2body[i].s2, g->s2,
                                                   g->ss2, eta_2body[i].ss2,
                                                   -eta_2body[i].me * g->me));
                // -P_34
                local_dgammads.push_back (TwoBody (g->s2, eta_2body[i].s2,
                                                   eta_2body[i].ss2, g->ss2,
                                                   -eta_2body[i].me * g->me));
                // P_12 P_34
                local_dgammads.push_back (TwoBody (eta_2body[i].s2, g->s2,
                                                   eta_2body[i].ss2, g->ss2,
                                                   eta_2body[i].me * g->me));
              }
          }
      }
  }

  auto stop_dgds = std::chrono::high_resolution_clock::now ();
  ms_double = stop_dgds - stop_maps;
  std::cout << "\t" << " done (" << ms_double.count () << "ms)" << std::endl;

  std::cout << "\t" << "Combining thread vectors" << std::flush;
  std::cout << "." << std::flush;

  dfds_vec.reserve (thread_dfds.size () * thread_dfds[0].size ());
  size_t dgammads_size = 0;
  for (const auto &local_vec : thread_dgammads)
    {
      dgammads_size += local_vec.size ();
    }

// #pragma omp parallel sections
//   {
// #pragma omp section
//     {
      for (const auto &local_map : thread_dfds)
        {
          dfds_vec.insert (dfds_vec.end (), local_map.begin (),
                           local_map.end ());
        }
      dfds_vec = sum_equal_elements (dfds_vec);
      dfds_map.reserve(dfds_vec.size());
    // }
    dgammads_vec.reserve (dgammads_size);
// #pragma omp section
//     {
      std::cout << "." << std::flush;
      for (const auto &local_map : thread_dgammads)
        {
          if (local_map.empty ())
            {
              std::cerr << "Empty local_map found!" << std::endl;
            }
          dgammads_vec.insert (dgammads_vec.end (), local_map.begin (),
                               local_map.end ());
        }
      dgammads_vec = sum_equal_elements (dgammads_vec);
    // }
  // }

  // dfds_vec = reduce_vec(thread_dfds);
  //     std::cout << "." << std::flush;
  // dgammads_vec = reduce_vec(thread_dgammads);

  std::cout << "." << std::flush;
  auto stop_comb = std::chrono::high_resolution_clock::now ();
  ms_double = stop_comb - stop_dgds;
  std::cout << "\t" << " done (" << ms_double.count () << "ms)" << std::endl;

  ms_double = stop_comb - start;
  std::cout << "\t" << "|-> done(" << ms_double.count () << "ms)" << std::endl;
  return GSL_SUCCESS;
}

// int
// Hamiltonian::flow_eq_gsl_wrapper (double s, const double y[], double dydt[],
//                                   void *params)
// {
//   // Cast the 'params' to a pointer to an instance of the class
//   Hamiltonian *h = static_cast<Hamiltonian *> (params);
//   return h->flow_eq_gsl (s, y, dydt, params);
// }

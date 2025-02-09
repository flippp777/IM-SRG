#include "Hamiltonian.hpp"
#include "TwoBodyQA.hpp"
#include "omp.h"
#include <iostream>
#include <omp.h>
#include <chrono>
#include <ostream>
#include <unordered_map>

void const
Hamiltonian::set_flow_eq ()
{
  dEds = 0;
  dfds_vec.clear ();
  dfds_map.clear ();
  dgammads_vec.clear ();
  dgammads_map.clear ();

  TwoBodyQA gamma_qa = TwoBodyQA (gamma, 1);
  TwoBodyQA eta2_qa = TwoBodyQA (eta_2body, 1);

  std::unordered_map<OneBodyKey, double> f_map;
  f_map.reserve (f.size ());
  for (const auto &ff : f)
    {
      f_map[OneBodyKey (ff.s, ff.ss)] = ff.me;
    }

  std::unordered_map<OneBodyKey, double> eta1_map;
  eta1_map.reserve (eta_1body.size ());
  for (const auto &e : eta_1body)
    {
      eta1_map[OneBodyKey (e.s, e.ss)] = e.me;
    }

  const int f_size = f.size ();
  const int eta1_size = eta_1body.size ();
  const int gamma_size = gamma.size ();
  const int eta2_size = eta_2body.size ();

  int max_threads = omp_get_max_threads ();
  std::vector<std::unordered_map<OneBodyKey, double> > thread_dfds (
      max_threads);
  std::vector<std::unordered_map<TwoBodyKey, double> > thread_dgammads (
      max_threads);
  for (int i = 0; i < max_threads; i++)
    {
      thread_dfds[i].reserve (5000);
      thread_dgammads[i].reserve (5000);
    }

#pragma omp parallel
  {
//////////////////////////////// dEds
#pragma omp for reduction(+ : dEds)
    for (int i = 0; i < f_size; ++i)
      {
        if (f[i].s > 3 && f[i].ss < 4
            && eta1_map.find ({ f[i].ss, f[i].s }) != eta1_map.end ())
          {
            dEds += f[i].me * eta1_map[OneBodyKey (f[i].ss, f[i].s)];
          }
        else if (f[i].s < 4 && f[i].s > 3
                 && eta1_map.find ({ f[i].ss, f[i].s }) != eta1_map.end ())
          {
            dEds += -f[i].me * eta1_map[OneBodyKey (f[i].ss, f[i].s)];
          }
      }
#pragma omp for reduction(+ : dEds) nowait
    for (int i = 0; i < eta2_size; ++i)
      {
        if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4 && eta_2body[i].ss1 > 3
            && eta_2body[i].ss2 > 3)
          {
            if (gamma_qa.map.find ({ eta_2body[i].ss1, eta_2body[i].ss2,
                                     eta_2body[i].s1, eta_2body[i].s2 })
                != gamma_qa.map.end ())
              {
                dEds += 0.5 * eta_2body[i].me
                        * gamma_qa
                              .map[{ eta_2body[i].ss1, eta_2body[i].ss2,
                                     eta_2body[i].s1, eta_2body[i].s2 }]
                              ->me;
              }
          }
      }

    // first sum
    int tid = omp_get_thread_num ();
    std::unordered_map<OneBodyKey, double> &local_dfds = thread_dfds[tid];

#pragma omp for collapse(2)
    for (int i = 0; i < f_size; ++i)
      {
        for (int j = 0; j < eta1_size; ++j)
          {
            // first sum
            if (eta_1body[j].ss == f[i].s)
              {
                local_dfds[{ eta_1body[j].s, f[i].ss }]
                    += eta_1body[j].me * f[i].me;
                local_dfds[{ f[i].ss, eta_1body[j].s }]
                    += eta_1body[j].me * f[i].me;
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
                local_dfds[{ g->s2, g->ss2 }] += eta_1body[i].me * g->me;
              }
          }
        else if (eta_1body[i].s > 3 && eta_1body[i].ss < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_1body[i].ss, eta_1body[i].s))
              {
                local_dfds[{ g->s2, g->ss2 }] += -eta_1body[i].me * g->me;
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
                local_dfds[{ e->s2, e->ss2 }] += -f[i].me * e->me;
              }
          }
        else if (f[i].s > 3 && f[i].ss < 4)
          {
            for (const auto &e : eta2_qa.check_s1ss1 (f[i].ss, f[i].s))
              {
                local_dfds[{ e->s2, e->ss2 }] += f[i].me * e->me;
              }
          }
      }

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
                local_dfds[{ eta_2body[i].s2, g->ss2 }]
                    += 0.5 * eta_2body[i].me * g->me;
                local_dfds[{ g->ss2, eta_2body[i].s2 }]
                    += 0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3
                 && eta_2body[i].s1 < 4)
          {
            for (const auto &g : gamma_qa.check_s1s2ss1 (
                     eta_2body[i].ss1, eta_2body[i].ss2, eta_2body[i].s1))
              {
                local_dfds[{ eta_2body[i].s2, g->ss2 }]
                    += -0.5 * eta_2body[i].me * g->me;
                local_dfds[{ g->ss2, eta_2body[i].s2 }]
                    += -0.5 * eta_2body[i].me * g->me;
              }
          }
      }
    ///////////////////////////// dGammads
    std::unordered_map<TwoBodyKey, double> &local_dgammads
        = thread_dgammads[tid];

// first sum
// gamma terms
#pragma omp for
    for (int i = 0; i < eta1_size; ++i)
      {
        for (const auto &g : gamma_qa.check_s1 (eta_1body[i].ss))
          {
            local_dgammads[{ eta_1body[i].s, g->s2, g->ss1, g->ss2 }]
                += eta_1body[i].me * g->me;
            local_dgammads[{ g->s2, eta_1body[i].s, g->ss1, g->ss2 }]
                += -eta_1body[i].me * g->me;
          }
        for (const auto &g : gamma_qa.check_ss1 (eta_1body[i].s))
          {
            local_dgammads[{ g->s1, g->s2, eta_1body[i].ss, g->ss2 }]
                += -eta_1body[i].me * g->me;
            local_dgammads[{ g->s1, g->s2, g->ss2, eta_1body[i].ss }]
                += eta_1body[i].me * g->me;
          }
      }
// eta2 terms
#pragma omp for
    for (int i = 0; i < f_size; ++i)
      {
        for (const auto &e : eta2_qa.check_s1 (f[i].ss))
          {
            local_dgammads[{ f[i].s, e->s2, e->ss1, e->ss2 }]
                += -f[i].me * e->me;
            local_dgammads[{ e->s2, f[i].s, e->ss1, e->ss2 }]
                += f[i].me * e->me;
          }
        for (const auto &e : eta2_qa.check_ss1 (f[i].s))
          {
            local_dgammads[{ e->s1, e->s2, f[i].ss, e->ss2 }]
                += f[i].me * e->me;
            local_dgammads[{ e->s1, e->s2, e->ss2, f[i].ss }]
                += -f[i].me * e->me;
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
                local_dgammads[{ eta_2body[i].s1, eta_2body[i].s2, g->ss1,
                                 g->ss2 }]
                    += -0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_ss1ss2 (eta_2body[i].s1, eta_2body[i].s2))
              {
                local_dgammads[{ g->s1, g->s2, eta_2body[i].ss1,
                                 eta_2body[i].ss2 }]
                    += 0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1s2 (eta_2body[i].ss1, eta_2body[i].ss2))
              {
                local_dgammads[{ eta_2body[i].s1, eta_2body[i].s2, g->ss1,
                                 g->ss2 }]
                    += 0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].s1 > 3 && eta_2body[i].s2 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_ss1ss2 (eta_2body[i].s1, eta_2body[i].s2))
              {
                local_dgammads[{ g->s1, g->s2, eta_2body[i].ss1,
                                 eta_2body[i].ss2 }]
                    += -0.5 * eta_2body[i].me * g->me;
              }
            // 3. sum
          }
        else if (eta_2body[i].ss1 < 4 && eta_2body[i].s1 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_2body[i].ss1, eta_2body[i].s1))
              {
                // 1
                local_dgammads[{ g->s2, eta_2body[i].s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += -eta_2body[i].me * g->me;
                // -P_12
                local_dgammads[{ eta_2body[i].s2, g->s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += eta_2body[i].me * g->me;
                // -P_34
                local_dgammads[{ g->s2, eta_2body[i].s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += eta_2body[i].me * g->me;
                // P_12 P_34
                local_dgammads[{ eta_2body[i].s2, g->s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += -eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].s1 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_2body[i].ss1, eta_2body[i].s1))
              {
                // 1
                local_dgammads[{ g->s2, eta_2body[i].s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += eta_2body[i].me * g->me;
                // -P_12
                local_dgammads[{ eta_2body[i].s2, g->s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += -eta_2body[i].me * g->me;
                // -P_34
                local_dgammads[{ g->s2, eta_2body[i].s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += -eta_2body[i].me * g->me;
                // P_12 P_34
                local_dgammads[{ eta_2body[i].s2, g->s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += eta_2body[i].me * g->me;
              }
          }
      }
  }

#pragma omp parallel sections
  {
#pragma omp section
    {
      for (const auto &local_map : thread_dfds)
        {
          for (const auto &[key, me] : local_map)
            {
              dfds_map[key] += me;
            }
        }
    }
#pragma omp section
    {
      for (const auto &local_map : thread_dgammads)
        {
          for (const auto &[key, me] : local_map)
            {
              dgammads_map[key] += me;
            }
        }
    }
  }
  // dfds_map = reduce_map(thread_dfds);
  // std::cout << "." << std::flush;
  // dgammads_map = reduce_map(thread_dgammads);
}

void const
Hamiltonian::set_flow_eq_print ()
{
  std::cout << "Constructing flow equations" << std::endl;
  auto start = std::chrono::high_resolution_clock::now ();

  std::cout << "\t" << "Clearing old flow eqs" << std::flush;
  dEds = 0;
  dfds_vec.clear ();
  dfds_map.clear ();
  dgammads_vec.clear ();
  dgammads_map.clear ();

  auto stop_clear = std::chrono::high_resolution_clock::now ();
  std::chrono::duration<double, std::milli> ms_double = stop_clear - start;
  std::cout << "..." << "\t" << " done (" << ms_double.count () << "ms)"
            << std::endl;

  std::cout << "\t" << "Creating maps" << std::flush;

  TwoBodyQA gamma_qa = TwoBodyQA (gamma, 1);
  std::cout << "." << std::flush;
  TwoBodyQA eta2_qa = TwoBodyQA (eta_2body, 1);

  std::unordered_map<OneBodyKey, double> f_map;
  f_map.reserve (f.size ());
  for (const auto &ff : f)
    {
      f_map[OneBodyKey (ff.s, ff.ss)] = ff.me;
    }

  std::cout << "." << std::flush;

  std::unordered_map<OneBodyKey, double> eta1_map;
  eta1_map.reserve (eta_1body.size ());
  for (const auto &e : eta_1body)
    {
      eta1_map[OneBodyKey (e.s, e.ss)] = e.me;
    }

  std::cout << "." << std::flush;

  const int f_size = f.size ();
  const int eta1_size = eta_1body.size ();
  const int gamma_size = gamma.size ();
  const int eta2_size = eta_2body.size ();

  int max_threads = omp_get_max_threads ();
  std::vector<std::unordered_map<OneBodyKey, double> > thread_dfds (
      max_threads);
  std::vector<std::unordered_map<TwoBodyKey, double> > thread_dgammads (
      max_threads);
  for (int i = 0; i < max_threads; i++)
    {
      thread_dfds[i].reserve (5000);
      thread_dgammads[i].reserve (5000);
    }
  auto stop_maps = std::chrono::high_resolution_clock::now ();
  ms_double = stop_maps - stop_clear;
  std::cout << "\t\t" << " done (" << ms_double.count () << "ms)" << std::endl;

  std::cout << "\t" << "Constructing new flow eqs" << std::flush;

#pragma omp parallel
  {
//////////////////////////////// dEds
#pragma omp for reduction(+ : dEds)
    for (int i = 0; i < f_size; ++i)
      {
        if (f[i].s > 3 && f[i].ss < 4 && eta1_map.find({f[i].ss, f[i].s}) != eta1_map.end())
          {
            dEds += f[i].me * eta1_map[OneBodyKey (f[i].ss, f[i].s)];
          }
        else if (f[i].s < 4 && f[i].s > 3 && eta1_map.find({f[i].ss, f[i].s}) != eta1_map.end())
          {
            dEds += -f[i].me * eta1_map[OneBodyKey (f[i].ss, f[i].s)];
          }
      }
#pragma omp for reduction(+ : dEds) nowait
    for (int i = 0; i < eta2_size; ++i)
      {
        if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4 && eta_2body[i].ss1 > 3
            && eta_2body[i].ss2 > 3)
          {
            if (gamma_qa.map.find ({ eta_2body[i].ss1, eta_2body[i].ss2,
                                     eta_2body[i].s1, eta_2body[i].s2 })
                != gamma_qa.map.end ())
              {
                dEds += 0.5 * eta_2body[i].me
                        * gamma_qa
                              .map[{ eta_2body[i].ss1, eta_2body[i].ss2,
                                     eta_2body[i].s1, eta_2body[i].s2 }]
                              ->me;
              }
          }
      }
#pragma omp single
    std::cout << "." << std::flush;
    ///////////////////////////////////////////////// dfds
    // first sum
    int tid = omp_get_thread_num ();
    std::unordered_map<OneBodyKey, double> &local_dfds = thread_dfds[tid];

#pragma omp for collapse(2)
    for (int i = 0; i < f_size; ++i)
      {
        for (int j = 0; j < eta1_size; ++j)
          {
            // first sum
            if (eta_1body[j].ss == f[i].s)
              {
                local_dfds[{ eta_1body[j].s, f[i].ss }]
                    += eta_1body[j].me * f[i].me;
                local_dfds[{ f[i].ss, eta_1body[j].s }]
                    += eta_1body[j].me * f[i].me;
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
                local_dfds[{ g->s2, g->ss2 }] += eta_1body[i].me * g->me;
              }
          }
        else if (eta_1body[i].s > 3 && eta_1body[i].ss < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_1body[i].ss, eta_1body[i].s))
              {
                local_dfds[{ g->s2, g->ss2 }] += -eta_1body[i].me * g->me;
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
                local_dfds[{ e->s2, e->ss2 }] += -f[i].me * e->me;
              }
          }
        else if (f[i].s > 3 && f[i].ss < 4)
          {
            for (const auto &e : eta2_qa.check_s1ss1 (f[i].ss, f[i].s))
              {
                local_dfds[{ e->s2, e->ss2 }] += f[i].me * e->me;
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
                local_dfds[{ eta_2body[i].s2, g->ss2 }]
                    += 0.5 * eta_2body[i].me * g->me;
                local_dfds[{ g->ss2, eta_2body[i].s2 }]
                    += 0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3
                 && eta_2body[i].s1 < 4)
          {
            for (const auto &g : gamma_qa.check_s1s2ss1 (
                     eta_2body[i].ss1, eta_2body[i].ss2, eta_2body[i].s1))
              {
                local_dfds[{ eta_2body[i].s2, g->ss2 }]
                    += -0.5 * eta_2body[i].me * g->me;
                local_dfds[{ g->ss2, eta_2body[i].s2 }]
                    += -0.5 * eta_2body[i].me * g->me;
              }
          }
      }
    ///////////////////////////// dGammads
    std::unordered_map<TwoBodyKey, double> &local_dgammads
        = thread_dgammads[tid];

// first sum
// gamma terms
#pragma omp for
    for (int i = 0; i < eta1_size; ++i)
      {
        for (const auto &g : gamma_qa.check_s1 (eta_1body[i].ss))
          {
            local_dgammads[{ eta_1body[i].s, g->s2, g->ss1, g->ss2 }]
                += eta_1body[i].me * g->me;
            local_dgammads[{ g->s2, eta_1body[i].s, g->ss1, g->ss2 }]
                += -eta_1body[i].me * g->me;
          }
        for (const auto &g : gamma_qa.check_ss1 (eta_1body[i].s))
          {
            local_dgammads[{ g->s1, g->s2, eta_1body[i].ss, g->ss2 }]
                += -eta_1body[i].me * g->me;
            local_dgammads[{ g->s1, g->s2, g->ss2, eta_1body[i].ss }]
                += eta_1body[i].me * g->me;
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
            local_dgammads[{ f[i].s, e->s2, e->ss1, e->ss2 }]
                += -f[i].me * e->me;
            local_dgammads[{ e->s2, f[i].s, e->ss1, e->ss2 }]
                += f[i].me * e->me;
          }
        for (const auto &e : eta2_qa.check_ss1 (f[i].s))
          {
            local_dgammads[{ e->s1, e->s2, f[i].ss, e->ss2 }]
                += f[i].me * e->me;
            local_dgammads[{ e->s1, e->s2, e->ss2, f[i].ss }]
                += -f[i].me * e->me;
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
                local_dgammads[{ eta_2body[i].s1, eta_2body[i].s2, g->ss1,
                                 g->ss2 }]
                    += -0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].s1 < 4 && eta_2body[i].s2 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_ss1ss2 (eta_2body[i].s1, eta_2body[i].s2))
              {
                local_dgammads[{ g->s1, g->s2, eta_2body[i].ss1,
                                 eta_2body[i].ss2 }]
                    += 0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].ss2 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1s2 (eta_2body[i].ss1, eta_2body[i].ss2))
              {
                local_dgammads[{ eta_2body[i].s1, eta_2body[i].s2, g->ss1,
                                 g->ss2 }]
                    += 0.5 * eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].s1 > 3 && eta_2body[i].s2 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_ss1ss2 (eta_2body[i].s1, eta_2body[i].s2))
              {
                local_dgammads[{ g->s1, g->s2, eta_2body[i].ss1,
                                 eta_2body[i].ss2 }]
                    += -0.5 * eta_2body[i].me * g->me;
              }
            // 3. sum
          }
        else if (eta_2body[i].ss1 < 4 && eta_2body[i].s1 > 3)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_2body[i].ss1, eta_2body[i].s1))
              {
                // 1
                local_dgammads[{ g->s2, eta_2body[i].s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += -eta_2body[i].me * g->me;
                // -P_12
                local_dgammads[{ eta_2body[i].s2, g->s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += eta_2body[i].me * g->me;
                // -P_34
                local_dgammads[{ g->s2, eta_2body[i].s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += eta_2body[i].me * g->me;
                // P_12 P_34
                local_dgammads[{ eta_2body[i].s2, g->s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += -eta_2body[i].me * g->me;
              }
          }
        else if (eta_2body[i].ss1 > 3 && eta_2body[i].s1 < 4)
          {
            for (const auto &g :
                 gamma_qa.check_s1ss1 (eta_2body[i].ss1, eta_2body[i].s1))
              {
                // 1
                local_dgammads[{ g->s2, eta_2body[i].s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += eta_2body[i].me * g->me;
                // -P_12
                local_dgammads[{ eta_2body[i].s2, g->s2, g->ss2,
                                 eta_2body[i].ss2 }]
                    += -eta_2body[i].me * g->me;
                // -P_34
                local_dgammads[{ g->s2, eta_2body[i].s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += -eta_2body[i].me * g->me;
                // P_12 P_34
                local_dgammads[{ eta_2body[i].s2, g->s2, eta_2body[i].ss2,
                                 g->ss2 }]
                    += eta_2body[i].me * g->me;
              }
          }
      }
  }
  auto stop_dgds = std::chrono::high_resolution_clock::now ();
  ms_double = stop_dgds - stop_maps;
  std::cout << "\t" << " done (" << ms_double.count () << "ms)" << std::endl;

  std::cout << "\t" << "Combining thread maps" << std::flush;
  std::cout << "." << std::flush;
#pragma omp parallel sections
  {
#pragma omp section
    {
      for (const auto &local_map : thread_dfds)
        {
          for (const auto &[key, me] : local_map)
            {
              dfds_map[key] += me;
            }
        }
    }
#pragma omp section
    {
      std::cout << "." << std::flush;
      for (const auto &local_map : thread_dgammads)
        {
          for (const auto &[key, me] : local_map)
            {
              dgammads_map[key] += me;
            }
        }
    }
  }
  // dfds_map = reduce_map(thread_dfds);
  // std::cout << "." << std::flush;
  // dgammads_map = reduce_map(thread_dgammads);

  std::cout << "." << std::flush;
  auto stop_comb = std::chrono::high_resolution_clock::now ();
  ms_double = stop_comb - stop_dgds;
  std::cout << "\t" << " done (" << ms_double.count () << "ms)" << std::endl;

  ms_double = stop_comb - start;
  std::cout << "\t" << "|-> done(" << ms_double.count () << "ms)" << std::endl;
}

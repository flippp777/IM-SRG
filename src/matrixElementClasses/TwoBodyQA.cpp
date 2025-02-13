#include "TwoBodyQA.hpp"
#include "TwoBody.hpp"
#include <cstdio>
#include <omp.h>
#include <vector>

TwoBodyQA::TwoBodyQA ()
    : index_s1 (), index_ss1 (), index_s1s2 (), index_ss1ss2 (),
      index_s1ss1 (), map ()
{
}

//   TwoBodyQA::~TwoBodyQA() {
//   map.clear();
//   index_s1.clear();
//   index_ss1.clear();
//   index_s1s2.clear();
//   index_s1ss1.clear();
//   index_ss1ss2.clear();
//   index_s1s2ss1.clear();
// }

TwoBodyQA::TwoBodyQA (const TwoBody element) { insert (element); }

TwoBodyQA::TwoBodyQA (const std::vector<TwoBody> &vec, int par)
{
  map.reserve (vec.size ());
  if (par == 0)
    {
      for (const auto &element : vec)
        {
          insert (element);
        }
    }
  else
    {
#pragma omp parallel sections
      {
#pragma omp section
        for (const auto &element : vec)
          {
            map[TwoBodyKey (element.s1, element.s2, element.ss1, element.ss2)]
                = element;
          }

#pragma omp section
        for (const auto &element : vec)
          {
            index_s1[element.s1].push_back (element);
          }
#pragma omp section
        for (const auto &element : vec)
          {
            index_ss1[element.ss1].push_back (element);
          }

#pragma omp section
        for (const auto &element : vec)
          {
            index_s1s2[element.s1][element.s2].push_back (element);
          }
#pragma omp section
        for (const auto &element : vec)
          {
            index_ss1ss2[element.ss1][element.ss2].push_back (element);
          }
#pragma omp section
        for (const auto &element : vec)
          {
            index_s1ss1[element.s1][element.ss1].push_back (element);
          }
#pragma omp section
        for (const auto &element : vec)
          {
            index_s1s2ss1[element.s1][element.s2][element.ss1].push_back (element);
          }
      }
    }
}

void
TwoBodyQA::insert (const TwoBody &element)
{
  map[{ element.s1, element.s2, element.ss1, element.ss2 }] = element;

  index_s1[element.s1].push_back (element);
  index_ss1[element.ss1].push_back (element);

  index_s1s2[element.s1][element.s2].push_back (element);
  index_ss1ss2[element.ss1][element.ss2].push_back (element);
  index_s1ss1[element.s1][element.ss1].push_back (element);
  index_s1s2ss1[element.s1][element.s2][element.ss1].push_back (element);
}

// check single indices
const std::vector<TwoBody>
TwoBodyQA::check_s1 (const int s1)
{
  if (index_s1.find (s1) != index_s1.end ())
    {
      return index_s1[s1];
    }
  return {};
}
const std::vector<TwoBody>
TwoBodyQA::check_ss1 (const int ss1)
{
  if (index_ss1.find (ss1) != index_ss1.end ())
    {
      return index_ss1[ss1];
    }
  return {};
}

// check two indices
const std::vector<TwoBody>
TwoBodyQA::check_s1s2 (const int s1, const int s2)
{
  if (index_s1s2.find (s1) != index_s1s2.end ()
      && index_s1s2[s1].find (s2) != index_s1s2[s1].end ())
    return index_s1s2[s1][s2];

  return {};
}
const std::vector<TwoBody>
TwoBodyQA::check_ss1ss2 (const int ss1, const int ss2)
{
  if (index_ss1ss2.find (ss1) != index_ss1ss2.end ()
      && index_ss1ss2[ss1].find (ss2) != index_ss1ss2[ss1].end ())
    return index_ss1ss2[ss1][ss2];

  return {};
}
const std::vector<TwoBody>
TwoBodyQA::check_s1ss1 (const int s1, const int ss1)
{
  if (index_s1ss1.find (s1) != index_s1ss1.end ()
      && index_s1ss1[s1].find (ss1) != index_s1ss1[s1].end ())
    return index_s1ss1[s1][ss1];

  return {};
}
const std::vector<TwoBody>
TwoBodyQA::check_s1s2ss1 (const int s1, const int s2, const int ss1)
{
  if (index_s1s2ss1.find (s1) != index_s1s2ss1.end ()
      && index_s1s2ss1[s1].find (s2) != index_s1s2ss1[s1].end ()
      && index_s1s2ss1[s1][s2].find (ss1) != index_s1s2ss1[s1][s2].end ())
    return index_s1s2ss1[s1][s2][ss1];

  return {};
}

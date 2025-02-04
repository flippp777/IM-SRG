#include "TwoBodyQA.hpp"
#include "TwoBody.hpp"
#include <cstdio>
#include <omp.h>
#include <vector>

TwoBodyQA::TwoBodyQA()
    : map(), index_s1(), index_ss1(), index_s1s2(), index_ss1ss2(),
      index_s1ss1() {}

TwoBodyQA::TwoBodyQA(const TwoBody element) { insert(element); }

TwoBodyQA::TwoBodyQA(const std::vector<TwoBody> &vec, int par) {
  map.reserve(vec.size());
  if (par == 0) {
    for (const auto &element : vec) {
      insert(element);
    }
  } else {
#pragma omp parallel sections
    {
#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        map[TwoBodyKey(element.s1, element.s2, element.ss1, element.ss2)] = ptr;
      }

#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        index_s1[element.s1].push_back(ptr);
      }
#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        index_ss1[element.ss1].push_back(ptr);
      }

#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        index_s1s2[element.s1][element.s2].push_back(ptr);
      }
#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        index_ss1ss2[element.ss1][element.ss2].push_back(ptr);
      }
#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        index_s1ss1[element.s1][element.ss1].push_back(ptr);
      }
#pragma omp section
      for (const auto &element : vec) {
        const TwoBody *ptr = &element;
        index_s1s2ss1[element.s1][element.s2][element.ss1].push_back(ptr);
      }
    }
  }
}

void TwoBodyQA::insert(const TwoBody &element) {
  const TwoBody *ptr = &element;
  map[{element.s1, element.s2, element.ss1, element.ss2}] = ptr;

  index_s1[element.s1].push_back(ptr);
  index_ss1[element.ss1].push_back(ptr);

  index_s1s2[element.s1][element.s2].push_back(ptr);
  index_ss1ss2[element.ss1][element.ss2].push_back(ptr);
  index_s1ss1[element.s1][element.ss1].push_back(ptr);
  index_s1s2ss1[element.s1][element.s2][element.ss1].push_back(ptr);
}

// check single indices
const std::vector<const TwoBody *> TwoBodyQA::check_s1(const int s1) {
  return index_s1[s1];
}
const std::vector<const TwoBody *> TwoBodyQA::check_ss1(const int ss1) {
  return index_ss1[ss1];
}

// check two indices
const std::vector<const TwoBody *> TwoBodyQA::check_s1s2(const int s1,
                                                         const int s2) {
  return index_s1s2[s1][s2];
}
const std::vector<const TwoBody *> TwoBodyQA::check_ss1ss2(const int ss1,
                                                           const int ss2) {
  return index_ss1ss2[ss1][ss2];
}
const std::vector<const TwoBody *> TwoBodyQA::check_s1ss1(const int s1,
                                                          const int ss1) {
  return index_s1ss1[s1][ss1];
}
const std::vector<const TwoBody *>
TwoBodyQA::check_s1s2ss1(const int s1, const int s2, const int ss1) {
  return index_s1s2ss1[s1][s2][ss1];
}

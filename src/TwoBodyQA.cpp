#include "TwoBodyQA.hpp"
#include <cstdio>

TwoBodyQA::TwoBodyQA()
    : gamma(), gamma_map(), index_s1(), index_ss1(), index_s1s2(),
      index_ss1ss2(), index_s1ss1() {}

TwoBodyQA::TwoBodyQA(const TwoBody element) { insert(element); }

TwoBodyQA::TwoBodyQA(const std::vector<TwoBody> vec) {
  gamma.reserve(vec.size());
  gamma_map.reserve(vec.size());
  for (const auto &element : vec) {
    insert(element);
  }
}

void TwoBodyQA::insert(TwoBody element) {
  TwoBody *ptr = &element;
  gamma.push_back(ptr);
  gamma_map[TwoBodyKey(element.s1, element.s2, element.ss1, element.ss2)] =
      &element.me;

  index_s1[element.s1].push_back(ptr);
  index_ss1[element.ss1].push_back(ptr);

  index_s1s2[element.s1][element.s2].push_back(ptr);
  index_ss1ss2[element.ss1][element.ss2].push_back(ptr);
  index_s1ss1[element.s1][element.ss1].push_back(ptr);
  index_s1s2ss1[element.s1][element.s2][element.ss1].push_back(ptr);
}

// check single indices
std::vector<TwoBody *> TwoBodyQA::check_s1(const int s1) {
  return index_s1[s1];
}
std::vector<TwoBody *> TwoBodyQA::check_ss1(const int ss1) {
  return index_ss1[ss1];
}

// check two indices
std::vector<TwoBody *> TwoBodyQA::check_s1s2(const int s1, const int s2) {
  return index_s1s2[s1][s2];
}
std::vector<TwoBody *> TwoBodyQA::check_ss1ss2(const int ss1, const int ss2) {
  return index_ss1ss2[ss1][ss2];
}
std::vector<TwoBody *> TwoBodyQA::check_s1ss1(const int s1, const int ss1) {
  return index_s1ss1[s1][ss1];
}
std::vector<TwoBody *> TwoBodyQA::check_s1s2ss1(const int s1, const int s2,
                                                const int ss1) {
  return index_s1s2ss1[s1][s2][ss1];
}
double * TwoBodyQA::check_s1s2ss1ss2(const int s1, const int s2, const int ss1,
                                   const int ss2) {
  if (gamma_map.count({s1, s2, ss1, ss2}) == 1) {
    return gamma_map[{s1, s2, ss1, ss2}];
  } else if (gamma_map.count({s1, s2, ss1, ss2}) == 0) {
    printf("set_flow_eq: No gamma");
    return 0;
  } else {
    printf("set_flow_eq: More than one Gamma");
    return 0;
  }
}

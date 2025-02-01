#include "TwoBody.hpp"
#include <algorithm>
#include <unordered_map>
#include <vector>

TwoBody::TwoBody() : s1(), s2(), ss1(), ss2(), me() {}

TwoBody::TwoBody(int s1, int s2, int ss1, int ss2, double me)
    : s1(s1), s2(s2), ss1(ss1), ss2(ss2), me(me) {}



















///////////////////////////////////////////////
// unordered_map implementation
///////////////////////////////////////////////


std::unordered_map<TwoBodyKey, double> vec_to_map(const std::vector<TwoBody>& vec) {
  std::unordered_map<TwoBodyKey, double> map;
  map.reserve(vec.size());

  for (const auto& el : vec) {
    map[TwoBodyKey(el.s1, el.s2, el.ss1, el.ss2)] += el.me;
  }
  
  return map;
}

std::vector<TwoBody> map_to_vec(const std::unordered_map<TwoBodyKey, double>& map) {
  std::vector<TwoBody> vec;
  vec.reserve(map.size());

  for (const auto& [key, me] : map) {
    vec.push_back(TwoBody(key.s1, key.s2, key.ss1, key.ss2, me));
  }

  return vec;
}

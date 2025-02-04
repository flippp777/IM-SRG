#include "TwoBody.hpp"
#include <algorithm>
#include <unordered_map>
#include <vector>

TwoBody::TwoBody() : s1(), s2(), ss1(), ss2(), me() {}

TwoBody::TwoBody(int s1, int s2, int ss1, int ss2, double me)
    : s1(s1), s2(s2), ss1(ss1), ss2(ss2), me(me) {}

std::vector<TwoBody> operator+(const std::vector<TwoBody> &v,
                               const TwoBody &a) {
  std::vector<TwoBody> new_vec = v;

  for (int i = 0; i < new_vec.size(); ++i) {
    if (new_vec[i].s1 == a.s1 && new_vec[i].s2 == a.s2 &&
        new_vec[i].ss1 == a.ss1 && new_vec[i].ss2 == a.ss2) {
      new_vec[i].me = new_vec[i].me + a.me;
      return new_vec;
    }
  }
  new_vec.push_back(a);
  return new_vec;
}
std::vector<TwoBody> operator+(const std::vector<TwoBody> &v1,
                               const std::vector<TwoBody> &v2) {
  std::vector<TwoBody> new_vec = v1;
  const int v2_size = v2.size();

  for (int j = 0; j < v2_size; ++j) {
    new_vec = new_vec + v2[j];
  }
  return new_vec;
}

void sort(std::vector<TwoBody> &v) {
  std::sort(v.begin(), v.end(), [](TwoBody a, TwoBody b) {
    if (a.s1 == b.s1 && a.s2 == b.s2 && a.ss1 == b.ss1) {
      return a.ss2 < b.ss2;
    } else if (a.s1 == b.s1 && a.s2 == b.s2) {
      return a.ss1 < b.ss1;
    } else if (a.s1 == b.s1) {
      return a.s2 < b.s2;
    } else {
      return a.s1 < b.s2;
    }
  });
}
///////////////////////////////////////////////
// unordered_map implementation
///////////////////////////////////////////////

std::unordered_map<TwoBodyKey, double>
vec_to_map(const std::vector<TwoBody> &vec) {
  std::unordered_map<TwoBodyKey, double> map;
  map.reserve(vec.size());

  for (const auto &el : vec) {
    map[TwoBodyKey(el.s1, el.s2, el.ss1, el.ss2)] += el.me;
  }

  return map;
}

std::vector<TwoBody>
map_to_vec(const std::unordered_map<TwoBodyKey, double> &map) {
  std::vector<TwoBody> vec;
  vec.reserve(map.size());

  for (const auto &[key, me] : map) {
    vec.push_back(TwoBody(key.s1, key.s2, key.ss1, key.ss2, me));
  }

  return vec;
}

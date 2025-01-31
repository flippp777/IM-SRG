#include "TwoBody.hpp"
#include <algorithm>
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

void sort(std::vector<TwoBody> &v) {
  std::sort(v.begin(), v.end(), [](TwoBody a, TwoBody b) {
    if (a.s1 == b.s1) {
      if (a.s2 == b.s2) {
        if (a.ss1 == b.ss1) {
          return a.ss2 < b.ss2;
        } else {
          return a.ss1 < b.ss1;
        }
      } else {
        return a.s2 < b.s2;
      }
    } else {
      return a.s1 < b.s1;
    }
  });
}

///////////////////////////////////////////////
// unordered_map implementation
///////////////////////////////////////////////

template <> struct std::hash<TwoBodyKey> {
  size_t operator()(const TwoBodyKey &k) const {
    size_t seed = 0;

    // Robust hash combining
    auto hash_combine = [](size_t &seed, auto val) {
      std::hash<std::decay_t<decltype(val)>> hasher;
      seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    hash_combine(seed, k.s1);
    hash_combine(seed, k.s2);
    hash_combine(seed, k.ss1);
    hash_combine(seed, k.ss2);

    return seed;
  }
};

#ifndef TWOBODY
#define TWOBODY

#include <functional>
#include <unordered_map>
#include <vector>

class TwoBody {
public:
  int s1, s2, ss1, ss2;
  double me;

  TwoBody();
  TwoBody(int s1, int s2, int ss1, int ss2, double me);
};

std::vector<TwoBody> operator+(const std::vector<TwoBody> &v, const TwoBody &a);

void sort(std::vector<TwoBody> &v);

///////////////////////////////////////////////
// unordered_map implementation
///////////////////////////////////////////////

struct TwoBodyKey {
  int s1;
  int s2;
  int ss1;
  int ss2;

  TwoBodyKey() {}
  TwoBodyKey(int s1, int s2, int ss1, int ss2)
      : s1(s1), s2(s2), ss1(ss1), ss2(ss2) {}

  bool operator==(const TwoBodyKey &other) const {
    return std::tie(s1, s2, ss1, ss2) ==
           std::tie(other.s1, other.s2, other.ss1, other.ss2);
  }
};

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

std::vector<TwoBody> map_to_vec(const std::unordered_map<TwoBodyKey, double>& map);

std::unordered_map<TwoBodyKey, double> vec_to_map(const std::vector<TwoBody>& vec);

#endif // !TWOBODY

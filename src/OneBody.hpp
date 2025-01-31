#ifndef ONEBODY
#define ONEBODY

#include <unordered_map>
#include <vector>
struct OneBody {
public:
  int s, ss;
  double me;

  OneBody();
  OneBody(int s, int ss, double me);
};

std::vector<OneBody> operator+(const std::vector<OneBody> &v,
                               const OneBody &a);
std::vector<OneBody> operator+(const std::vector<OneBody> &v1,
                               const std::vector<OneBody> &v2);

void sort(std::vector<OneBody> &v);

///////////////////////////////////////////////
// unordered_map implementation
///////////////////////////////////////////////

struct OneBodyKey {
  int s;
  int ss;

  OneBodyKey() {}
  OneBodyKey(int s, int ss)
      : s(s), ss(ss) {}

  bool operator==(const OneBodyKey &other) const {
    return std::tie(s, ss) ==
           std::tie(other.s, other.ss);
  }
};

template <> struct std::hash<OneBodyKey> {
  size_t operator()(const OneBodyKey &k) const {
    size_t seed = 0;

    // Robust hash combining
    auto hash_combine = [](size_t &seed, auto val) {
      std::hash<std::decay_t<decltype(val)>> hasher;
      seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };

    hash_combine(seed, k.s);
    hash_combine(seed, k.ss);

    return seed;
  }
};
std::unordered_map<OneBodyKey, double> vec_to_map(const std::vector<OneBody>& vec);
std::vector<OneBody> map_to_vex(const std::unordered_map<OneBodyKey, double>& map);

#endif // !ONEBODY

#ifndef TWOBODY
#define TWOBODY

#include <functional>
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

template <> struct std::hash<TwoBodyKey>;

#endif // !TWOBODY

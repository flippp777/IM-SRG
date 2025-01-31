#ifndef ONEBODY
#define ONEBODY

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

template <> struct std::hash<OneBodyKey>;

#endif // !ONEBODY

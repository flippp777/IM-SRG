#ifndef TWOBODY
#define TWOBODY

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

#endif // !TWOBODY

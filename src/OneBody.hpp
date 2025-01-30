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

void sort(std::vector<OneBody> &v);

#endif // !ONEBODY

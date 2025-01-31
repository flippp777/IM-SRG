#include "OneBody.hpp"
#include <algorithm>
#include <vector>

OneBody::OneBody() : s(), ss(), me() {}

OneBody::OneBody(int s, int ss, double me) : s(s), ss(ss), me(me) {}

std::vector<OneBody> operator+(const std::vector<OneBody> &v,
                               const OneBody &a) {
  std::vector<OneBody> new_vec = v;
  for (int i = 0; i < new_vec.size(); ++i) {
    if (new_vec[i].s == a.s && new_vec[i].ss == a.ss) {
      new_vec[i].me = new_vec[i].me + a.me;
      return new_vec;
    }
  }
  new_vec.push_back(a);
  return new_vec;
}
std::vector<OneBody> operator+(const std::vector<OneBody> &v1,
                               const std::vector<OneBody> &v2) {
  std::vector<OneBody> new_vec = v1;
  const int v2_size = v2.size();

  for (int j = 0; j < v2_size; ++j) {
    new_vec = new_vec + v2[j];
  }
  return new_vec;
}

void sort(std::vector<OneBody> &v) {
  std::sort(v.begin(), v.end(), [](OneBody a, OneBody b) {
    if (a.s == b.s) {
      return a.ss < b.ss;
    } else {
      return a.s < b.s;
    }
  });
}

#include "TwoBody.hpp"

TwoBody::TwoBody() : s1(), s2(), ss1(), ss2(), me() {}

TwoBody::TwoBody(int s1, int s2, int ss1, int ss2, double me)
    : s1(s1), s2(s2), ss1(ss1), ss2(ss2), me(me) {}

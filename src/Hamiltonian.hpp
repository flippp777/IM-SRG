#include <string>
#include <vector>

struct OneBody {
  int s, ss;
  double me;
};

struct TwoBody {
  int s1, s2, ss1, ss2;
  double me;
};

class Hamiltonian {
private:
  double E;
  std::vector<OneBody> f;
  std::vector<TwoBody> Gamma;

public:
  Hamiltonian(std::string filepath);

  const std::vector<TwoBody> &get_two_body() const {
    return Gamma;
  }
  const std::vector<OneBody> &get_one_body() const {
    return f;
  }
  const double &get_zero_body() const {
    return E;
  }

  void f0();
  void Gamma0();
};

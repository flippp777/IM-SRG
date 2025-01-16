#include <string>
#include <vector>

struct MatrixElement {
  int s1, s2, ss1, ss2;
  double me;
};

class MatrixElements {
private:
  std::vector<int> states;
  std::vector<MatrixElement> elements;

public:
  double E;
  std::vector<std::vector<double>> f;
  std::vector<std::vector<std::vector<std::vector<double>>>> Gamma;
  int sqrt_matrix_size = 139;

  MatrixElements(std::string filename);
  const std::vector<int> &getStates() const { return states; }
  const std::vector<MatrixElement> &getElements() const {
    return elements;
  }

  std::vector<MatrixElement> truncateElementsS(int Smax) const;
  void setTruncatedElementsS(int Nmax);

  void E0();
  void f0();
  void Gamma0();
};

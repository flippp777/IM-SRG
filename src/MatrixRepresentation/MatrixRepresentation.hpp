#ifndef MATRIX
#define MATRIX

#include "../Hamiltonian/Hamiltonian.hpp"
#include "../matrixElementClasses/OneBody.hpp"
#include "../matrixElementClasses/TwoBody.hpp"
#include <unordered_map>

class MatrixRepresentation
{
private:
  std::unordered_map<OneBodyKey, double> f;
  std::unordered_map<TwoBodyKey, double> G;

  std::unordered_map<OneBodyKey, int> index_f;
  std::unordered_map<OneBodyKey, int> index_G;

public:
  MatrixRepresentation (const Hamiltonian &H);

  const void export_one_body ();
  const void export_two_body ();
};

#endif // !MATRIX

#ifndef MATRIX
#define MATRIX

#include "../Hamiltonian/Hamiltonian.hpp"
#include "../matrixElementClasses/OneBody.hpp"
#include "../matrixElementClasses/TwoBody.hpp"

class MatrixRepresentation
{
private:
  int num_states;

  std::vector<OneBody> f;
  std::vector<TwoBody> Gamma;

  std::unordered_map<OneBodyKey, int> index_two_body;

public:
  MatrixRepresentation (Hamiltonian &H);

  const void set_one_body(std::vector<OneBody> new_f) { f = new_f; }
  const void set_two_body(std::vector<TwoBody> new_Gamma) { Gamma = new_Gamma; }

  const void export_one_body (std::string name);
  const void export_two_body (std::string name);
};

#endif // !MATRIX

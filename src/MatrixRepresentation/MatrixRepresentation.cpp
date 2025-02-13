#include "MatrixRepresentation.hpp"
#include <fstream>

MatrixRepresentation::MatrixRepresentation(Hamiltonian &H) {
  num_states = H.get_num_states();

  f = H.get_one_body();
  Gamma = H.get_two_body();

  int idx = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = i; i < 4; ++i) {
      index_two_body.insert({{i, j}, idx});
      idx++;
    }
  }
  for (int i = 4; i < num_states; ++i) {
    for (int j = i; i < num_states; ++i) {
      index_two_body.insert({{i, j}, idx});
      idx++;
    }
  }
}

const void MatrixRepresentation::export_one_body(std::string name) {
  std::ofstream f_file(name);
  for (const auto &ff : f) {
    f_file << ff.s << "\t" << ff.ss << "\t" << ff.me << "\n";
  }
  f_file.close();
}

const void MatrixRepresentation::export_two_body(std::string name) {
  std::ofstream Gamma_file(name);
  for (const auto &gg : Gamma) {
    if (index_two_body.find({gg.s1, gg.s2}) != index_two_body.end() &&
        index_two_body.find({gg.ss1, gg.ss2}) != index_two_body.end()) {
      Gamma_file << index_two_body.at({gg.s1, gg.s2}) << "\t"
                 << index_two_body.at({gg.ss1, gg.ss2}) << "\t" << gg.me
                 << "\n";
    }
  }
  Gamma_file.close();
}

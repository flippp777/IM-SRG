#include <iostream>
#include <vector>
#include <string>

#include "MatrixElements.hpp"

int main() {
    MatrixElements H0 = MatrixElements("matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m");

    H0.setTruncatedElementsS(19);

    H0.E0();
    std::cout << "E0 = " << H0.E << std::endl;

    H0.f0();
    for (int i = 0; i < H0.f.size(); i++) {
        for (int j = 0; j < H0.f[i].size(); j++) {
            std::cout << H0.f[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    H0.Gamma0();
    for (int i = 0; i < H0.Gamma.size(); i++) {
        for (int j = 0; j < H0.Gamma[i].size(); j++){
            std::cout << H0.Gamma[i][j][0][0] << " ";
        }
        std::cout << std::endl;
    }


    return 0;
}




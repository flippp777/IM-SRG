#include <iostream>
#include <vector>
#include <string>

#include "MatrixElements.hpp"

int main() {
    MatrixElements H0 = MatrixElements("matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m");
    H0.E0();
    std::cout << "E0 = " << H0.E << std::endl;

    H0.f0();
    for (const auto& row : H0.f) {
        for (const auto& element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}




#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

#include "../include/E.hpp"

int main() {
    std::string path = "../matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m";
    double E_1 = E(path);
    std::cout << E_1 << std::endl;
    return 0;
}




#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

#include "../include/E.hpp"
#include "../include/f.hpp"

int main() {
    std::ifstream file("matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m");
    if (!file) {
        std::cerr << "error: could not open file" << std::endl;
        return 1;
    }

    std::vector<std::vector<double>> states;
    std::vector<std::vector<double>> elements;
    std::string line;

    int i = 0;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double num;
        if (2 < i && i < 142) {
            while (iss >> num) {
                row.push_back(num);
            }
            states.push_back(row);
        }
        else if (145 < i) {
            while (iss >> num) {
                row.push_back(num);
            }
            elements.push_back(row);
        }
        else {}
        i++;
    }
    file.close();

    double E = E0(elements);
    std::cout << E << std::endl;

    std::vector<std::vector<double>> f = f0(elements);
    std::cout << f.size() << "x" <<  f[0].size() << std::endl;


    return 0;
}




#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

int main() {
    std::ifstream file("../matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m");
    if (!file) {
        std::cerr << "error: could not open file" << std::endl;
        return 1;
    }

    std::vector<std::vector<int>> states;
    std::string line;

    while (std::getline(file, line)) {
        std::vector<int> row;
        std::istringstream iss(line);
        int num;
        while (iss >> num) {
            row.push_back(num);
        }
        states.push_back(row);
    }

    file.close();

    for (const auto& row : states) {
        for (const auto& num : row) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}




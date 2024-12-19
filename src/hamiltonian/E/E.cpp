#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

int main() {
    std::ifstream file("../../../matrix_elements/hint_He4_chi2bM3A_srg0800_eMax04_nMax02_lMax04_hwHO020.me2m");
    if (!file) {
        std::cerr << "error: could not open file" << std::endl;
        return 1;
    }

    std::vector<std::vector<float>> states;
    std::vector<std::vector<float>> V;
    std::string line;

    int i = 0;
    while (std::getline(file, line)) {
        std::vector<float> row;
        std::istringstream iss(line);
        float num;
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
            V.push_back(row);
        }
        else {}
        i++;
    }

    file.close();


    
    std::cout << states[0][0] << std::endl;
    std::cout << states[0][0] << std::endl;
    std::cout << V[0][0] << std::endl;
    std::cout << V[0][0] << std::endl;
    
    // for (const auto& row : states) {
    //     for (const auto& num : row) {
    //         std::cout << num << " ";
    //     }
    //     std::cout << std::endl;
    // }



    return 0;
}
#include "MatrixElements.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

MatrixElements::MatrixElements(std::string filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "error: could not open file" << std::endl;
        return;
    }

    std::string line;
    int i = 0;
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::istringstream iss(line);
        double num;
        if (i > 2 && i < 142) {
            while (iss >> num) {
                row.push_back(num);
            }
            states.push_back(static_cast<int>(row[1]));
        } else if (i > 145) {
            MatrixElement element;
            iss >> element.s1 >> element.s2 >> element.ss1 >> element.ss2 >> element.me;
            elements.push_back(element);
        }
        i++;
    }
    file.close();
}

std::vector<MatrixElement> MatrixElements::truncateElements(int Nmax) const {
    std::vector<MatrixElement> truncated_elements;
    for (const auto& element : elements) {
        if (states[element.s1] + states[element.s2] <= Nmax && states[element.ss1] + states[element.ss2] <= Nmax) {
            truncated_elements.push_back(element);
        }
    }
    return truncated_elements;
}

void MatrixElements::setTruncatedElements(int Nmax) {
    elements = truncateElements(Nmax);
}

void MatrixElements::E0() {
    double E0 = 0.0;
    for (const auto& element : elements) {
        E0 += element.me;
    }
    E = E0;
}
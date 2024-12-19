#include "../include/f.hpp"

#include <vector>

using namespace std;


double f_ij(const int i, int j, const vector<vector<double>> elements) {
    double f_ij = 0; 
    
    for (auto & element : elements) {
        if (element[0] == i && element[1] < 4 && element[2] == j && element[3] == element[1]) {
            f_ij += element[4];
        }
        else if (element[1] == i && element[0] < 4 && element[3] == j && element[2] == element[1]) {
            f_ij += element[4];
        }
    }
    
    return f_ij;
}

vector<vector<double>> f0(vector<vector<double>> elements) {
    vector<vector<double>> f;
    vector<double> row;

    int N = 0;

    for (int i = 0; i < elements.size(); i++) {
        for (int j = 0; j < elements.size(); j++) {
            row.push_back(f_ij(i, j, elements));
        }
        f.push_back(row);
    }
    return f;
} 
#include "../include/E.hpp"
#include <vector>
#include <numeric>


double E_0(std::vector<std::vector<double>> states, std::vector<std::vector<double>> elements) {
   double E = 0;

   for (auto & element : elements) {
    if (element[0] < 4 && element[1] < 4 && element[2] < 4 && element[3] < 4) {
        E += element[4];
    }
   }
    return E;
}
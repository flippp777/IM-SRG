#ifndef REDUCELIST
#define REDUCELIST 
#include "../matrixElementClasses/OneBody.hpp"
#include "../matrixElementClasses/TwoBody.hpp"
#include <unordered_map>
#include <vector>

std::unordered_map<OneBodyKey, double> reduce_map(std::vector<std::unordered_map<OneBodyKey, double>> &vec);
std::unordered_map<TwoBodyKey, double> reduce_map(std::vector<std::unordered_map<TwoBodyKey, double>> &vec);

std::vector<OneBody> sum_equal_elements(std::vector<OneBody> &vec);
std::vector<OneBody> reduce_vec(std::vector<std::vector<OneBody>> &vec);

std::vector<TwoBody> sum_equal_elements(std::vector<TwoBody> &vec);
std::vector<TwoBody> reduce_vec(std::vector<std::vector<TwoBody>> &vec);

#endif // !REDUCELIST

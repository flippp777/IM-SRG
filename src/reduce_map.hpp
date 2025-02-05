#ifndef REDUCEMAP
#define REDUCEMAP 
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include <unordered_map>

std::unordered_map<OneBodyKey, double> reduce_map(const std::vector<std::unordered_map<OneBodyKey, double>> &vec);

std::unordered_map<TwoBodyKey, double> reduce_map(const std::vector<std::unordered_map<TwoBodyKey, double>> &vec);

#endif // !REDUCEMAP

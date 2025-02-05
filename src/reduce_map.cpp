#include "reduce_map.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include <unordered_map>
#include <vector>

std::unordered_map<OneBodyKey, double> reduce_map(
     std::vector<std::unordered_map<OneBodyKey, double>> &old_maps) {
  while (old_maps.size() > 1) {
    std::vector<std::unordered_map<OneBodyKey, double>> new_maps;

#pragma omp parallel
    {
      std::vector<std::unordered_map<OneBodyKey, double>> local_maps;

#pragma omp for nowait
      for (int i = 0; i < old_maps.size() / 2; ++i) {
        std::unordered_map<OneBodyKey, double> new_map = old_maps[2 * i];
        for (const auto &[key, el] : old_maps[2 * i + 1]) {
          new_map[key] += el;
        }
        local_maps.push_back(new_map);
      }
      #pragma omp critical
      new_maps.insert(new_maps.end(), local_maps.begin(), local_maps.end());
    }

    if (old_maps.size() % 2 == 1 ){
      new_maps.push_back(old_maps.back());
    }

    old_maps = std::move(new_maps);
  }
  return old_maps[0];
}

std::unordered_map<TwoBodyKey, double> reduce_map(
     std::vector<std::unordered_map<TwoBodyKey, double>> &old_maps) {
  while (old_maps.size() > 1) {
    std::vector<std::unordered_map<TwoBodyKey, double>> new_maps;

#pragma omp parallel
    {
      std::vector<std::unordered_map<TwoBodyKey, double>> local_maps;

#pragma omp for nowait
      for (int i = 0; i < old_maps.size() / 2; ++i) {
        std::unordered_map<TwoBodyKey, double> new_map = old_maps[2 * i];
        for (const auto &[key, el] : old_maps[2 * i + 1]) {
          new_map[key] += el;
        }
        local_maps.push_back(new_map);
      }
      #pragma omp critical
      new_maps.insert(new_maps.end(), local_maps.begin(), local_maps.end());
    }

    if (old_maps.size() % 2 == 1 ){
      new_maps.push_back(old_maps.back());
    }

    old_maps = std::move(new_maps);
  }
  return old_maps[0];
}

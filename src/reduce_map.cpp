#include "reduce_map.hpp"
#include "OneBody.hpp"
#include "TwoBody.hpp"
#include <omp.h>
#include <unordered_map>
#include <vector>

std::unordered_map<OneBodyKey, double> reduce_map(
     std::vector<std::unordered_map<OneBodyKey, double>> &old_maps) {

  while (old_maps.size() > 1) {
    std::vector<std::unordered_map<OneBodyKey, double>> new_maps((int)old_maps.size()/2);

#pragma omp parallel num_threads((int)old_maps.size()/2)
    {
     std::unordered_map<OneBodyKey, double> local_map;

#pragma omp for nowait
      for (int i = 0; i < old_maps.size() / 2; ++i) {
        local_map = old_maps[2 * i];
        for (const auto &[key, el] : old_maps[2 * i + 1]) {
          local_map[key] += el;
        }
      }
      new_maps[omp_get_thread_num()] = local_map;
    }

    if (old_maps.size() % 2 == 1 ){
      new_maps.push_back(std::move(old_maps.back()));
    }

    old_maps = std::move(new_maps);
  }
  return old_maps[0];
}

std::unordered_map<TwoBodyKey, double> reduce_map(
     std::vector<std::unordered_map<TwoBodyKey, double>> &old_maps) {

  while (old_maps.size() > 1) {
    std::vector<std::unordered_map<TwoBodyKey, double>> new_maps((int)old_maps.size()/2);

#pragma omp parallel num_threads((int)old_maps.size()/2)
    {
      std::unordered_map<TwoBodyKey, double> local_map;

#pragma omp for nowait
      for (int i = 0; i < old_maps.size() / 2; ++i) {
        local_map = old_maps[2 * i];
        for (const auto &[key, el] : old_maps[2 * i + 1]) {
          local_map[key] += el;
        }
      }
      new_maps[omp_get_thread_num()] = local_map;
    }

    if (old_maps.size() % 2 == 1 ){
      new_maps.push_back(std::move(old_maps.back()));
    }

    old_maps = std::move(new_maps);
  }
  return old_maps[0];
}

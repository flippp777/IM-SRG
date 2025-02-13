#include "reduce_lists.hpp"
#include "../matrixElementClasses/OneBody.hpp"
#include "../matrixElementClasses/TwoBody.hpp"
#include <iostream>
#include <omp.h>
#include <unordered_map>
#include <vector>

std::unordered_map<OneBodyKey, double>
reduce_map(std::vector<std::unordered_map<OneBodyKey, double>> &old_maps) {

  while (old_maps.size() > 1) {
    std::vector<std::unordered_map<OneBodyKey, double>> new_maps(
        (int)old_maps.size() / 2);

#pragma omp parallel num_threads((int)old_maps.size() / 2)
    {
      std::unordered_map<OneBodyKey, double> local_map;

#pragma omp for nowait
      for (size_t i = 0; i < old_maps.size() / 2; ++i) {
        local_map = old_maps[2 * i];
        for (const auto &[key, el] : old_maps[2 * i + 1]) {
          local_map[key] += el;
        }
      }
      new_maps[omp_get_thread_num()] = local_map;
    }

    if (old_maps.size() % 2 == 1) {
      new_maps.push_back(std::move(old_maps.back()));
    }

    old_maps = std::move(new_maps);
  }
  return old_maps[0];
}

std::unordered_map<TwoBodyKey, double>
reduce_map(std::vector<std::unordered_map<TwoBodyKey, double>> &old_maps) {

  while (old_maps.size() > 1) {
    std::vector<std::unordered_map<TwoBodyKey, double>> new_maps(
        (int)old_maps.size() / 2);

#pragma omp parallel num_threads((int)old_maps.size() / 2)
    {
      std::unordered_map<TwoBodyKey, double> local_map;

#pragma omp for nowait
      for (size_t i = 0; i < old_maps.size() / 2; ++i) {
        local_map = old_maps[2 * i];
        for (const auto &[key, el] : old_maps[2 * i + 1]) {
          local_map[key] += el;
        }
      }
      new_maps[omp_get_thread_num()] = local_map;
    }

    if (old_maps.size() % 2 == 1) {
      new_maps.push_back(std::move(old_maps.back()));
    }

    old_maps = std::move(new_maps);
  }
  return old_maps[0];
}

std::vector<OneBody> sum_equal_elements(std::vector<OneBody> &vec) {
  sort(vec);
  int vec_size = vec.size();
  std::vector<OneBody> new_vec;
  new_vec.reserve(vec_size);

  for (int i = 0; i < vec_size; ++i) {
    double sum = vec[i].me;

    while (i + 1 < vec_size && vec[i].s == vec[i + 1].s && vec[i].ss == vec[i + 1].ss) {
      sum += vec[i + 1].me;
      ++i;
    }
    new_vec.push_back(OneBody(vec[i].s, vec[i].ss, sum));
  }
  return new_vec;
}

std::vector<OneBody> reduce_vec(std::vector<std::vector<OneBody>> &vecs) {

  while (vecs.size() > 1) {
    std::vector<std::vector<OneBody>> new_vecs((int)vecs.size() / 2);

#pragma omp parallel for
    for (size_t i = 0; i < vecs.size() / 2; ++i) {
      std::vector<OneBody> local_vec = vecs[i * 2];
      local_vec.insert(local_vec.end(), vecs[2 * i + 1].begin(),
                       vecs[2 * i + 1].end());
      local_vec = sum_equal_elements(local_vec);
      new_vecs[i] = local_vec;
    }

    if (vecs.size() % 2 == 1) {
      new_vecs.push_back(std::move(vecs.back()));
    }

    vecs = std::move(new_vecs);
  }
  return vecs[0];
}

std::vector<TwoBody> sum_equal_elements(std::vector<TwoBody> &vec) {
  // sort(vec);
  int vec_size = vec.size();
  std::vector<TwoBody> new_vec;
  new_vec.reserve(vec_size);

  for (int i = 0; i < vec_size; ++i) {
    double sum = vec[i].me;

    while (i + 1 < vec_size && vec[i].s1 == vec[i + 1].s1 && vec[i].s2 == vec[i + 1].s2 &&
           vec[i].ss1 == vec[i + 1].ss1 && vec[i].ss2 == vec[i + 1].ss2) {
      sum += vec[i + 1].me;
      ++i;
    }
    new_vec.push_back(
        TwoBody(vec[i].s1, vec[i].s2, vec[i].ss1, vec[i].ss2, sum));
  }
  return new_vec;
}

std::vector<TwoBody> reduce_vec(std::vector<std::vector<TwoBody>> &vecs) {

  while (vecs.size() > 1) {
    std::vector<std::vector<TwoBody>> new_vecs((int)vecs.size() / 2);

#pragma omp parallel for
    for (size_t i = 0; i < vecs.size() / 2; ++i) {
      std::cout << omp_get_thread_num() << std::endl;
      std::vector<TwoBody> local_vec = vecs[i * 2];
      local_vec.insert(local_vec.end(), vecs[2 * i + 1].begin(),
                       vecs[2 * i + 1].end());
      local_vec = sum_equal_elements(local_vec);
      new_vecs[i] = local_vec;
    }

    if (vecs.size() % 2 == 1) {
      new_vecs.push_back(std::move(vecs.back()));
    }

    vecs = std::move(new_vecs);
  }
  return vecs[0];
}

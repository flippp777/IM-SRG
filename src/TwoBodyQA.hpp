#include "TwoBody.hpp"
#include <unordered_map>
#include <vector>
#ifndef TWOBODYQA
#define TWOBODYQA

#include "TwoBody.hpp"
class TwoBodyQA {
private:
  // single index maps
  std::unordered_map<int, std::vector<TwoBody *>> index_s1;
  std::unordered_map<int, std::vector<TwoBody *>> index_ss1;

  // two index maps
  std::unordered_map<int, std::unordered_map<int, std::vector<TwoBody *>>>
      index_s1s2;
  std::unordered_map<int, std::unordered_map<int, std::vector<TwoBody *>>>
      index_ss1ss2;
  std::unordered_map<int, std::unordered_map<int, std::vector<TwoBody *>>>
      index_s1ss1;
  std::unordered_map<
      int,
      std::unordered_map<int, std::unordered_map<int, std::vector<TwoBody *>>>>
      index_s1s2ss1;

public:
  std::unordered_map<TwoBodyKey, TwoBody *> map;

  TwoBodyQA();
  TwoBodyQA(const TwoBody element);
  TwoBodyQA(const std::vector<TwoBody> vec);

  void insert(const TwoBody element);

  // check single indices
  std::vector<TwoBody *> check_s1(int s1);
  std::vector<TwoBody *> check_ss1(int ss1);

  // check two indices
  std::vector<TwoBody *> check_s1s2(int s1, int s2);
  std::vector<TwoBody *> check_ss1ss2(int ss1, int ss2);
  std::vector<TwoBody *> check_s1ss1(int s1, int ss1);
  std::vector<TwoBody *> check_s1s2ss1(int s1, int s2, int ss1);
};

#endif // !TWOBODYQA

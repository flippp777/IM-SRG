#include "time.hpp"

#include <iostream>

const void print_duration(
    const std::string name,
    const std::chrono::time_point<std::chrono::high_resolution_clock> start,
    const std::chrono::time_point<std::chrono::high_resolution_clock> stop) {

  std::chrono::duration<double, std::milli> ms_double = stop - start;
  std::cout << name << ms_double.count() << "ms\n";
}

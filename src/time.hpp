#include <chrono>
#include <string>

#ifndef TIME
#define TIME

const void print_duration(
    const std::string name,
    const std::chrono::time_point<std::chrono::high_resolution_clock> start,
    const std::chrono::time_point<std::chrono::high_resolution_clock> stop);

#endif // !TIME

#ifndef RTWEEKEND_UTILITY_H
#define RTWEEKEND_UTILITY_H

#include <limits>
#include <random>

const double pi = 3.1415926535897932385;

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

double random_int(int min, int max) {
    return int(random_double(min, max + 1));
}

#endif //RTWEEKEND_UTILITY_H

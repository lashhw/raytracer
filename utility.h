#ifndef RAYTRACER_UTILITY_H
#define RAYTRACER_UTILITY_H

#include <limits>
#include <random>

const double PI = 3.1415926535897932385;

double degrees_to_radians(double degrees) {
    return degrees * PI / 180.0;
}

double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

#endif //RAYTRACER_UTILITY_H

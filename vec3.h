#ifndef RTWEEKEND_VEC3_H
#define RTWEEKEND_VEC3_H

#include <cmath>
#include <iostream>
#include "utility.h"
using std::sqrt;
using std::round;
using std::fabs;

class vec3 {
private:
    double e[3];
public:
    vec3() { }
    vec3(double x, double y, double z) : e{x, y, z} { }
    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }
    vec3 operator-() const {
        return vec3(-e[0], -e[1], -e[2]);
    }
    double operator[](int x) const {
        return e[x];
    }
    double length_squared() {
        return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
    }
    double length() {
        return sqrt(length_squared());
    }
    vec3 unit_vector() {
        double len = length();
        return vec3(e[0] / len, e[1] / len, e[2] / len);
    }
    bool near_zero() {
        const auto s = 1e-8;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }
    static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }
    static vec3 random(double min, double max) {
        return vec3(random_double(min, max), random_double(min, max), random_double(min, max));
    }
    static vec3 random_in_unit_sphere() {
        while (true) {
            vec3 p = vec3::random(-1, 1);
            if (p.length_squared() <= 1) return p;
        }
    }
    static vec3 random_in_unit_disk() {
        while (true) {
            vec3 p = vec3(random_double(-1,1), random_double(-1,1), 0);
            if (p.length_squared() >= 1) continue;
            return p;
        }
    }
};

vec3 operator+(const vec3 &v1, const vec3 &v2) {
    return vec3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

vec3 operator-(const vec3 &v1, const vec3 &v2) {
    return v1 + (-v2);
}

vec3 operator*(const vec3 &v, double t) {
    return vec3(v.x() * t, v.y() * t, v.z() * t);
}

vec3 operator*(double t, const vec3 &v) {
    return v * t;
}

vec3 operator*(const vec3 &v1, const vec3 &v2) {
    return vec3(v1.x() * v2.x(), v1.y() * v2.y(), v1.z() * v2.z());
}

vec3 operator/(const vec3 &v, double t) {
    return v * (1.0 / t);
}


double dot(const vec3 &u, const vec3 &v) {
    return u.x() * v.x() + u.y() * v.y() + u.z() * v.z();
}

vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.y() * v.z() - u.z() * v.y(),
                u.z() * v.x() - u.x() * v.z(),
                u.x() * v.y() - u.y() * v.x());
}

vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2 * dot(v, n) * n;
}

vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

void write_color(std::ostream &out, const vec3 &color) {
    // color: [0, 1]
    double r = sqrt(color.x());
    double g = sqrt(color.y());
    double b = sqrt(color.z());
    out << round(255.0 * r) << ' '
        << round(255.0 * g) << ' '
        << round(255.0 * b) << '\n';
}

#endif //RTWEEKEND_VEC3_H

#ifndef RTWEEKEND_RAY_H
#define RTWEEKEND_RAY_H

#include "vec3.h"

class ray {
private:
    vec3 orig;
    vec3 dir;
public:
    ray() { }
    ray(const vec3 &orig, const vec3 &dir) : orig(orig), dir(dir) {}
    vec3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    vec3 at(double t) const {
        return orig + t * dir;
    }
};

#endif //RTWEEKEND_RAY_H

#ifndef RAYTRACER_HIT_RECORD_H
#define RAYTRACER_HIT_RECORD_H

#include "ray.h"

using std::shared_ptr;

struct hit_record {
    double t;
    vec3 p;
    bool front_face;
    vec3 normal;
    double u;
    double v;

    void set_normal(const ray &r, const vec3 &outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

#endif //RAYTRACER_HIT_RECORD_H

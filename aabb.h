#ifndef RAYTRACER_AABB_H
#define RAYTRACER_AABB_H

#include "vec3.h"
#include "ray.h"
using std::swap;
using std::min;
using std::max;

class aabb {
public:
    aabb() : minimum(vec3(0, 0, 0)), maximum(vec3(0, 0, 0)) { }
    aabb(const vec3 &minimum, const vec3 &maximum) : minimum(minimum), maximum(maximum) { }
    bool hit(const ray &r, double t_min, double t_max) const {
        for (int i = 0; i < 3; i++) {
            double t1 = (minimum[i] - r.origin()[i]) / r.direction()[i];
            double t2 = (maximum[i] - r.origin()[i]) / r.direction()[i];
            if (t1 > t2) swap(t1, t2);
            t_min = max(t_min, t1);
            t_max = min(t_max, t2);
            if (t_min > t_max) return false;
        }
        return true;
    }
public:
    vec3 minimum, maximum;
};

aabb surrounding_box(aabb box1, aabb box2) {
    vec3 minimum(
         min(box1.minimum.x(), box2.minimum.x()),
         min(box1.minimum.y(), box2.minimum.y()),
         min(box1.minimum.z(), box2.minimum.z())
    );
    vec3 maximum(
         max(box1.maximum.x(), box2.maximum.x()),
         max(box1.maximum.y(), box2.maximum.y()),
         max(box1.maximum.z(), box2.maximum.z())
    );
    return aabb(minimum, maximum);
}

#endif //RAYTRACER_AABB_H

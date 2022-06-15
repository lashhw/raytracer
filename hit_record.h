#ifndef RAYTRACER_HIT_RECORD_H
#define RAYTRACER_HIT_RECORD_H
using std::shared_ptr;

class vec3;
class ray;
class material;
struct hit_record {
    double t;
    vec3 p;
    vec3 normal;
    bool front_face;
    shared_ptr<material> mat_ptr;
    double u;
    double v;

    void set_normal(const ray &r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

#endif //RAYTRACER_HIT_RECORD_H

#ifndef RAYTRACER_HITTABLE_H
#define RAYTRACER_HITTABLE_H

#include <memory>
#include <vector>
#include "ray.h"
#include "hit_record.h"

using std::shared_ptr;
using std::vector;

class hittable {
public:
    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const = 0;
};

class sphere : public hittable {
public:
    sphere(const vec3 &center, double radius, shared_ptr<material> mat_ptr) : center(center), radius(radius),
                                                                              mat_ptr(mat_ptr) {}

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        vec3 amc = r.origin() - center;
        double a = dot(r.direction(), r.direction());
        double b = 2.0 * dot(amc, r.direction());
        double c = dot(amc, amc) - radius * radius;
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) return false;

        double root = (-b - sqrt(discriminant)) / (2 * a);
        if (!(t_min < root && root < t_max)) {
            root = (-b + sqrt(discriminant)) / (2 * a);
            if (!(t_min < root && root < t_max)) {
                return false;
            }
        }

        rec.t = root;
        rec.p = r.at(root);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        get_sphere_uv(outward_normal, rec.u, rec.v);

        return true;
    }

private:
    vec3 center;
    double radius;
    shared_ptr<material> mat_ptr;

    static void get_sphere_uv(const vec3 &p, double &u, double &v) {
        // p: a given p on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        auto theta = acos(-p.y());
        auto phi = atan2(-p.z(), p.x()) + PI;

        u = phi / (2 * PI);
        v = theta / PI;
    }
};

class hittable_list : public hittable {
public:
    vector<shared_ptr<hittable>> objects;
public:
    void add(shared_ptr<hittable> object) { objects.push_back(object); }

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        bool hit_anything = false;
        for (auto &object: objects) {
            if (object->hit(r, t_min, t_max, rec)) {
                hit_anything = true;
                t_max = rec.t;
            }
        }
        return hit_anything;
    }
};

class xy_rect : public hittable {
public:
    xy_rect() {}

    xy_rect(double x0, double x1, double y0, double y1, double z, shared_ptr<material> mat_ptr)
            : x0(x0), x1(x1), y0(y0), y1(y1), z(z), mat_ptr(mat_ptr) {};

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        auto t = (z - r.origin().z()) / r.direction().z();
        if (t < t_min || t > t_max)
            return false;
        auto x = r.origin().x() + t * r.direction().x();
        auto y = r.origin().y() + t * r.direction().y();
        if (x < x0 || x > x1 || y < y0 || y > y1)
            return false;

        rec.p = r.at(t);
        auto outward_normal = vec3(0, 0, 1);
        rec.set_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (y - y0) / (y1 - y0);

        return true;
    }

public:
    double x0, x1, y0, y1, z;
    shared_ptr<material> mat_ptr;
};

#endif //RAYTRACER_HITTABLE_H

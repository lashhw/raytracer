#ifndef RAYTRACER_HITTABLE_H
#define RAYTRACER_HITTABLE_H

#include <memory>
#include <vector>
#include <algorithm>
#include "ray.h"
#include "hit_record.h"
#include "aabb.h"

using std::shared_ptr;
using std::vector;
using std::cerr;
using std::nth_element;
using std::make_shared;

class hittable {
public:
    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const = 0;

    virtual aabb bounding_box() const = 0;
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
        get_sphere_uv(outward_normal, rec.u, rec.v);
        rec.mat_ptr = mat_ptr;

        return true;
    }

    virtual aabb bounding_box() const override {
        return aabb(
                center - vec3(radius, radius, radius),
                center + vec3(radius, radius, radius)
        );
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
    void clear() { objects.clear(); }

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

    virtual aabb bounding_box() const override {
        aabb output_box;
        for (auto x: objects) {
            output_box = surrounding_box(output_box, x->bounding_box());
        }
        return output_box;
    }
};

bool compare_box(const shared_ptr<hittable> &a, const shared_ptr<hittable> &b, size_t axis) {
    aabb box_a, box_b;
    return box_a.minimum[axis] < box_b.minimum[axis];
}

bool compare_box_x(const shared_ptr<hittable> &a, const shared_ptr<hittable> &b) {
    return compare_box(a, b, 0);
}

bool compare_box_y(const shared_ptr<hittable> &a, const shared_ptr<hittable> &b) {
    return compare_box(a, b, 1);
}

bool compare_box_z(const shared_ptr<hittable> &a, const shared_ptr<hittable> &b) {
    return compare_box(a, b, 2);
}

class bvh_node : public hittable {
public:
    bvh_node(vector<shared_ptr<hittable>> &objects, size_t start, size_t end) {
        size_t object_size = end - start;

        if (object_size == 1) {
            left = right = objects[start];
        } else {
            int dice = random_int(0, 2);
            auto comparator = dice == 0 ? compare_box_x : dice == 1 ? compare_box_y : compare_box_z;
            size_t mid = start + object_size / 2;
            sort(objects.begin() + start, objects.begin() + end, comparator);

            if (object_size == 2) {
                left = objects[start];
                right = objects[start + 1];
            } else if (object_size == 3) {
                left = objects[start];
                right = make_shared<bvh_node>(objects, mid, end);
            } else {
                left = make_shared<bvh_node>(objects, start, mid);
                right = make_shared<bvh_node>(objects, mid, end);
            }
        }

        aabb box_left = left->bounding_box();
        aabb box_right = right->bounding_box();

        box = surrounding_box(box_left, box_right);
    }

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        if (!box.hit(r, t_min, t_max)) return false;
        bool hit_left = left->hit(r, t_min, t_max, rec);
        bool hit_right = right->hit(r, t_min, hit_left ? rec.t : t_max, rec);
        return hit_left || hit_right;
    }

    virtual aabb bounding_box() const override {
        return box;
    }

private:
    shared_ptr<hittable> left;
    shared_ptr<hittable> right;
    aabb box;
};

class xy_rect : public hittable {
public:
    xy_rect() {}

    xy_rect(double _x0, double _x1, double _y0, double _y1, double _k,
            shared_ptr<material> mat)
            : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        auto t = (k - r.origin().z()) / r.direction().z();
        if (t < t_min || t > t_max)
            return false;
        auto x = r.origin().x() + t * r.direction().x();
        auto y = r.origin().y() + t * r.direction().y();
        if (x < x0 || x > x1 || y < y0 || y > y1)
            return false;
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (y - y0) / (y1 - y0);
        rec.t = t;
        auto outward_normal = vec3(0, 0, 1);
        rec.set_normal(r, outward_normal);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return true;
    }

    virtual aabb bounding_box() const override {
        // The bounding box must have non-zero width in each dimension, so pad the Z
        // dimension a small amount.
        return aabb(vec3(x0, y0, k - 0.0001), vec3(x1, y1, k + 0.0001));
    }

public:
    shared_ptr<material> mp;
    double x0, x1, y0, y1, k;
};

class xz_rect : public hittable {
public:
    xz_rect() {}

    xz_rect(double _x0, double _x1, double _z0, double _z1, double _k,
            shared_ptr<material> mat)
            : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        auto t = (k - r.origin().y()) / r.direction().y();
        if (t < t_min || t > t_max)
            return false;
        auto x = r.origin().x() + t * r.direction().x();
        auto z = r.origin().z() + t * r.direction().z();
        if (x < x0 || x > x1 || z < z0 || z > z1)
            return false;
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        auto outward_normal = vec3(0, 1, 0);
        rec.set_normal(r, outward_normal);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return true;
    }

    virtual aabb bounding_box() const override {
        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        return aabb(vec3(x0, k - 0.0001, z0), vec3(x1, k + 0.0001, z1));
    }

public:
    shared_ptr<material> mp;
    double x0, x1, z0, z1, k;
};

class yz_rect : public hittable {
public:
    yz_rect() {}

    yz_rect(double _y0, double _y1, double _z0, double _z1, double _k,
            shared_ptr<material> mat)
            : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec) const override {
        auto t = (k - r.origin().x()) / r.direction().x();
        if (t < t_min || t > t_max)
            return false;
        auto y = r.origin().y() + t * r.direction().y();
        auto z = r.origin().z() + t * r.direction().z();
        if (y < y0 || y > y1 || z < z0 || z > z1)
            return false;
        rec.u = (y - y0) / (y1 - y0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        auto outward_normal = vec3(1, 0, 0);
        rec.set_normal(r, outward_normal);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return true;
    }

    virtual aabb bounding_box() const override {
        // The bounding box must have non-zero width in each dimension, so pad the X
        // dimension a small amount.
        return aabb(vec3(k - 0.0001, y0, z0), vec3(k + 0.0001, y1, z1));
    }

public:
    shared_ptr<material> mp;
    double y0, y1, z0, z1, k;
};

#endif //RAYTRACER_HITTABLE_H

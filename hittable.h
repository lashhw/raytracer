#ifndef RAYTRACER_HITTABLE_H
#define RAYTRACER_HITTABLE_H

#include <memory>
#include <vector>
#include "hit_record.h"
#include "material.h"

using std::shared_ptr;
using std::vector;

class hittable {
public:
    virtual bool
    hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const = 0;
};

class sphere : public hittable {
public:
    sphere(const vec3 &center, double radius, shared_ptr<material> mat_ptr) : center(center), radius(radius),
                                                                              mat_ptr(mat_ptr) {}

    virtual bool
    hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
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
        hit_mat = mat_ptr;

        return true;
    }

private:
    vec3 center;
    double radius;
    shared_ptr<material> mat_ptr;

    static void get_sphere_uv(const vec3 &p, double &u, double &v) {
        // p: a given p on the sphere of radius one, centered at the origin.
        // u: returned value [0, 1] of angle from X=-1 to X=+1.
        // v: returned value [0, 1] of angle around the X axis from Y=-1.

        auto phi = atan2(p.z(), -p.y());
        auto theta = acos(-p.x());

        u = theta / PI;
        if (phi < 0) v = (phi + 2 * PI) / (2 * PI);
        else v = phi / (2 * PI);
    }
};

class hittable_list : public hittable {
public:
    vector<shared_ptr<hittable>> objects;
public:
    void add(shared_ptr<hittable> object) { objects.push_back(object); }

    virtual bool
    hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        bool hit_anything = false;
        for (auto &object: objects) {
            if (object->hit(r, t_min, t_max, rec, hit_mat)) {
                hit_anything = true;
                t_max = rec.t;
            }
        }
        return hit_anything;
    }
};

class xy_rect : public hittable {
public:
    xy_rect(double x0, double x1, double y0, double y1, double z, shared_ptr<material> mat_ptr)
            : x0(x0), x1(x1), y0(y0), y1(y1), z(z), mat_ptr(mat_ptr) {};

    virtual bool
    hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        double t = (z - r.origin().z()) / r.direction().z();
        if (t < t_min || t > t_max)
            return false;
        double x = r.origin().x() + t * r.direction().x();
        double y = r.origin().y() + t * r.direction().y();
        if (x < x0 || x > x1 || y < y0 || y > y1)
            return false;

        rec.t = t;
        rec.p = r.at(t);
        vec3 outward_normal = vec3(0, 0, 1);
        rec.set_normal(r, outward_normal);
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (y - y0) / (y1 - y0);
        hit_mat = mat_ptr;

        return true;
    }

private:
    double x0, x1, y0, y1, z;
    shared_ptr<material> mat_ptr;
};

class xz_rect : public hittable {
public:
    xz_rect(double x0, double x1, double z0, double z1, double y, shared_ptr<material> mat_ptr)
            : x0(x0), x1(x1), z0(z0), z1(z1), y(y), mat_ptr(mat_ptr) {};

    virtual bool
    hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        double t = (y - r.origin().y()) / r.direction().y();
        if (t < t_min || t > t_max)
            return false;
        double x = r.origin().x() + t * r.direction().x();
        double z = r.origin().z() + t * r.direction().z();
        if (x < x0 || x > x1 || z < z0 || z > z1)
            return false;

        rec.t = t;
        rec.p = r.at(t);
        vec3 outward_normal = vec3(0, 1, 0);
        rec.set_normal(r, outward_normal);
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (z - z0) / (z1 - z0);
        hit_mat = mat_ptr;

        return true;
    }

private:
    double x0, x1, z0, z1, y;
    shared_ptr<material> mat_ptr;
};

class yz_rect : public hittable {
public:
    yz_rect(double y0, double y1, double z0, double z1, double x, shared_ptr<material> mat_ptr)
            : y0(y0), y1(y1), z0(z0), z1(z1), x(x), mat_ptr(mat_ptr) {};

    virtual bool
    hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        double t = (x - r.origin().x()) / r.direction().x();
        if (t < t_min || t > t_max)
            return false;
        double y = r.origin().y() + t * r.direction().y();
        double z = r.origin().z() + t * r.direction().z();
        if (y < y0 || y > y1 || z < z0 || z > z1)
            return false;

        rec.t = t;
        rec.p = r.at(t);
        vec3 outward_normal = vec3(1, 0, 0);
        rec.set_normal(r, outward_normal);
        rec.u = (y - y0) / (y1 - y0);
        rec.v = (z - z0) / (z1 - z0);
        hit_mat = mat_ptr;

        return true;
    }

private:
    double y0, y1, z0, z1, x;
    shared_ptr<material> mat_ptr;
};

class box : public hittable {
public:
    box(const vec3 &p0, const vec3 &p1, shared_ptr<material> mat_ptr) {
        sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), mat_ptr));
        sides.add(make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), mat_ptr));

        sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), mat_ptr));
        sides.add(make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), mat_ptr));

        sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), mat_ptr));
        sides.add(make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), mat_ptr));
    }

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        return sides.hit(r, t_min, t_max, rec, hit_mat);
    }

private:
    hittable_list sides;
};

class translate : public hittable {
public:
    translate(shared_ptr<hittable> hittable_ptr, const vec3& offset)
            : hittable_ptr(hittable_ptr), offset(offset) {}

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        ray moved_r(r.origin() - offset, r.direction());
        if (!hittable_ptr->hit(moved_r, t_min, t_max, rec, hit_mat))
            return false;

        rec.p = rec.p + offset;
        return true;
    }

public:
    shared_ptr<hittable> hittable_ptr;
    vec3 offset;
};

class rotate_x : public hittable {
public:
    rotate_x(shared_ptr<hittable> hittable_ptr, double angle) : hittable_ptr(hittable_ptr) {
        double radians = degrees_to_radians(angle);
        sin_theta = sin(radians);
        cos_theta = cos(radians);
    }

    virtual bool hit(const ray &r, double t_min, double t_max, hit_record &rec, shared_ptr<material> &hit_mat) const override {
        vec3 origin = r.origin();
        vec3 direction = r.direction();

        origin.e[1] = cos_theta * r.origin()[1] - sin_theta * r.origin()[2];
        origin.e[2] = sin_theta * r.origin()[1] + cos_theta * r.origin()[2];

        direction.e[1] = cos_theta * r.direction()[1] - sin_theta * r.direction()[2];
        direction.e[2] = sin_theta * r.direction()[1] + cos_theta * r.direction()[2];

        ray rotated_r(origin, direction);

        if (!hittable_ptr->hit(rotated_r, t_min, t_max, rec, hit_mat))
            return false;

        vec3 p = rec.p;
        vec3 normal = rec.normal;

        p.e[1] =  cos_theta * rec.p[1] + sin_theta * rec.p[2];
        p.e[2] = -sin_theta * rec.p[1] + cos_theta * rec.p[2];

        normal.e[1] =  cos_theta * rec.normal[1] + sin_theta * rec.normal[2];
        normal.e[2] = -sin_theta * rec.normal[1] + cos_theta * rec.normal[2];

        rec.p = p;
        rec.normal = normal;

        return true;
    }

public:
    shared_ptr<hittable> hittable_ptr;
    double sin_theta;
    double cos_theta;
};

#endif //RAYTRACER_HITTABLE_H

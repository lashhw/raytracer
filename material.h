#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H

#include "texture.h"
#include "hittable.h"
#include "hit_record.h"
using std::shared_ptr;
using std::make_shared;

class material {
public:
    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const = 0;
    virtual vec3 emitted(double u, double v, const vec3& p) const {
        return vec3(0,0,0);
    }
};

class lambertian : public material {
public:
    lambertian(const vec3& a) : albedo(make_shared<solid_color>(a)) { }
    lambertian(shared_ptr<texture> a) : albedo(a) { }

    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const override {
        vec3 scatter_direction = rec.normal + vec3::random_in_unit_sphere().unit_vector();
        if (scatter_direction.near_zero()) scatter_direction = rec.normal;
        scattered = ray(rec.p, scatter_direction);
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

    shared_ptr<texture> albedo;
};

class metal : public material {
private:
    vec3 albedo;
    double fuzz;

public:
    metal(const vec3& albedo, double fuzz) : albedo(albedo), fuzz(fuzz) {}

    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const override {
        vec3 reflected = reflect(r_in.direction().unit_vector(), rec.normal);
        attenuation = albedo;
        scattered = ray(rec.p, reflected + fuzz * vec3::random_in_unit_sphere());
        return (dot(scattered.direction(), rec.normal) > 0);
    }
};

class dielectric : public material {
private:
    double ir; // Index of Refraction
    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1-ref_idx) / (1+ref_idx);
        r0 = r0*r0;
        return r0 + (1-r0)*pow((1 - cosine),5);
    }

public:
    dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const override {
        attenuation = vec3(1.0, 1.0, 1.0);
        double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

        vec3 unit_direction = r_in.direction().unit_vector();
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction);
        return true;
    }
};

class diffuse_light : public material  {
public:
    diffuse_light(shared_ptr<texture> a) : emit(a) {}
    diffuse_light(vec3 c) : emit(make_shared<solid_color>(c)) {}

    virtual bool scatter(
            const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered
    ) const override {
        return false;
    }

    virtual vec3 emitted(double u, double v, const vec3& p) const override {
        return emit->value(u, v, p);
    }

public:
    shared_ptr<texture> emit;
};

#endif //RAYTRACER_MATERIAL_H

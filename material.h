#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H

#include "texture.h"
#include "hit_record.h"

using std::shared_ptr;
using std::make_shared;

class material {
public:
    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const = 0;

    virtual vec3 emitted(double u, double v) const {
        return vec3(0, 0, 0);
    }
};

class lambertian : public material {
public:
    lambertian(const vec3 &albedo) : texture_ptr(make_shared<solid_color>(albedo)) {}

    lambertian(shared_ptr<texture> texture_ptr) : texture_ptr(texture_ptr) {}

    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const override {
        vec3 scatter_direction = rec.normal + vec3::random_in_unit_sphere().unit_vector();
        if (scatter_direction.near_zero()) scatter_direction = rec.normal;
        attenuation = texture_ptr->value(rec.u, rec.v);
        scattered = ray(rec.p, scatter_direction);
        return true;
    }

private:
    shared_ptr<texture> texture_ptr;
};

class metal : public material {
public:
    metal(const vec3 &albedo, double fuzz) : texture_ptr(make_shared<solid_color>(albedo)), fuzz(fuzz) {}

    metal(shared_ptr<texture> texture_ptr, double fuzz) : texture_ptr(texture_ptr), fuzz(fuzz) {}

    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const override {
        vec3 reflected = reflect(r_in.direction().unit_vector(), rec.normal);
        attenuation = texture_ptr->value(rec.u, rec.v);
        scattered = ray(rec.p, reflected + fuzz * vec3::random_in_unit_sphere());
        return true;
    }

private:
    shared_ptr<texture> texture_ptr;
    double fuzz;
};

class dielectric : public material {
public:
    dielectric(double index_of_refraction) : index_of_refraction(index_of_refraction) {}

    virtual bool
    scatter(const ray &r_in, const hit_record &rec, vec3 &attenuation,
            ray &scattered) const override {
        attenuation = vec3(1.0, 1.0, 1.0);
        double refraction_ratio = rec.front_face ? (1.0 / index_of_refraction) : index_of_refraction;

        vec3 unit_direction = r_in.direction().unit_vector();
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;

        vec3 direction;
        if (cannot_refract || random_double() < reflectance(cos_theta))
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction);
        return true;
    }

private:
    double index_of_refraction;
    double reflectance(double cos_theta) const {
        // Use Schlick's approximation for reflectance
        auto r0 = (1 - index_of_refraction) / (1 + index_of_refraction);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cos_theta), 5);
    }
};

class diffuse_light : public material {
public:
    diffuse_light(shared_ptr<texture> emit) : emit(emit) {}

    diffuse_light(vec3 c) : emit(make_shared<solid_color>(c)) {}

    virtual bool scatter(
            const ray &r_in, const hit_record &rec, vec3 &attenuation, ray &scattered
    ) const override {
        return false;
    }

    virtual vec3 emitted(double u, double v) const override {
        return emit->value(u, v);
    }

public:
    shared_ptr<texture> emit;
};

#endif //RAYTRACER_MATERIAL_H

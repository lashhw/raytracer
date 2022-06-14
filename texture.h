#ifndef RTWEEKEND_TEXTURE_H
#define RTWEEKEND_TEXTURE_H

#include "vec3.h"

class texture {
public:
    virtual vec3 value(double u, double v, const vec3& p) const = 0;
};

class solid_color : public texture {
public:
    solid_color() {}
    solid_color(vec3 c) : color_value(c) {}

    solid_color(double red, double green, double blue)
            : solid_color(vec3(red,green,blue)) {}

    virtual vec3 value(double u, double v, const vec3& p) const override {
        return color_value;
    }

private:
    vec3 color_value;
};

class checker_texture : public texture {
public:
    checker_texture() {}

    checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
            : even(_even), odd(_odd) {}

    checker_texture(vec3 c1, vec3 c2)
            : even(make_shared<solid_color>(c1)) , odd(make_shared<solid_color>(c2)) {}

    virtual vec3 value(double u, double v, const vec3& p) const override {
        auto sines = sin(10*p.x())*sin(10*p.y())*sin(10*p.z());
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }

public:
    shared_ptr<texture> odd;
    shared_ptr<texture> even;
};

#endif //RTWEEKEND_TEXTURE_H

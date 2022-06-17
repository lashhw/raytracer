#ifndef RAYTRACER_TEXTURE_H
#define RAYTRACER_TEXTURE_H

#include <memory>

#include "vec3.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using std::make_shared;
using std::shared_ptr;
using std::cerr;

class texture {
public:
    virtual vec3 value(double u, double v) const = 0;
};

class solid_color : public texture {
public:
    solid_color(vec3 c) : color_value(c) {}

    solid_color(double red, double green, double blue)
            : solid_color(vec3(red, green, blue)) {}

    virtual vec3 value(double u, double v) const override {
        return color_value;
    }

private:
    vec3 color_value;
};

class checker_texture : public texture {
public:
    checker_texture(shared_ptr<texture> even, shared_ptr<texture> odd, int squares_per_axis)
            : even(even), odd(odd), squares_per_axis(squares_per_axis) {}

    checker_texture(vec3 c1, vec3 c2, int squares_per_axis)
            : checker_texture(make_shared<solid_color>(c1), make_shared<solid_color>(c2), squares_per_axis) {}

    virtual vec3 value(double u, double v) const override {
        int u_pos = u * squares_per_axis;
        int v_pos = v * squares_per_axis;
        if ((u_pos + v_pos) % 2 == 1) return odd->value(u, v);
        else return even->value(u, v);
    }

private:
    shared_ptr<texture> odd;
    shared_ptr<texture> even;
    int squares_per_axis;
};

class image_texture : public texture {
public:
    image_texture(const char* filename) {
        int n;
        data = stbi_load(filename, &width, &height, &n, 3);

        if (!data) {
            cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
        }

        bytes_per_scanline = 3 * width;
    }

    ~image_texture() {
        stbi_image_free(data);
    }

    virtual vec3 value(double u, double v) const override {
        // If we have no texture data, then return solid cyan as a debugging aid.
        if (data == nullptr)
            return vec3(0, 1, 1);

        auto i = int(u * height);
        auto j = int(v * width);

        // Clamp integer mapping
        clamp(i, 0, height - 1);
        clamp(j, 0, width - 1);

        unsigned char *pixel = data + i * bytes_per_scanline + j * 3;

        vec3 color(pixel[0], pixel[1], pixel[2]);
        return color / 255.0;
    }

private:
    unsigned char *data;
    int width, height;
    int bytes_per_scanline;
};

#endif //RAYTRACER_TEXTURE_H

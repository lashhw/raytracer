#include "raytracer.h"

#include <iostream>
#include <memory>
#include <fstream>
using std::make_shared;

vec3 ray_color(const ray& r, const vec3& background, const hittable& world, int depth) {
    hit_record rec;
    if (depth <= 0) return vec3(0, 0, 0);
    if (world.hit(r, 0.001, INFINITY, rec)) {
        vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
        vec3 attenuation;
        ray scattered;
        if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return emitted;
        return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
    }
    return background;
}

void add_random_sphere(hittable_list &world) {
    double x = 0;
    double y = random_double(-7, 7);
    double z = random_double(-7, 7);
    int dice = rand() % 10;
    shared_ptr<material> material;
    if (dice <= 1) material = make_shared<dielectric>(1.5);
    else if (dice <= 5) material = make_shared<metal>(vec3(random_double(), random_double(), random_double()), random_double(0, 0.5));
    else material = make_shared<lambertian>(vec3(random_double(), random_double(), random_double()));
    world.add(make_shared<sphere>(vec3(x, y, z), 0.5, material));
}

int main() {
    std::ofstream file("image.ppm");

    const double ASPECT_RATIO = 16.0 / 9.0;
    const int IMAGE_WIDTH = 500;
    const int IMAGE_HEIGHT = IMAGE_WIDTH / ASPECT_RATIO;
    const int SAMPLES_PER_PIXEL = 50;
    const int MAX_DEPTH = 50;

    hittable_list world;
    auto checker = make_shared<checker_texture>(vec3(0.2, 0.3, 0.1), vec3(0.9, 0.9, 0.9));
    shared_ptr<material> material_ground = make_shared<lambertian>(checker);
    shared_ptr<material> material_big_sphere = make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(vec3(1000.5,  0.0, 0.0), 1000, material_ground));
    world.add(make_shared<sphere>(vec3(-1.5, 3, 3), 2, material_big_sphere));
    for (int i = 1; i <= 30; i++) add_random_sphere(world);
    shared_ptr<material> difflight = make_shared<diffuse_light>(vec3(4, 4, 4));
    world.add(make_shared<xy_rect>(-1, 1, -5, -2, -5, difflight));

    vec3 lookfrom(-6, 0, 25), lookat(0, 0, 0);
    camera cam(
        lookfrom, lookat, vec3(-1, 0, 0),
        20, ASPECT_RATIO, 0.1, (lookat - lookfrom).length()
    );

    vec3 background(0.1, 0.1, 0.1);

    file << "P3\n" << IMAGE_WIDTH << ' ' << IMAGE_HEIGHT << "\n255\n";
    for (int i = 0; i < IMAGE_HEIGHT; i++) {
        std::cout << "Processing row " << i << '\n';
        for (int j = 0; j < IMAGE_WIDTH; j++) {
            vec3 pixel_color = vec3(0, 0, 0);
            for (int k = 0; k < SAMPLES_PER_PIXEL; k++) {
                double u = (i + random_double(-1, 1)) / (IMAGE_HEIGHT - 1);
                double v = (j + random_double(-1, 1)) / (IMAGE_WIDTH - 1);
                pixel_color = pixel_color + ray_color(cam.get_ray(u, v), background, world, MAX_DEPTH);
            }
            pixel_color = pixel_color / SAMPLES_PER_PIXEL;
            pixel_color.write_color(file);
        }
    }
}

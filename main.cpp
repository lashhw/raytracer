#include "raytracer.h"

#include <iostream>
#include <memory>
#include <fstream>
using std::make_shared;

vec3 ray_color(const ray& r, const vec3& background, const hittable& world, int depth) {
    hit_record rec;
    shared_ptr<material> hit_mat;
    if (depth <= 0) return vec3(0, 0, 0);
    if (world.hit(r, 0.001, INFINITY, rec, hit_mat)) {
        vec3 emitted = hit_mat->emitted(rec.u, rec.v, rec.p);
        vec3 attenuation;
        ray scattered;
        if (!hit_mat->scatter(r, rec, attenuation, scattered))
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

    /*
    const double ASPECT_RATIO = 16.0 / 9.0;
    const int image_width = 500;
    const int image_height = image_width / ASPECT_RATIO;
    const int samples_per_pixel = 2;
    const int max_depth = 8;

    hittable_list world;
    shared_ptr<material> ground_material = make_shared<lambertian>(vec3(0.9, 0.9, 0.9));
    shared_ptr<material> big_sphere_material = make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0);
    shared_ptr<texture> earth_texture = make_shared<image_texture>("../earthmap.jpg");
    shared_ptr<material> earth_material = make_shared<lambertian>(earth_texture);
    world.add(make_shared<sphere>(vec3(1000.5,  0.0, 0.0), 1000, ground_material));
    world.add(make_shared<sphere>(vec3(-1.5, 3, 3), 2, big_sphere_material));
    world.add(make_shared<sphere>(vec3(-0.5, 0, -10), 1, earth_material));
    for (int i = 1; i <= 30; i++) add_random_sphere(world);
    shared_ptr<material> difflight = make_shared<diffuse_light>(vec3(10, 10, 10));
    world.add(make_shared<xy_rect>(-1, 1, -5, -2, -5, difflight));

    vec3 lookfrom(-6, 0, 25), lookat(0, 0, 0);
    camera cam(
        lookfrom, lookat, vec3(-1, 0, 0),
        20, ASPECT_RATIO, 0.1, (lookat - lookfrom).length()
    );

    vec3 background(0.1, 0.1, 0.1);
     */

    double aspect_ratio = 1.0;
    int image_width = 600;
    int image_height = image_width / aspect_ratio;
    int samples_per_pixel = 200;
    int max_depth = 8;
    vec3 background = vec3(0, 0, 0);
    vec3 lookfrom(278, 278, -800), lookat(278, 278, 0);
    double vfov = 40.0;
    double aperture = 0.0;

    shared_ptr<material> red   = make_shared<lambertian>(vec3(.65, .05, .05));
    shared_ptr<material> white = make_shared<lambertian>(vec3(.73, .73, .73));
    shared_ptr<material> green = make_shared<lambertian>(vec3(.12, .45, .15));
    shared_ptr<material> light = make_shared<diffuse_light>(vec3(15, 15, 15));

    hittable_list world;
    world.add(make_shared<xz_rect>(0, 555, 0, 555, 555, green));
    world.add(make_shared<xz_rect>(0, 555, 0, 555, 0, red));
    world.add(make_shared<yz_rect>(213, 343, 227, 332, 1, light));
    world.add(make_shared<yz_rect>(0, 555, 0, 555, 0, white));
    world.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
    world.add(make_shared<yz_rect>(0, 555, 0, 555, 555, white));
    shared_ptr<hittable> box1 = make_shared<box>(vec3(0, 0, 0), vec3(330, 165, 165), white);
    box1 = make_shared<rotate_x>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(225, 265, 295));
    world.add(box1);
    shared_ptr<hittable> box2 = make_shared<box>(vec3(0, 0, 0), vec3(165, 165, 165), white);
    box2 = make_shared<rotate_x>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(390, 130, 65));
    world.add(box2);

    camera cam(
            lookfrom, lookat, vec3(-1, 0, 0),
            vfov, aspect_ratio, aperture, (lookat - lookfrom).length()
    );

    file << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    for (int i = 0; i < image_height; i++) {
        std::cout << "Processing row " << i << '\n';
        for (int j = 0; j < image_width; j++) {
            vec3 pixel_color = vec3(0, 0, 0);
            for (int k = 0; k < samples_per_pixel; k++) {
                double u = (i + random_double(-1, 1)) / (image_height - 1);
                double v = (j + random_double(-1, 1)) / (image_width - 1);
                pixel_color = pixel_color + ray_color(cam.get_ray(u, v), background, world, max_depth);
            }
            pixel_color = pixel_color / samples_per_pixel;
            pixel_color.write_color(file);
        }
    }
}

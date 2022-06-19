#include <iostream>
#include <memory>
#include <fstream>

#include "camera.h"
#include "hittable.h"

using std::make_shared;

vec3 ray_color(const ray& r, const vec3& background, const hittable& world, int depth) {
    hit_record rec;
    shared_ptr<material> hit_mat;
    if (depth <= 0) return vec3(0, 0, 0);
    if (world.hit(r, 0.001, INFINITY, rec, hit_mat)) {
        vec3 emitted = hit_mat->emitted(rec.u, rec.v);
        vec3 attenuation;
        ray scattered;
        if (!hit_mat->scatter(r, rec, attenuation, scattered))
            return emitted;
        return emitted + attenuation * ray_color(scattered, background, world, depth - 1);
    }
    return background;
}

int main() {
    double aspect_ratio = 1.0;
    int image_width = 600;
    int image_height = image_width / aspect_ratio;
    int samples_per_pixel = 10000;
    int max_depth = 16;
    vec3 background = vec3(0, 0, 0);
    vec3 lookfrom(278, 278, -800), lookat(278, 278, 0);
    double vfov = 40.0;
    double aperture = 0.0;

    shared_ptr<material> red   = make_shared<lambertian>(vec3(.65, .05, .05));
    shared_ptr<material> white = make_shared<lambertian>(vec3(.73, .73, .73));
    shared_ptr<material> green = make_shared<lambertian>(vec3(.12, .45, .15));
    shared_ptr<material> light = make_shared<diffuse_light>(vec3(15, 15, 15));

    hittable_list world;
    world.add(make_shared<xz_rect>(0, 555, 0, 555, 555, red));
    world.add(make_shared<xz_rect>(0, 555, 0, 555, 0, green));
    world.add(make_shared<yz_rect>(213, 343, 227, 332, 1, light));
    world.add(make_shared<yz_rect>(0, 555, 0, 555, 0, white));
    world.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));
    world.add(make_shared<yz_rect>(0, 555, 0, 555, 555, white));
    /*
    shared_ptr<hittable> box1 = make_shared<box>(vec3(0, 0, 0), vec3(330, 165, 165), white);
    box1 = make_shared<rotate_x>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(225, 265, 295));
    world.add(box1);
    shared_ptr<hittable> box2 = make_shared<box>(vec3(0, 0, 0), vec3(165, 165, 165), white);
    box2 = make_shared<rotate_x>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(390, 130, 65));
    world.add(box2);
     */
    shared_ptr<material> mirror_material = make_shared<metal>(vec3(0.8, 0.8, 0.9), 0.0);
    shared_ptr<material> ball_material = make_shared<metal>(vec3(0.9, 0.73, 0.05), 1.0);
    shared_ptr<material> glass_material = make_shared<dielectric>(1.5);
    shared_ptr<material> earth_material = make_shared<lambertian>(make_shared<image_texture>("../earthmap.jpg"));
    world.add(make_shared<sphere>(vec3(475, 175, 300), 80, mirror_material));
    world.add(make_shared<sphere>(vec3(475, 380, 150), 80, ball_material));
    world.add(make_shared<sphere>(vec3(515, 100, 80), 40, earth_material));
    world.add(make_shared<sphere>(vec3(200, 350, 100), 80, glass_material));

    camera cam(
            lookfrom, lookat, vec3(-1, 0, 0),
            vfov, aspect_ratio, aperture, (lookat - lookfrom).length()
    );

    std::ofstream file("image.ppm");
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

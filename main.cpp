#include <iostream>
#include <memory>
#include <fstream>
#include <bvh/bvh.hpp>
#include <bvh/triangle.hpp>
#include <bvh/sweep_sah_builder.hpp>
#include <bvh/primitive_intersectors.hpp>
#include <bvh/single_ray_traverser.hpp>

#include "camera.h"
#include "hittable.h"
#include "happly.h"

using std::make_shared;
using Vector3 = bvh::Vector3<double>;
using Triangle = bvh::Triangle<double>;
using Ray = bvh::Ray<double>;
using Bvh = bvh::Bvh<double>;

vec3 ray_color(const ray &r,
               const vec3 &background,
               const std::vector<Triangle> &world,
               const bvh::ClosestPrimitiveIntersector<Bvh, Triangle> &intersector,
               const bvh::SingleRayTraverser<Bvh> &traverser,
               int depth) {
    /*
    hit_record rec;
    shared_ptr<material> hit_mat;
     */
    if (depth <= 0) return vec3(0, 0, 0);
    Ray R(Vector3(r.origin()[0], r.origin()[1], r.origin()[2]),
            Vector3(r.direction()[0], r.direction()[1], r.direction()[2]),
            0.001,
            INFINITY);
    if (auto hit = traverser.traverse(R, intersector)) {
        Vector3 e1 = world[hit->primitive_index].e1;
        Vector3 e2 = world[hit->primitive_index].e2;
        vec3 e1_vec3(e1[0], e1[1], e1[2]);
        vec3 e2_vec3(e2[0], e2[1], e2[2]);

        hit_record rec;
        rec.t = hit->intersection.t;
        rec.p = r.at(rec.t);
        rec.front_face = true;
        rec.normal = cross(e1_vec3, e2_vec3).unit_vector();
        rec.u = hit->intersection.u;
        rec.v = hit->intersection.v;
        vec3 attenuation;
        ray scattered;
        lambertian(vec3(0.9, 0.73, 0.05)).scatter(r, rec, attenuation, scattered);
        return attenuation * ray_color(scattered, background, world,
                                       intersector, traverser, depth - 1);
    }
    return background;
}

int main() {
    double aspect_ratio = 1.0;
    int image_width = 600;
    int image_height = image_width / aspect_ratio;
    int samples_per_pixel = 10;
    int max_depth = 16;
    vec3 background = vec3(1, 1, 1);
    vec3 lookfrom(0, 0.1, 0.5), lookat(0, 0.1, 0);
    double vfov = 40.0;
    vec3 vup = vec3(0, 1, 0);
    double aperture = 0.0;

    happly::PLYData plyIn("../bunny.ply");
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();

    std::vector<Triangle> world;
    for (const auto &idx : fInd) {
        Vector3 points[3];
        for (int i = 0; i < 3; i++) {
            points[i] = Vector3(vPos[idx[i]][0], vPos[idx[i]][1], vPos[idx[i]][2]);
        }
        world.emplace_back(points[0], points[1], points[2]);
    }

    Bvh bvh;
    auto [bboxes, centers] = bvh::compute_bounding_boxes_and_centers(world.data(), world.size());
    auto global_bbox = bvh::compute_bounding_boxes_union(bboxes.get(), world.size());

    bvh::SweepSahBuilder<Bvh> builder(bvh);
    builder.build(global_bbox, bboxes.get(), centers.get(), world.size());

    bvh::ClosestPrimitiveIntersector<Bvh, Triangle> intersector(bvh, world.data());
    bvh::SingleRayTraverser<Bvh> traverser(bvh);

    /* cornell box config
    double aspect_ratio = 1.0;
    int image_width = 600;
    int image_height = image_width / aspect_ratio;
    int samples_per_pixel = 1;
    int max_depth = 16;
    vec3 background = vec3(1, 1, 1);
    vec3 lookfrom(278, 278, -800), lookat(278, 278, 0);
    double vfov = 40.0;
    vec3 vup = vec3(-1, 0, 0);
    double aperture = 0.0;
     */

    /* walls
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
     */

    /* two boxes
    shared_ptr<hittable> box1 = make_shared<box>(vec3(0, 0, 0), vec3(330, 165, 165), white);
    box1 = make_shared<rotate_x>(box1, 15);
    box1 = make_shared<translate>(box1, vec3(225, 265, 295));
    world.add(box1);
    shared_ptr<hittable> box2 = make_shared<box>(vec3(0, 0, 0), vec3(165, 165, 165), white);
    box2 = make_shared<rotate_x>(box2, -18);
    box2 = make_shared<translate>(box2, vec3(390, 130, 65));
    world.add(box2);
     */

    /* balls
    shared_ptr<material> mirror_material = make_shared<metal>(vec3(0.8, 0.8, 0.9), 0.0);
    shared_ptr<material> ball_material = make_shared<metal>(vec3(0.9, 0.73, 0.05), 1.0);
    shared_ptr<material> glass_material = make_shared<dielectric>(1.5);
    shared_ptr<material> earth_material = make_shared<lambertian>(make_shared<image_texture>("../earthmap.jpg"));
    world.add(make_shared<sphere>(vec3(475, 175, 300), 80, mirror_material));
    world.add(make_shared<sphere>(vec3(475, 380, 150), 80, ball_material));
    world.add(make_shared<sphere>(vec3(515, 100, 80), 40, earth_material));
    world.add(make_shared<sphere>(vec3(200, 350, 100), 80, glass_material));
     */

    camera cam(
            lookfrom, lookat, vup,
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
                pixel_color = pixel_color + ray_color(cam.get_ray(u, v), background, world,
                                                      intersector, traverser, max_depth);
            }
            pixel_color = pixel_color / samples_per_pixel;
            pixel_color.write_color(file);
        }
    }
}

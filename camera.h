#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H

#include "ray.h"

class camera {
public:
    camera(
            vec3 lookfrom,
            vec3 lookat,
            vec3 vup,
            double vfov_deg,
            double aspect_ratio,
            double aperture,
            double focus_dist
    ) {
        double vfov_rad = degrees_to_radians(vfov_deg);
        double viewpoint_height = 2.0 * tan(vfov_rad / 2);
        double viewpoint_width = viewpoint_height * aspect_ratio;

        w = (lookfrom - lookat).unit_vector();
        u = -(vup - dot(vup, w) * w).unit_vector();
        v = cross(w, u);

        origin = lookfrom;
        vertical = focus_dist * viewpoint_height * u;
        horizontal = focus_dist * viewpoint_width * v;
        upper_left_corner = origin - vertical / 2 - horizontal / 2 - focus_dist * w;
        lens_radius = aperture / 2;
    }

    ray get_ray(double s, double t) {
        vec3 rd = lens_radius * vec3::random_in_unit_disk();
        vec3 offset = u * rd.x() + v * rd.y();
        return ray(
                origin + offset,
                upper_left_corner + s * vertical + t * horizontal - (origin + offset)
        );
    }

private:
    vec3 origin;
    vec3 vertical;
    vec3 horizontal;
    vec3 upper_left_corner;
    vec3 u;
    vec3 v;
    vec3 w;
    double lens_radius;
};

#endif //RAYTRACER_CAMERA_H

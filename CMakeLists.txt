cmake_minimum_required(VERSION 3.21)
project(raytracer)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE O3)

add_executable(raytracer main.cpp aabb.h camera.h hit_record.h hittable.h material.h ray.h rtweekend.h texture.h utility.h vec3.h)

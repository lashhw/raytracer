project(raytracer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE O3)

add_executable(raytracer main.cpp camera.h hit_record.h hittable.h material.h ray.h stb_image.h texture.h utility.h vec3.h)

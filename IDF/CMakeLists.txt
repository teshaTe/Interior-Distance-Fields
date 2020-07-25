project(Interior_df_lib)
cmake_minimum_required(VERSION 3.0)

set(DMAPS_SRC external/dmaps/src/diffusion_map.cpp
              external/dmaps/src/distance_matrix.cpp
              external/dmaps/src/metrics.cpp)

set(DMAPS_HRC external/dmaps/include/diffusion_map.h
              external/dmaps/include/distance_matrix.h
              external/dmaps/include/metrics.h
              external/dmaps/include/types.h)

set(HEADERS   include/interiorDistanceField.h
              include/barycentricCoords.h
              include/frep2D.h
              include/frep3D.h
              include/render.h
              include/timer.hpp )

set(SOURCES   src/interiorDistanceField.cpp
              src/barycentricCoords.cpp
              src/frep2D.cpp
              src/frep3D.cpp
              src/render.cpp )

add_library(IDF ${DMAPS_SRC} ${DMAPS_HRC} ${HEADERS} ${SOURCES})
target_include_directories(IDF PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(IDF Eigen3::Eigen ${OpenGL_LIBARARIES} ${GTS_LIBRARIES} ${GLIB_LDFLAGS} gomp omp omptarget
                          Python::Python Python::NumPy igl_tetgen glad glfw igl igl_opengl igl_opengl_glfw triangle tetgen
                          sfml-window sfml-system sfml-graphics)

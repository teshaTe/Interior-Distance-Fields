project(Interior_dist_fields)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(PROJECT_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include_directories( BEFORE ./ ./include )

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS " -lstdc++ -lm ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3")

find_package(Eigen3 3.3 REQUIRED)
find_package(Python COMPONENTS Development NumPy)
find_package(OpenGL)
find_package(OpenMP)
find_package(GTS)
find_package(PkgConfig REQUIRED)
pkg_search_module(GLIB REQUIRED glib-2.0)
include_directories(${GLIB_INCLUDE_DIRS})

add_definitions(-DUSE_MATPLOTLIB)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    message("clang compiler is used!")
    set(CMAKE_CXX_FLAGS "-fopenmp=libomp ${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message("intel compiler is used!")
    set(CMAKE_CXX_FLAGS "-qopenmp ${CMAKE_CXX_FLAGS}")
else()
    message("gcc compiler is used!")
    set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -qopt-report=5 -xCORE-AVX2")
endif()

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
              include/timer.hpp )

set(SOURCES   src/interiorDistanceField.cpp
              src/barycentricCoords.cpp
              src/frep2D.cpp
              src/frep3D.cpp)

link_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/libIgl/lib")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external/libIgl/include")

add_executable(diffMapTest samples/diffMapsTest.cpp ${DMAPS_SRC} ${DMAPS_HRC})
target_include_directories(diffMapTest PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(diffMapTest Eigen3::Eigen gomp omp omptarget Python::Python Python::NumPy )

add_executable(idfClassTest samples/idfClassTest.cpp ${DMAPS_SRC} ${DMAPS_HRC} ${SOURCES} ${HEADERS})
target_include_directories(idfClassTest PRIVATE ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS})
target_link_libraries(idfClassTest Eigen3::Eigen gomp omp omptarget igl_tetgen ${GTS_LIBRARIES}
                                   glad glfw igl igl_opengl igl_opengl_glfw triangle tetgen
                                   ${OpenGL_LIBARARIES} Python::Python Python::NumPy ${GLIB_LDFLAGS})

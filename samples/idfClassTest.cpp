#include "include/interiorDistanceField.h"
#include "include/frep2D.h"
#include "include/frep3D.h"

#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>

#include <iostream>
#include <functional>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <gtsconfig.h>
#include <gts.h>

static void heart(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj3D frep(g.nx, g.ny, g.nz, 50.0);

   for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
        f[i][j] = -frep.heart3D(glm::f32vec3(x, y, z), glm::f32vec3(0,0,0));
}

static void sphere (gdouble ** f, GtsCartesianGrid g, guint k, gpointer data)
{
  gdouble x, y, z = g.z;
  guint i, j;

  for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
    for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
      f[i][j] = x*x + y*y + z*z;
}

int main(int argc, char *argv[])
{
    // Create the boundary of a square
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    int sizeZeroLS = 46;
    V.resize(sizeZeroLS, 2);
    E.resize(sizeZeroLS, 2);

    V << -2,-2,  -1.5,-2,  -1,-2,  -0.5,-2,  -0.5,-1.5,  -0.5,-1,  -0.25,-1,  0,-1,
          0.25,-1,  0.5,-1,  0.5,-1.5,  0.5,-2,  1,-2,  1.5,-2,  2,-2,
          2,-1.5,  2,-1,  2,-0.5,  2,0,  2,0.5,  2,1,  2,1.5,  2,2,
          1.5,2,  1,2,  0.5,2,  0.5,1.5,  0.5,1,  0.5,0.25,  0.25,0.25,  0,0.25,
          -0.25,0.25,  -0.5,0.25,  -0.5,1,  -0.5,1.5,  -0.5,2,
          -1.0,2,  -1.5,2,  -2,2,  -2,1.5,  -2,1,  -2,0.5, -2,0,  -2,-0.5,  -2,-1,  -2,-1.5;
    E << 0,1, 1,2, 2,3, 3,4,   4,5,   5,6,
         6,7, 7,8, 8,9, 9,10, 10,11, 11,12,
         12,13, 13,14, 14,15, 15,16, 16,17,
         17,18, 18,19, 19,20, 20,21, 21,22,
         22,23, 23,24, 24,25, 25,26, 26,27,
         27,28, 28,29, 29,30, 30,31, 31,32,
         32,33, 33,34, 34,35, 35,36, 36,37,
         37,38, 38,39, 39,40, 40,41, 41,42,
         42,43, 43,44, 44,45, 45,0;

    //V << -2,-2, -1.0,-2, 0,-2, 1.0,-2, 2,-2, 2,-1, 2,0, 2,1,
    //      2,2, 1,2, 0,2, -1,2, -2,2, -2,1, -2,0, -2,-1;
    //E << 0,1, 1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,10,
    //        10,11, 11,12, 12,13, 13,14, 14, 15, 15, 0;

    //Eigen::Vector2d srcP(-1.5, 1.25);

    idf::IDFdiffusion IDF;
    /*idf.computeIDF_polygon2D(V, E, srcP, 5, 0.1);
    idf.plotDiffusionMap();
    idf.plotIDF2D(30);*/

    Eigen::MatrixXd Vm;
    Eigen::MatrixXi Fm;
    IDF.triangulateFunction3D(heart, Eigen::Vector3i(128, 128, 128), 0, Vm, Fm);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(Vm, Fm);
    viewer.launch();

    Eigen::Vector3d srcP1(0.5, 1.0, 0);
    //idf.computeIDF_mesh3D(Vm, Fm, srcP1, 3, 4.0);
    //idf.computeIDF_mesh3D(sphere, srcP1, Eigen::Vector3i(35, 35, 35), 4.0, 15, 4.0);
    //idf.computeIDF_mesh3D(heart, srcP1, Eigen::Vector3i(35, 35, 35), 0.0, 10, 4.0 );
    IDF.computeIDF_mesh3D(heart, srcP1, Eigen::Vector3i(128, 128, 128), 0.0, 50, 4.0 );
    IDF.plotIDF3D(70);

    IDF.computeIDF_slice(heart, srcP1, Eigen::Vector3i(128, 128, 128), 0.0, 50, 4.0 );
    IDF.plotDiffusionMap();
    return 0;
}

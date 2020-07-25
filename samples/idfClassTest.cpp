#include "interiorDistanceField.h"
#include "frep2D.h"
#include "frep3D.h"
#include "render.h"

#include <igl/readOFF.h>
#include <igl/centroid.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/remove_duplicates.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>

#include <iostream>
#include <functional>
#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <gtsconfig.h>
#include <gts.h>

static void suriken(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj2D frep0(g.nx, g.ny, 20);

   for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
        f[i][j] = -frep0.suriken(glm::f32vec2(x, y), glm::f32vec2(0, 0));
}

static void bat(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj2D frep0(g.nx, g.ny, 20);

   for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
        f[i][j] = -frep0.bat(glm::f32vec2(x, y), glm::f32vec2(0, 0));
}

static void circle(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj2D frep0(g.nx, g.ny, 6.4);

   for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
        f[i][j] = -frep0.circle(glm::f32vec2(x, y), glm::f32vec2(0,0), 1);
}

static void heart2d(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj2D frep0(g.nx, g.ny, 6.4);

   for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
        f[i][j] = -frep0.heart2D(glm::f32vec2(x, y), glm::f32vec2(0,0));
}

static void heart3d(gdouble **f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj3D frep(g.nx, g.ny, g.nz, 60);

   for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
      for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
        f[i][j] = -frep.heart3D(glm::f32vec3(x, y, z), glm::f32vec3(0, 0, 0));
}

static void sphere(gdouble ** f, GtsCartesianGrid g, guint k, gpointer data)
{
    gdouble x, y, z = g.z;
    guint i, j;
    frep::FRepObj3D frep(g.nx, g.ny, g.nz, 6.4);

    for (i = 0, x = g.x; i < g.nx; i++, x += g.dx)
        for (j = 0, y = g.y; j < g.ny; j++, y += g.dy)
            f[i][j] = -frep.sphere(glm::f32vec3(x, y, z), glm::f32vec3(0, 0, 0), 60);
}

int main(int argc, char *argv[])
{
    //*************************************************************************
    //I. First example dedicated to computing IDF for the polygon defined using
    //   vertices and edges;
    //*************************************************************************
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

    //computing IDF
    idf::IDFdiffusion IDF;
    Eigen::Vector2d srcP0(-1.5, 1.25);
    //IDF.computeIDF_polygon2D(V, E, srcP0, 5, 0.1);
    //IDF.plotDiffusionMap();
    //IDF.plotIDF2D();

    Eigen::MatrixXd Vm_in;
    Eigen::MatrixXi Fm_in;
    igl::readOFF("bunny_small.off", Vm_in, Fm_in);
    Eigen::Vector3d srcP1(-2.0, 0.0, 0.0), center;
    igl::centroid(Vm_in, Fm_in, center);

    //IDF.computeIDF_mesh3D(10*Vm_in, Fm_in, srcP1, 50, 4.0);
    //IDF.computeIDF_slice(10*Vm_in, Fm_in, center, 70, 4.0);
    //IDF.plotDiffusionMap();
    //IDF.plotIDF3D(100);

    //******************************************************************************************
    //II. This example is dedicated to computing IDFs using function representation (FRep) in 2D;
    //******************************************************************************************

    Eigen::Vector2d srcP2(0.4, 0.0);
    IDF.computeIDF_polygon2D(suriken, Eigen::Vector3i(256, 256, 256), srcP2, 80, 0.005);
    IDF.plotDiffusionMap();
    IDF.plotIDF2D(60);

    //************************************************************************
    //III. This example is dedicated to computing of the IDFs using FRep in 3D
    //************************************************************************

    Eigen::Vector3d srcP3(0.0, 0.0, 0.0);
    //IDF.computeIDF_mesh3D(heart3d, srcP3, Eigen::Vector3i(128, 128, 128), 0.0, 140, 0.1);
    //IDF.plotDiffusionMap();
    //IDF.plotIDF3D(400);

    //IDF.computeIDF_slice(heart3d, srcP2, Eigen::Vector3i(128, 128, 128), 0.0, 82, 4.0);
    //IDF.plotIDF3D(70);

    return 0;
}

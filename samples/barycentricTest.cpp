#include "include/barycentricCoords.h"

#include <Eigen/Core>
#include <SFML/Graphics.hpp>
#include <igl/triangle/triangulate.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/isolines.h>

#include <vector>
#include <iostream>

bool findValsinMatr(Eigen::Vector2d search, Eigen::MatrixXd in)
{
    for(size_t i = 0; i < in.rows(); i++)
    {
        if(in(i, 0) == search.x() && in(i, 1) == search.y())
            return true;
    }
    return false;
}

int main(int argc, char const *argv[])
{
    // Create the boundary of a square
    Eigen::MatrixXd V2, H, Vm;
    Eigen::MatrixXi E, Fm;
    int sizeZeroLS = 12;
    V2.resize(sizeZeroLS, 2);
    E.resize(sizeZeroLS, 2);

    V2 << -2,-2, -1.2,-2, -0.5,-1, 0.5,-1, 1.2,-2, 2,-2,
           2,2,   1.2,2,   0.5,1, -0.5,1, -1.2,2, -2,2;

    E << 0,1, 1,2, 2,3, 3,4,   4,5,   5,6,
         6,7, 7,8, 8,9, 9,10, 10,11, 11,0;

    igl::triangle::triangulate(V2, E, H, "a0.005q", Vm, Fm);
    igl::opengl::glfw::Viewer viewer;

    std::vector<Eigen::Vector2d> vals;

    for(size_t i = 0; i < Vm.rows(); i++)
    {
        if(!findValsinMatr(Vm.row(i), V2))
            vals.push_back(Vm.row(i));
    }
    Eigen::MatrixXd V_inter; V_inter.resize(vals.size(), 2);
    for(size_t i = 0; i < vals.size(); i++)
        V_inter.row(i) = vals[i];

    std::vector<Eigen::Vector2d> meshPoints;
    for(size_t j = 0; j < V2.rows(); j++)
        meshPoints.push_back(V2.row(j));

    hfrep::baryCoords mvc;
    /*std::vector<Eigen::Vector2d> bary;
    for(size_t i = 0; i < V_inter.rows(); i++)
    {
        Eigen::VectorXd baryC;
        mvc.computeMVC(meshPoints, V_inter.row(i), baryC);
        bary.push_back(baryC);
    }*/

    Eigen::VectorXd baryC;
    Eigen::Vector2d p(1.5, 2.0);
    mvc.computeMVC(meshPoints, p, baryC);

    std::cout << baryC << std::endl;



    return 0;
}

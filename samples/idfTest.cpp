#include <igl/copyleft/marching_cubes.h>
#include <igl/sparse_voxel_grid.h>
#include <igl/triangle/triangulate.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>
#include <igl/isolines_map.h>

#include <Eigen/Core>

#include "dmaps/include/diffusion_map.h"
#include "dmaps/include/distance_matrix.h"
#include "dmaps/include/metrics.h"

#include "include/barycentricCoords.h"
#include "include/frep2D.h"

#include <vector>
#include <functional>
#include <iostream>
#include <matplotlib-cpp/matplotlibcpp.h>

namespace plt = matplotlibcpp;

bool findValsinMatr(Eigen::Vector2d search, Eigen::MatrixXd in)
{
    for(size_t i = 0; i < in.rows(); i++)
    {
        if(in(i, 0) == search.x() && in(i, 1) == search.y())
            return true;
    }
    return false;
}

void set_colormap(igl::opengl::glfw::Viewer & viewer)
{
  const int num_intervals = 30;
  Eigen::MatrixXd CM(num_intervals,3);
  // Colormap texture
  for(int i = 0;i<num_intervals;i++)
  {
    double t = double(num_intervals - i - 1)/double(num_intervals-1);
    CM(i,0) = std::max(std::min(2.0*t-0.0,1.0),0.0);
    CM(i,1) = std::max(std::min(2.0*t-1.0,1.0),0.0);
    CM(i,2) = std::max(std::min(6.0*t-5.0,1.0),0.0);
  }
  igl::isolines_map(Eigen::MatrixXd(CM),CM);
  viewer.data().set_colormap(CM);
}

int main(int argc, char *argv[])
{
    /*std::function<double(const Eigen::RowVector3d&)> scalar_func = [](const Eigen::RowVector3d& pt)->double
    {
        return pt.norm() - 1.0;
    };

    Eigen::RowVector3d p0(0.0, 0.0, 1.0);
    const double eps = 0.1; // parameter for construcing sparse voxel grid with size eps

    Eigen::VectorXd Cs;     // holding scalar value at each cube vertex corresponding scalar field
    Eigen::MatrixXd Cv;     // position of the corners of the sparse vocel grid
    Eigen::MatrixXi Ci;     // cubes x 8 matrix size; storing indices of 8 corners of a cube in each row

    igl::sparse_voxel_grid(p0, scalar_func, eps, 2024, Cs, Cv, Ci);
    Eigen::MatrixXi F;      // faces
    Eigen::MatrixXf V;      // vertices
    igl::copyleft::marching_cubes(Cs, Cv, Ci, V, F);*/

    //simple 2D example
    //defining frep object
   /* hfrep2D::FRepObj2D frep( 512, 512, 4.0f );
    auto fun = std::bind(&hfrep2D::FRepObj2D::heart2D, frep, std::placeholders::_1, std::placeholders::_2);
    std::vector<float> heart = frep.getFRep2D( glm::vec2( 250.0f, 150.0f ), fun );

    auto fun1 = std::bind(&hfrep2D::FRepObj2D::rectangle, frep, std::placeholders::_1, std::placeholders::_2,
                                                                std::placeholders::_3, std::placeholders::_4);
    std::vector<float> frep1H = frep.getFRep2D( glm::vec2( 180.0f, 256.0f ),
                                                       50.0f, 150.0f, fun1 );
    std::vector<float> frep2H = frep.getFRep2D( glm::vec2( 320.0f, 256.0f ),
                                                       50.0f, 150.0f, fun1 );
    std::vector<float> frep3H = frep.getFRep2D( glm::vec2( 256.0f, 256.0f ),
                                                       90.0f, 40.0f, fun1 );
    std::vector<float> frepH;
    std::vector<Eigen::Vector2d> frepZeroCs;

    for( size_t i = 0; i < 512*512; i++)
        frepH.push_back( frep.union_function(frep1H[i], frep.union_function( frep2H[i], frep3H[i], 0.0f, 0.0f ), 0.0f, 0.0f));

    //hfrep2D::render2D render;
    //sf::Image frepImg = render.drawIsolines(&frepH, 512, 512, 0.005);
    //render.displayImage(frepImg);

    glm::vec2 sLims = frep.findZeroLevelSetInterval(frepH, 70);

    //  newvalue= (max'-min')/(max-min)*(value-min)+min'.
    for(size_t y = 0; y < 512; y++)
    {
        for(size_t x = 0; x < 512; x++)
        {
            if( frepH[x+y*512] >= sLims.x && frepH[x+y*512] <= sLims.y )
            {
                Eigen::Vector2d point;
                point.x() = (8.0/512.0)*x - 4.0;
                point.y() = (8.0/512.0)*y - 4.0;
                frepZeroCs.push_back(point);
            }
        }
    }

    std::sort(frepZeroCs.begin(), frepZeroCs.end(), [](Eigen::Vector2d const &l, Eigen::Vector2d const &r)
                                                                                    { return l.x() < r.x(); });
    std::sort(frepZeroCs.begin(), frepZeroCs.end(), [](Eigen::Vector2d const &l, Eigen::Vector2d const &r)
                                                                                    { return l.y() < r.y(); });
*/

    // Create the boundary of a square
    Eigen::MatrixXd V2, H, Vm;
    Eigen::MatrixXi E, Fm;
    int sizeZeroLS = 8;
    V2.resize(sizeZeroLS, 2);
    E.resize(sizeZeroLS, 2);

    /*V2 << -2,-2, -1.2,-2, -0.5,-1, 0.5,-1, 1.2,-2, 2,-2,
           2,2,   1.2,2,   0.5,1, -0.5,1, -1.2,2, -2,2;
    E << 0,1, 1,2, 2,3, 3,4,   4,5,   5,6,
         6,7, 7,8, 8,9, 9,10, 10,11, 11,0;*/

    V2 << -2,-2, 0,-2, 2,-2, 2,0,
           2,2, 0,2, -2,2, -2,0;
    E << 0,1, 1,2, 2,3, 3,4,
         4,5, 5,6, 6,7, 7,0;

    /*int sizeZeroLS = frepZeroCs.size();
    V2.resize(sizeZeroLS, 2);
    E.resize(sizeZeroLS, 2);

    for(size_t i = 0; i < frepZeroCs.size(); i++)
    {
        V2.row(i) = frepZeroCs[i];
        if(i == frepZeroCs.size()-1)
            E.row(i) << i, 0;
        else
            E.row(i) << i, i+1;

        std::cout << V2.row(i) << std::endl;
    }*/

    dmaps::matrix_t matr = V2.cast<float>();

    //computing distances
    int num_threads = omp_get_thread_num();
    dmaps::distance_matrix distM(matr, num_threads);
    auto metrics = std::bind(&dmaps::euclidean, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    distM.compute(metrics);
    dmaps::matrix_t dist = distM.get_distances();

    //computing diffusion map and eigen vectors
    dmaps::vector_t weights;
    dmaps::diffusion_map difMap(dist, weights, num_threads);

    difMap.set_kernel_bandwidth(4.0);
    difMap.compute(4, 1.0, 0.0);
    dmaps::matrix_t eigVec = difMap.get_eigenvectors();
    dmaps::vector_t eigVal = difMap.get_eigenvalues();

    /*for(int i = 0; i < eigVec.rows(); i++)
        for(int j = 0; j < eigVec.cols(); j++)
        {
            eigVec(i, j) = eigVec(i, j) / eigVec(i,0);
        }*/

    //visualizing two eigenvectors 1st and 2nd normalized by values of the 0
    std::string answer;
    std::cout << "Show diffusion map [y/n]: ";
    std::cin >> answer;
    if(answer == "y" )
    {
        std::vector<float> vecXn, vecYn, vecX, vecY, vecZ;
        for(int i = 0; i < eigVec.rows(); i++)
        {
            vecXn.push_back(eigVec(i, 1) / eigVec(i, 0));
            vecYn.push_back(eigVec(i, 2) / eigVec(i, 0));
        }
        plt::scatter(vecXn, vecYn);
        plt::show();
    }

    //then we can use obtained eigen vectors and eigen values for computing L2 distance
    // 1st: compute boundary distances using diffusion map and computed eigen functions and eigen values
    dmaps::matrix_t D_ij;
    D_ij.resize(eigVec.rows(), eigVec.rows());
    float distSq = 0;
    float t = 1.0f/(8.0f*eigVal[1]);

    for(int i = 0; i < eigVec.rows(); i++)
        for(int j = 0; j < eigVec.rows(); j++)
        {
            for(int k = 1; k < eigVec.cols(); k++)
            {
                float eigDif = eigVec(i, k) - eigVec(j, k);
                distSq += std::exp(-2.0*eigVal[k]*t)*eigDif*eigDif;
            }
            D_ij(i, j) = distSq;
            distSq = 0;
        }

    //2nd: create grid with points in interior before computing barycentric interpolation
    igl::triangle::triangulate(V2, E, H, "a0.005q", Vm, Fm);

    //defining source point
    int sInd = 200;
    Eigen::Vector2d srcP = Vm.row(sInd);

    //5th: generating barycoords on the surface
    std::vector<Eigen::Vector2d> meshPoints;
    Eigen::VectorXd baryW1, baryW2;
    for(size_t j = 0; j < V2.rows(); j++)
        meshPoints.push_back(V2.row(j));

    //6th: compute distance in interior using D_ij and barycentric interpolation
    Eigen::VectorXd D; D.resize(Vm.rows());
    double dSum1 = 0.0, dSum2 = 0.0;
    hfrep::baryCoords mvc;

    for(size_t Jj = 0; Jj < Vm.rows(); Jj++)
    {
        for(size_t i = 0; i < D_ij.rows(); i++ )
            for(size_t j = 0; j < D_ij.cols(); j++)
            {
                mvc.computeMVC(meshPoints, srcP, baryW1);
                mvc.computeMVC(meshPoints, Vm.row(Jj), baryW2);

                dSum1 += D_ij(i, j) * baryW1[i] * baryW2[j];
                dSum2 += D_ij(i, j) * (baryW1[i] * baryW1[j] + baryW2[i] * baryW2[j]);
            }
        D[Jj] = std::sqrt(dSum1 - 0.5*dSum2);
        dSum1 = 0.0;
        dSum2 = 0.0;
    }

    Eigen::MatrixXd C, isoV;
    Eigen::MatrixXi isoE;

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(Vm, Fm);
    viewer.data().set_data(D);
    set_colormap(viewer);
    viewer.data().show_lines = false;
    viewer.launch();

    return 0;
}

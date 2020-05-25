#ifndef H_INTERIOR_DISTANCE_FIELDS_CLASS
#define H_INTERIOR_DISTANCE_FIELDS_CLASS

#include "external/dmaps/include/diffusion_map.h"
#include "external/dmaps/include/distance_matrix.h"

#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <vector>

#include <gtsconfig.h>
#include <gts.h>

namespace idf {

class IDFdiffusion
{
public:
    IDFdiffusion()
    {
#ifdef _OPENMP
    Eigen::initParallel();
#endif
    slice_z = 0.5;
    }
    ~IDFdiffusion() {}

    void computeIDF_polygon2D(const Eigen::MatrixXd &polyVerts, const Eigen::MatrixXi &meshEdges, const Eigen::Vector2d &srcP, const int eigVecNumber, double kernelBandW);
    void computeIDF_polygon2D(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes, const int eigVecNumber, double kernelBandW);

    void computeIDF_mesh3D(const Eigen::MatrixXd &meshVerts, const Eigen::MatrixXi &meshFaces, const Eigen::Vector3d &srcP, const int eigVecNumber, double kernelBandW);
    void computeIDF_mesh3D(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes, double iso,
                           const int eigVecNumber, double kernelBandW);
    void computeIDF_slice(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes, double iso,
                          const int eigVecNumber, double kernelBandW);

    void plotIDF2D(int isoNum = 30);
    void plotIDF3D(int isoNum = 30);

    void triangulateFunciton2D(GtsIsoCartesianFunc func, const Eigen::Vector3i gridRes, double iso, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
    void triangulateFunction3D(GtsIsoCartesianFunc func, const Eigen::Vector3i gridRes, double iso, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
    inline Eigen::MatrixXf getIDF() { return IDF; }

#ifdef USE_MATPLOTLIB
    void plotDiffusionMap();
#endif

private:
    void resetParams();
    void setColormap(igl::opengl::glfw::Viewer & viewer, int isoNum);
    void computeDiffusionMap(const dmaps::matrix_t &inPoints, const int eigVecNum, double kernelBandWidth );
    void computeInteriorDF2D(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXd &inVerts, const Eigen::VectorXd &srcP);
    void computeInteriorDF3D(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXd &inVerts, const Eigen::MatrixXi faces, const Eigen::VectorXd &srcP);
    void checkMeanValueCoords();
    void update_visualization(igl::opengl::glfw::Viewer &viewer);
    bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);
    void iterateIDF_slice(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXi faces, const Eigen::VectorXd &srcP );

private:
    dmaps::matrix_t eigVecs;
    dmaps::vector_t eigVals;
    dmaps::matrix_t dist, kernelM;

    Eigen::MatrixXd V, Vm, TVm;
    Eigen::MatrixXi F, Fm, TFm, Tm;
    Eigen::VectorXf IDF;

    //parameters for slicer
    Eigen::MatrixXd V_surf;
    Eigen::MatrixXi F_surf;
    Eigen::Vector3d sP;

    int isoNumb;
    double slice_z;
};

} // namespace hfrep
#endif
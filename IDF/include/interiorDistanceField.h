#ifndef H_INTERIOR_DISTANCE_FIELDS_CLASS
#define H_INTERIOR_DISTANCE_FIELDS_CLASS

/*
 * This class is based on the following papers:
 * Rustamov, R.M., Lipman, Y. and Funkhouser, T. (2009), Interior Distance Using Barycentric Coordinates.
 * Computer Graphics Forum, 28: 1279-1288. doi:10.1111/j.1467-8659.2009.01505.x
 *
 * The external library - dmaps is based on the following paper:
 * Ronald R. Coifman, Stéphane Lafon, Diffusion maps, Applied and Computational Harmonic Analysis,
 * Volume 21, Issue 1, 2006, Pages 5-30, https://doi.org/10.1016/j.acha.2006.04.006.
 *
 * GTS library is used for extracting polygonised surface of the object defined by the function
 * for its further processing in dmaps to obtain diffusion map.
 *
 * libIGL library is used for cleaning the obtained mesh after triangulation/tetrahidralisation
 * and rendering the computed IDF.
 *
 */

#include "dmaps/include/diffusion_map.h"
#include "dmaps/include/distance_matrix.h"

#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <vector>

#include <gtsconfig.h>
#include <gts.h>

namespace idf {

class IDFdiffusion
{
public:
    IDFdiffusion() : slice_z(0.5)
    {
#ifdef _OPENMP
    Eigen::initParallel();
#endif
    }
    ~IDFdiffusion() = default;

    void computeIDF_polygon2D(const Eigen::MatrixXd &polyVerts, const Eigen::MatrixXi &meshEdges, const Eigen::Vector2d &srcP, const int eigVecNumber, double kernelBandW);
    void computeIDF_polygon2D(GtsIsoCartesianFunc func, const Eigen::Vector3i resGr, const Eigen::Vector2d &srcP, const int eigVecNumber, double kernelBandW);

    void computeIDF_mesh3D(const Eigen::MatrixXd &meshVerts, const Eigen::MatrixXi &meshFaces, const Eigen::Vector3d &srcP, const int eigVecNumber, double kernelBandW);
    void computeIDF_mesh3D(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes, double iso,
                           const int eigVecNumber, double kernelBandW);
    void computeIDF_slice(const Eigen::MatrixXd &meshVerts, const Eigen::MatrixXi &meshFaces, const Eigen::Vector3d &srcP,
                          const int eigVecNumber, double kernelBandW);
    void computeIDF_slice(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes, double iso,
                          const int eigVecNumber, double kernelBandW);

    void plotIDF2D(int isoNum = 30);
    void plotIDF3D(int isoNum = 30);

    void getSurfaceComponents2D(GtsIsoCartesianFunc func, const Eigen::Vector3i gridRes,
                                double iso, Eigen::MatrixXd &V, Eigen::MatrixXi &E);
    void getSurfaceComponents3D(GtsIsoCartesianFunc func, const Eigen::Vector3i gridRes,
                                double iso, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

    inline Eigen::VectorXf getIDF(){ return IDF; }

#ifdef USE_MATPLOTLIB
    void plotDiffusionMap();
#endif

private:
    void resetParams();
    void computeDiffusionMap(const dmaps::matrix_t &inPoints, const int eigVecNum, double kernelBandWidth );
    void computeInteriorDF2D(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXd &inVerts, const Eigen::VectorXd &srcP);
    void computeInteriorDF3D(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXd &inVerts, const Eigen::MatrixXi faces, const Eigen::VectorXd &srcP);
    void iterateIDF_slice(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXi faces, const Eigen::VectorXd &srcP );
    void autoSelectEigenVectors(Eigen::VectorXf &eigVals_s, Eigen::MatrixXf &eigVecs_s, const float t);
    void update_visualization(igl::opengl::glfw::Viewer &viewer);
    void setColormap(igl::opengl::glfw::Viewer & viewer, int isoNum);
    bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);

    Eigen::MatrixXf computePairwiseDist();

private:
    dmaps::matrix_t eigVecs;
    dmaps::vector_t eigVals;
    dmaps::matrix_t dist, kernelM;

    Eigen::MatrixXd V, Vm, TVm;
    Eigen::MatrixXi F, Fm, TFm, Tm, Em;
    Eigen::VectorXf IDF;

    //parameters for slicer
    Eigen::MatrixXd V_surf;
    Eigen::MatrixXi F_surf;
    Eigen::Vector3d sP;

    int isoNumb;
    double slice_z;
};

} // namespace idf
#endif

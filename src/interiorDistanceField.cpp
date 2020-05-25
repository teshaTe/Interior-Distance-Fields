#include "include/interiorDistanceField.h"
#include "include/barycentricCoords.h"
#include "include/timer.hpp"

#include "external/dmaps/include/metrics.h"

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/remove_duplicates.h>
#include <igl/triangle/triangulate.h>
#include <igl/isolines_map.h>
#include <igl/barycenter.h>
#include <igl/marching_tets.h>
#include <igl/parula.h>

#include <gts.h>
#include <iostream>
#include <functional>
#include <omp.h>


#ifdef USE_MATPLOTLIB
    #include <matplotlib-cpp/matplotlibcpp.h>
#endif

namespace idf {

void IDFdiffusion::computeIDF_polygon2D(const Eigen::MatrixXd &polyVerts, const Eigen::MatrixXi &meshEdges,
                                        const Eigen::Vector2d &srcP, const int eigVecNumber, double kernelBandW)
{
    assert(polyVerts.size() > 0);
    assert(meshEdges.size() > 0);

    Eigen::MatrixXd H;
    resetParams();
    computeDiffusionMap(polyVerts.cast<float>(), eigVecNumber, kernelBandW);
    igl::triangle::triangulate(polyVerts, meshEdges, H, "a0.005q", Vm, Fm);
    computeInteriorDF2D(polyVerts, Vm, srcP);
}

void IDFdiffusion::computeInteriorDF2D(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXd &inVerts, const Eigen::VectorXd &srcP)
{
    assert(eigVecs.size() > 0);
    assert(eigVals.size() > 0);

    std::cout << "Stage: starting computing IDF." << std::endl;

    //computing distances on the boundary of the mesh/polygon using diffusion map;
    Eigen::MatrixXf D_ij; D_ij.resize(eigVecs.rows(), eigVecs.rows());
    float t =1.0f/ (8.0f * eigVals[1]);
    float dSq = 0.0f;

    std::cout << "Stage: starting computing diffusion distances on the boundary." << std::endl;
    for(int i = 0; i < eigVecs.rows(); i++)
        for(int j = 0; j < eigVecs.rows(); j++)
        {
            for(int k = 1; k < eigVecs.cols(); k++)
            {
                float eigDiff = eigVecs(i, k) - eigVecs(j , k);
                dSq += std::pow(eigVals[k], 2.0f * t) * eigDiff * eigDiff;
            }
            D_ij(i, j) = dSq;
            dSq = 0.0f;
        }
    std::cout << "Stage: finished.\n" << std::endl;

    std::vector<Eigen::Vector2d> meshPoints;
    for(size_t i = 0; i < surfMeshV.rows(); i++)
        meshPoints.push_back(surfMeshV.row(i));

    float dSum1 = 0.0f, dSum2 = 0.0f;
    IDF.resize(inVerts.rows());

    //computing mean-value coordinates; 1st: compute them for the source point srcP;
    idf::baryCoords mvc;
    Eigen::VectorXf baryW1 = mvc.meanValueCoords2D(meshPoints, srcP);

    std::cout << "\nStage: starting computing mean value interpoaltion." << std::endl;
    std::cout << "Total points to process: " << inVerts.rows() << std::endl;
    prof::timer time;
    Eigen::VectorXf baryW2;
    //2nd: computing mean value coords for the rest interior points and points along the boundary
    for(size_t l = 0; l < inVerts.rows(); l++)
    {
        time.Start();
        baryW2 = mvc.meanValueCoords2D(meshPoints, inVerts.row(l));

        for(int i = 0; i < D_ij.rows(); i++)
            for(int j = 0; j < D_ij.cols(); j++)
            {
                dSum1 += D_ij(i, j) * baryW1[i] * baryW2[j];
                dSum2 += D_ij(i, j) * (baryW1[i] * baryW1[j] + baryW2[i] * baryW2[j]);
            }
        IDF[l] = std::sqrt(dSum1 - 0.5f * dSum2);
        dSum1 = dSum2 = 0.0f;
        time.End("1 point: ");
    }
    std::cout << "\nStage: finished.\n" << std::endl;
}

void IDFdiffusion::plotIDF2D(int isoNum)
{
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(Vm, Fm);
    viewer.data().set_data(IDF.cast<double>());
    setColormap(viewer, isoNum);
    viewer.data().show_lines = false;
    viewer.launch();
}

void IDFdiffusion::computeIDF_mesh3D(const Eigen::MatrixXd &meshVerts, const Eigen::MatrixXi &meshFaces, const Eigen::Vector3d &srcP, const int eigVecNumber, double kernelBandW)
{
    resetParams();
    computeDiffusionMap(meshVerts.cast<float>(), eigVecNumber, kernelBandW);
    igl::copyleft::tetgen::tetrahedralize(meshVerts, meshFaces, "pq1.414Y", TVm, Tm, TFm);
    computeInteriorDF3D(meshVerts, TVm, meshFaces, srcP);
}

void IDFdiffusion::computeIDF_mesh3D(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes,
                                     double iso, const int eigVecNumber, double kernelBandW)
{
    resetParams();
    triangulateFunction3D(func, gridRes, iso, V, F);
    Eigen::VectorXi I;
    igl::remove_duplicates(V, F, Vm, Fm, I);

    computeDiffusionMap(Vm.cast<float>(), eigVecNumber, kernelBandW);
    igl::copyleft::tetgen::tetrahedralize(Vm, Fm, "pq1.414Y", TVm, Tm, TFm);
    computeInteriorDF3D(Vm, TVm, Fm, srcP);
}

void IDFdiffusion::computeIDF_slice(GtsIsoCartesianFunc func, const Eigen::Vector3d &srcP, const Eigen::Vector3i gridRes, double iso, const int eigVecNumber, double kernelBandW)
{
    resetParams();
    triangulateFunction3D(func, gridRes, iso, V, F);
    Eigen::VectorXi I;
    igl::remove_duplicates(V, F, Vm, Fm, I);

    computeDiffusionMap(Vm.cast<float>(), eigVecNumber, kernelBandW);
    igl::copyleft::tetgen::tetrahedralize(Vm, Fm, "pq1.414Y", TVm, Tm, TFm);
    iterateIDF_slice(Vm, Fm, srcP);
}

void IDFdiffusion::iterateIDF_slice(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXi faces, const Eigen::VectorXd &srcP)
{
    V_surf = surfMeshV;
    F_surf = faces;
    sP = srcP;

    igl::opengl::glfw::Viewer viewer;
    update_visualization(viewer);
    viewer.callback_key_down = std::bind(&IDFdiffusion::key_down, this, std::placeholders::_1,
                                                std::placeholders::_2, std::placeholders::_3);
    viewer.launch();
}

void IDFdiffusion::computeInteriorDF3D(const Eigen::MatrixXd &surfMeshV, const Eigen::MatrixXd &inVerts,
                                       const Eigen::MatrixXi faces, const Eigen::VectorXd &srcP)
{
    assert(eigVecs.size() > 0);
    assert(eigVals.size() > 0);

    std::cout << "Stage: starting computing IDF." << std::endl;

    //computing distances on the boundary of the mesh/polygon using diffusion map;
    Eigen::MatrixXf D_ij; D_ij.resize(eigVecs.rows(), eigVecs.rows());
    float t = 1.0f / (8.0f * eigVals[1]);
    float dSq = 0.0f;

    std::cout << "Stage: starting computing diffusion distances on the boundary." << std::endl;
    for(int i = 0; i < eigVecs.rows(); i++)
        for(int j = 0; j < eigVecs.rows(); j++)
        {
            for(int k = 1; k < eigVecs.cols(); k++)
            {
                float eigDiff = eigVecs(i, k) - eigVecs(j , k);
                dSq += std::pow(eigVals[k], 2.0f * t) * eigDiff * eigDiff;
            }
            D_ij(i, j) = dSq;
            dSq = 0.0f;
        }
    std::cout << "Stage: finished.\n" << std::endl;

    std::vector<Eigen::Vector3d> meshPoints;
    for(size_t i = 0; i < surfMeshV.rows(); i++)
        meshPoints.push_back(surfMeshV.row(i));

    IDF.resize(inVerts.rows());

    //computing mean-value coordinates; 1st: compute them for the source point srcP;
    idf::baryCoords mvc;
    Eigen::VectorXf baryW1 = mvc.meanValueCoords3D(meshPoints, faces, srcP);

    std::cout << "\nStage: starting computing mean value interpoaltion." << std::endl;
    std::cout << "Total points to process: " << inVerts.rows() << std::endl;

    prof::timer time;
    Eigen::VectorXf baryW2;
    float dSum1 = 0.0f, dSum2 = 0.0f;

    //2nd: computing mean value coords for the rest interior points and points along the boundary
    for(size_t l = 0; l < inVerts.rows(); l++)
    {
        time.Start();
        baryW2 = mvc.meanValueCoords3D(meshPoints, faces, inVerts.row(l));
#ifdef _OPENMP
#pragma omp parallel for reduction(+:dSum1, dSum2) shared(baryW1, D_ij, baryW2) schedule(static)
#endif
        for(int i = 0; i < D_ij.rows(); i++)
            for(int j = 0; j < D_ij.cols(); j++)
            {
                dSum1 += D_ij(i, j) * baryW1[i] * baryW2[j];
                dSum2 += D_ij(i, j) * (baryW1[i] * baryW1[j] + baryW2[i] * baryW2[j]);
            }
        IDF[l] = std::sqrt(dSum1 - 0.5f * dSum2);
        dSum1 = dSum2 = 0.0f;

        time.End("1 point: ");
    }

    std::cout << "\nStage: finished. \n" << std::endl;
}

static void pick_first_face (GtsFace * f, GtsFace ** first)
{
  if (*first == NULL)
    *first = f;
}

void IDFdiffusion::triangulateFunction3D(GtsIsoCartesianFunc func, const Eigen::Vector3i gridRes,
                                         double iso, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    GtsCartesianGrid g;
    GtsSurface * surface;

    g.nx = gridRes.x();
    g.ny = gridRes.y();
    g.nz = gridRes.z();

    /* interval is [-10:10][-10:10][-10:10] */
    g.x = -10.0; g.dx = 20./(gdouble) (g.nx - 1);
    g.y = -10.0; g.dy = 20./(gdouble) (g.ny - 1);
    g.z = -10.0; g.dz = 20./(gdouble) (g.nz - 1);

    surface = gts_surface_new (gts_surface_class (),
                               gts_face_class (),
                               gts_edge_class (),
                               gts_vertex_class ());
    //gts_isosurface_tetra (surface, g, func, NULL, iso);
    gts_isosurface_cartesian (surface, g, func, NULL, iso);
    gts_surface_print_stats (surface, stderr);

    GtsFace *first = NULL;
    gts_surface_foreach_face (surface, (GtsFunc) pick_first_face, &first);
    GtsRange depth_range; gts_range_init (&depth_range);
    int j = 0;

    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3i> faces;

    if (first)
    {
        GtsSurfaceTraverse * t = gts_surface_traverse_new(surface, first);
        GtsFace * f;
        guint level;
        while ((f = gts_surface_traverse_next (t, &level)))
        {
            GtsVertex *v1, *v2, *v3;
            gts_triangle_vertices(GTS_TRIANGLE(f), &v1, &v2, &v3);

            verts.push_back(Eigen::Vector3d(v1->p.x, v1->p.y, v1->p.z));
            verts.push_back(Eigen::Vector3d(v2->p.x, v2->p.y, v2->p.z));
            verts.push_back(Eigen::Vector3d(v3->p.x, v3->p.y, v3->p.z));

            faces.push_back(Eigen::Vector3i(j, j+1, j+2));
            j += 3;
            gts_range_add_value (&depth_range, level);
        }
        gts_surface_traverse_destroy (t);
    }

    V.resize(verts.size(), 3);
    F.resize(faces.size(), 3);

    for(int i = 0; i < faces.size(); i++)
        F.row(i) = faces[i];

    for(int i = 0; i < verts.size(); i++)
        V.row(i) = verts[i];
}

void IDFdiffusion::computeDiffusionMap(const dmaps::matrix_t &inPoints, const int eigVecNum, double kernelBandWidth)
{
    std::cout << "Stage: starting computing diffusion map." << std::endl;
    //computing distances
    int num_threads = omp_get_num_threads();
    dmaps::distance_matrix dMatr(inPoints, num_threads);
    auto metrics = std::bind(&dmaps::euclidean, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    dMatr.compute(metrics);
    dist = dMatr.get_distances();

    //computing diffusion map: eigen values and eigen functions of the Laplace-Beltrammi operator
    dmaps::vector_t weights;
    dmaps::diffusion_map diffMap(dist, weights, num_threads);
    diffMap.set_kernel_bandwidth(kernelBandWidth);
    diffMap.compute(eigVecNum, 1.0, 0.0);

    eigVals = diffMap.get_eigenvalues();
    eigVecs = diffMap.get_eigenvectors();
    std::cout << "Stage: finished.\n" << std::endl;
}

void IDFdiffusion::resetParams()
{
    dist.resize(0, 0);    dist.setZero();
    eigVals.resize(0);    eigVals.setZero();
    eigVecs.resize(0, 0); eigVecs.setZero();
    IDF.resize(0);        IDF.setZero();
    Vm.resize(0, 0);      Vm.setZero();
    Fm.resize(0, 0);      Fm.setZero();
    Tm.resize(0, 0);      Tm.setZero();
}

void IDFdiffusion::plotIDF3D(int isoNum)
{
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(TVm, TFm);
    viewer.data().set_data(IDF.cast<double>());
    setColormap(viewer, isoNum);
    viewer.data().show_faces = true;
    viewer.launch();
}

void IDFdiffusion::setColormap(igl::opengl::glfw::Viewer &viewer, int isoNum)
{
    int num_intervals = isoNum;
    Eigen::MatrixXd CM(num_intervals, 3);
    // Colormap texture
    for(int i = 0; i<num_intervals; i++)
    {
        double t = double(num_intervals - i - 1)/double(num_intervals - 1);
        CM(i, 0) = std::max(std::min(2.0 * t - 0.0, 1.0), 0.0);
        CM(i, 1) = std::max(std::min(2.0 * t - 1.0, 1.0), 0.0);
        CM(i, 2) = std::max(std::min(6.0 * t - 5.0, 1.0), 0.0);
    }
    igl::isolines_map(Eigen::MatrixXd(CM), CM);
    viewer.data().set_colormap(CM);
}

void IDFdiffusion::update_visualization(igl::opengl::glfw::Viewer &viewer)
{
    Eigen::Vector4d plane(0, 0, 1,-((1 - slice_z) * TVm.col(2).minCoeff() + slice_z * TVm.col(2).maxCoeff()));
    Eigen::MatrixXd V_vis;
    Eigen::MatrixXi F_vis;

    Eigen::VectorXi J;
    {
        Eigen::SparseMatrix<double> bary;
        const Eigen::VectorXd IV = ( TVm.col(0)*plane(0) +
                                     TVm.col(1)*plane(1) +
                                     TVm.col(2)*plane(2)).array() + plane(3);

        igl::marching_tets(TVm, Tm, IV, V_vis, F_vis, J, bary);
    }

    Eigen::MatrixXd C_vis;

    computeInteriorDF3D(V_surf, V_vis, F_surf, sP);
    setColormap(viewer, 50);

    viewer.data().clear();
    viewer.data().set_mesh(V_vis, F_vis);
    viewer.data().set_data(IDF.cast<double>());
    viewer.data().set_face_based(true);
}

bool IDFdiffusion::key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
    switch(key)
      {
        default:
          return false;
        case '.':
          slice_z = std::min(slice_z+0.01,0.99);
          break;
        case ',':
          slice_z = std::max(slice_z-0.01,0.01);
          break;
      }
      update_visualization(viewer);
      return true;
}

#ifdef USE_MATPLOTLIB
void IDFdiffusion::plotDiffusionMap()
{
    std::vector<float> vecXn, vecYn;
    for(int i = 0; i < eigVecs.rows(); i++)
    {
        vecXn.push_back(eigVecs(i, 1) / eigVecs(i, 0));
        vecYn.push_back(eigVecs(i, 2) / eigVecs(i, 0));
    }
    matplotlibcpp::scatter(vecXn, vecYn);
    matplotlibcpp::show();
}
#endif

} // namespace hfrep

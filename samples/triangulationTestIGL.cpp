#include <igl/copyleft/marching_cubes.h>
#include <igl/sparse_voxel_grid.h>
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <iostream>

int main(int argc, char *argv[])
{
    std::function<double(const Eigen::RowVector3d&)> scalar_func = [](const Eigen::RowVector3d& pt)->double
    {
        return pt.norm() - 1.0;
    };

    Eigen::RowVector3d p0(0.0, 0.0, 1.0);
    const double eps = 0.1; // parameter for construcing sparse voxel grid with size eps

    Eigen::VectorXd Cs;     // holding scalar value at each cube vertex corresponding scalar field
    Eigen::MatrixXd Cv;     // position of the corners of the sparse vocel grid
    Eigen::MatrixXi Ci;     // cubes x 8 matrix size; storing indices of 8 corners of a cube in each row

    igl::sparse_voxel_grid(p0, scalar_func, eps, 1024, Cs, Cv, Ci);
    Eigen::MatrixXi F;      // faces
    Eigen::MatrixXd V;      // vertices
    igl::copyleft::marching_cubes(Cs, Cv, Ci, V, F);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.launch();

    return 0;
}

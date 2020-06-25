#ifndef H_BARYCENTRIC_COORDS_CLASS
#define H_BARYCENTRIC_COORDS_CLASS

#include <Eigen/Core>
#include <vector>

namespace idf
{
class baryCoords
{
    using eigVector_t = Eigen::Matrix<float, Eigen::Dynamic, 1>;
    using eigMatrix_t = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

public:
    baryCoords() { eps = 1.0e-8f; }
    ~baryCoords(){}

    std::vector<Eigen::Vector2d> getTriangleVerts(Eigen::MatrixXd V, Eigen::Vector3i F);
    std::vector<Eigen::Vector2d> uniformPointDistrTriangle(Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3, int nP, int &seed);

    Eigen::VectorXf meanValueCoords2D(const std::vector<Eigen::Vector2d> &polyCoords, const Eigen::Vector2d &p );
    Eigen::VectorXf meanValueCoords3D(const std::vector<Eigen::Vector3d> &polyCoords,  const Eigen::MatrixXi &faces, const Eigen::Vector3d &p);

private:
    float eps;
    bool computeBoundaryCoordinates2D(const std::vector<Eigen::Vector2d> &polyCoords, const Eigen::Vector2d &p, Eigen::VectorXf &baryCoords);
    void r8vec_uniform_01(int nP, int &seed, eigVector_t &v);
};

} //namespace hfrep
#endif

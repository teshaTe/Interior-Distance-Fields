#include "include/barycentricCoords.h"

#include <Eigen/LU>

#include <cstdlib>
#include <ctime>
#include <iostream>

namespace idf {

std::vector<Eigen::Vector2d> baryCoords::getTriangleVerts(Eigen::MatrixXd V, Eigen::Vector3i F)
{
    std::vector<Eigen::Vector2d> triangleVert;
    for(size_t j = 0; j < 3; j++)
    {
        int vInd = F[j];
        Eigen::Vector2d p = V.row(vInd);
        triangleVert.push_back(p);
    }
    return triangleVert;
}

Eigen::VectorXf baryCoords::meanValueCoords2D(const std::vector<Eigen::Vector2d> &polyCoords, const Eigen::Vector2d &p)
{
    Eigen::VectorXf baryCoords;
    size_t size = polyCoords.size();
    baryCoords.resize(size);
    baryCoords.setZero();

    //computing boundary coordinates
    if(computeBoundaryCoordinates2D(polyCoords, p, baryCoords)) return baryCoords;

    //computing point in interior of the polygon
    eigMatrix_t s; s.resize(polyCoords.size(), 2);
    eigVector_t ri; ri.resize(size);

    for(size_t i = 0; i < size; i++)
    {
        s.row(i) = (polyCoords[i] - p).cast<float>();
        ri[i] = std::sqrt(s.row(i).squaredNorm());
    }

    //check boundary conditions according to the article:
    //"Mean Value Coordinates for Arbitrary Planar Polygons", Kai Hormann and Michael Floater;
    //algorithm in Figure 10, p. 12

    eigVector_t Ai; Ai.resize(size);
    eigVector_t Di; Di.resize(size);
    eigVector_t Ti; Ti.resize(size);

    for(size_t i = 0; i < size; i++)
    {
        const size_t iP = (i+1) % size;
        Ai[i] = 0.5f * (s(i, 0) * s(iP, 1) - s(iP, 0) * s(i, 1)); // square area of the triangle through vector multiplication and its determinant def
        Di[i] = s.row(i).dot(s.row(iP)); //dot product of two vectors

        //computing tan(a_i/2) = sin(a_i) / (1 + cos(a_i));
        //assert(std::abs(ri[i] * ri[iP] + Di[i]) > 0.0);
        Ti[i] = Ai[i] / (ri[i] * ri[iP] + Di[i]);
        assert(std::abs(Ai[i]) > 0.0f);
        //Ti[i] = (ri[i] * ri[iP] - Di[i]) / Ai[i];
    }

    eigVector_t w; w.resize(size);
    float W = 0.0f;

    //computing barycentric weights
    for(size_t i = 0; i < size; i++)
    {
        const size_t iM = (i + size - 1) % size;
        assert(std::abs(ri[i]) > 0.0f);
        w[i] = (Ti[iM] + Ti[i]) / ri[i];
        W += w[i];
    }

    //computing barycentric coordinates
    assert(std::abs(W) > 0.0f);
    const float invW = 1.0 / W;

    for(size_t i = 0; i < size; i++)
        baryCoords[i] = w[i] * invW; // equation 13 from the paper, p 10

    return baryCoords;
}

bool baryCoords::computeBoundaryCoordinates2D(const std::vector<Eigen::Vector2d> &polyCoords, const Eigen::Vector2d &p,
                                              Eigen::VectorXf &baryCoords)
{
    size_t size = polyCoords.size();
    baryCoords.resize(size);
    baryCoords.setZero();

    eigMatrix_t s; s.resize(size, 2);
    eigVector_t ri; ri.resize(size);
    for(size_t i = 0; i < size; i++)
    {
        s.row(i) = (polyCoords[i] - p).cast<float>();
        ri[i] = std::sqrt(s.row(i).squaredNorm());
        if(std::abs(ri[i]) < eps)
        {
            baryCoords[i] = 1.0f;
            return true;
        }
    }

    eigVector_t Ai; Ai.resize(size);
    eigVector_t Di; Di.resize(size);

    for(size_t i = 0; i < size; i++)
    {
        const size_t iP   = (i + 1) % size;
        Ai[i] = 0.5f * (s(i, 0) * s(iP, 1) - s(i, 1) * s(iP, 0)); // square area of the triangle through vector multiplication and its determinant def
        Di[i] = s.row(i).dot(s.row(iP)); //s(iP, 0)*s(i, 0) + s(iP, 1)*s(i, 1); //dot product of two vectors

        if(std::abs(Ai[i]) < eps && Di[i] < 0.0f)
        {
            Eigen::Vector2f s1 = (p - polyCoords[iP]).cast<float>();
            Eigen::Vector2f s2 = (polyCoords[i] - polyCoords[iP]).cast<float>();

            assert(std::abs(s2.squaredNorm()) > 0.0f);

            const float opScalar = s1.dot(s2);
            const float b1 = opScalar / s2.squaredNorm();

            //storing weighting coefficients for linear interpolation along polygon edges
            baryCoords[i] = b1;
            baryCoords[iP] = 1.0f - b1;
            return true;
        }
    }
    return false;
}

Eigen::VectorXf baryCoords::meanValueCoords3D(const std::vector<Eigen::Vector3d> &polyCoords,
                                              const Eigen::MatrixXi &faces, const Eigen::Vector3d &p)
{
    size_t size = polyCoords.size();
    Eigen::VectorXf baryCoords;
    baryCoords.resize(size);
    baryCoords.setZero();

    eigMatrix_t s; s.resize(size, 3);
    eigVector_t ri; ri.resize(size);

    for(size_t i = 0; i < size; i++)
    {
        s.row(i) = (polyCoords[i] - p).cast<float>();
        ri[i] = std::sqrt(s.row(i).squaredNorm());

        if(std::abs(ri[i]) < eps)
        {
             baryCoords[i] = 1.0f;
             return baryCoords;
        }
        //assert(std::abs(ri[i]) > 0.0f);
        s.row(i) /= ri[i];
    }

    std::vector<float> theta0, theta1, theta2, halfSum;
    theta0.resize(faces.rows());
    theta1.resize(faces.rows());
    theta2.resize(faces.rows());
    halfSum.resize(faces.rows());

    std::vector<int> pid0, pid1, pid2;
    pid0.resize(faces.rows());
    pid1.resize(faces.rows());
    pid2.resize(faces.rows());

    //cycle over all triangles to compute weights
    //1st: compute all angles

#ifdef _OPENMP
#pragma omp parallel for simd schedule(static)
#endif
    for(int i = 0; i < faces.rows(); i++)
    {
        //getting vertex ids
        pid0[i] = faces.row(i).x();
        pid1[i] = faces.row(i).y();
        pid2[i] = faces.row(i).z();
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int i = 0; i < faces.rows(); i++)
    {
        //obtaining unit vectors using ids
        Eigen::Vector3f u0 = s.row(pid0[i]);
        Eigen::Vector3f u1 = s.row(pid1[i]);
        Eigen::Vector3f u2 = s.row(pid2[i]);

        //computing edge length
        float l0 = std::sqrt((u1 - u2).squaredNorm());
        float l1 = std::sqrt((u2 - u0).squaredNorm());
        float l2 = std::sqrt((u0 - u1).squaredNorm());

        //computing angles
        theta0[i] = 2.0f * std::asin(l0 / 2.0f);
        theta1[i] = 2.0f * std::asin(l1 / 2.0f);
        theta2[i] = 2.0f * std::asin(l2 / 2.0f);
        halfSum[i] = (theta0[i] + theta1[i] + theta2[i]) / 2.0f;
    }

    //2nd: compute border baryCoords
    for(int i = 0; i < faces.rows(); i++)
    {
        if(M_PI - halfSum[i] < eps)
        {
            baryCoords.resize(size);
            baryCoords.setZero();

            baryCoords[pid0[i]] = std::sin(theta0[i]) * ri[pid1[i]] * ri[pid2[i]];
            baryCoords[pid1[i]] = std::sin(theta1[i]) * ri[pid2[i]] * ri[pid0[i]];
            baryCoords[pid2[i]] = std::sin(theta2[i]) * ri[pid0[i]] * ri[pid1[i]];
            float sumW = baryCoords[pid0[i]] + baryCoords[pid1[i]] + baryCoords[pid2[i]];

            baryCoords[pid0[i]] /= sumW;
            baryCoords[pid1[i]] /= sumW;
            baryCoords[pid2[i]] /= sumW;
            return baryCoords;
        }
    }

    //3rd: computing the rest
#ifdef _OPENMP
#pragma omp parallel for schedule(static) shared(theta0, theta1, theta2, halfSum)
#endif
    for(int i = 0; i < faces.rows(); i++)
    {
        // coefficient
        float sinHalfSum = std::sin(halfSum[i]);
        float sinHalfSumSubTheta0 = std::sin(halfSum[i] - theta0[i]);
        float sinHalfSumSubTheta1 = std::sin(halfSum[i] - theta1[i]);
        float sinHalfSumSubTheta2 = std::sin(halfSum[i] - theta2[i]);
        float sinTheta0 = std::sin(theta0[i]),
              sinTheta1 = std::sin(theta1[i]),
              sinTheta2 = std::sin(theta2[i]);

        //assert(std::abs(sinTheta0) > 0.0f);
        //assert(std::abs(sinTheta1) > 0.0f);
       // assert(std::abs(sinTheta2) > 0.0f);

        float c0 = 2.0f * sinHalfSum * sinHalfSumSubTheta0 / (sinTheta1 * sinTheta2) - 1.0f;
        float c1 = 2.0f * sinHalfSum * sinHalfSumSubTheta1 / (sinTheta2 * sinTheta0) - 1.0f;
        float c2 = 2.0f * sinHalfSum * sinHalfSumSubTheta2 / (sinTheta0 * sinTheta1) - 1.0f;

        if(std::abs(c0) > 1.0f) c0 = c0 > 0 ? 1 : -1;
        if(std::abs(c1) > 1.0f) c1 = c1 > 0 ? 1 : -1;
        if(std::abs(c2) > 1.0f) c2 = c2 > 0 ? 1 : -1;

        //checking sign of the determinant of three unit vectors
        Eigen::Matrix3f uniMatr;
        uniMatr.row(0) = s.row(pid0[i]);
        uniMatr.row(1) = s.row(pid1[i]);
        uniMatr.row(2) = s.row(pid2[i]);
        float det = uniMatr.determinant();

        if(std::abs(det) < eps)
        {
            i++; continue;
        }

        float detSign = det > 0 ? 1 : -1;
        float sign0 = detSign * std::sqrt(1.0f - c0*c0);
        float sign1 = detSign * std::sqrt(1.0f - c1*c1);
        float sign2 = detSign * std::sqrt(1.0f - c2*c2);

        // if 'p' lies on the plane of current triangle but outside it, ignore the current triangle
        if (std::abs(sign0) < eps || std::abs(sign1) < eps || std::abs(sign2) < eps)
        {
            i++; continue;
        }

        // weight
        baryCoords[pid0[i]] += (theta0[i] - c1*theta2[i] - c2*theta1[i]) / (ri[pid0[i]] * sinTheta1 * sign2);
        baryCoords[pid1[i]] += (theta1[i] - c2*theta0[i] - c0*theta2[i]) / (ri[pid1[i]] * sinTheta2 * sign0);
        baryCoords[pid2[i]] += (theta2[i] - c0*theta1[i] - c1*theta0[i]) / (ri[pid2[i]] * sinTheta0 * sign1);
    }

    // normalize weight
    float sumWeight = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sumWeight)
#endif
    for (int i = 0; i < size; ++i)
        sumWeight += baryCoords[i];

    if(!sumWeight)
        std::cerr << "WARNING: zero weights." << std::endl;

#ifdef _OPENMP
#pragma omp parallel for simd
#endif
    for (int i = 0; i < size; ++i)
        baryCoords[i] /= sumWeight;

    return baryCoords;
}

}//namespace hfrep

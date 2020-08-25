#include "frep3D.h"
#include <algorithm>

namespace frep {

FRepObj3D::FRepObj3D(int resX, int resY, int resZ, float scaleF) : resolutionX(resX),
                                                                   resolutionY(resY),
                                                                   resolutionZ(resZ),
                                                                   scale(scaleF) {}

float FRepObj3D::sphere(glm::f32vec3 pos, glm::f32vec3 center, float R)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);

    float shX = scale * (str.x - c0.x);
    float shY = scale * (str.y - c0.y);
    float shZ = scale * (str.z - c0.z);

    return rad*rad - shX*shX - shY*shY - shZ*shZ;
}

float FRepObj3D::blobby(glm::f32vec3 pos, glm::f32vec3 center, float R)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return -(shX*shX + shY*shY +shZ*shZ + std::sin(4.0f*shX) + std::sin(4.0f*shY) + std::sin(4.0f*shZ) - rad);
}

float FRepObj3D::cylinderX(glm::f32vec3 pos, glm::f32vec3 center, float R, float h)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);
    float H = convertToSTR(h);

    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    float infCylX = rad*rad - shY*shY - shZ*shZ;
    float cyl = intersect_function(infCylX, c0.y + H/2.0f);
    return cyl;
}

float FRepObj3D::cylinderY(glm::f32vec3 pos, glm::f32vec3 center, float R, float h)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);
    float H = convertToSTR(h);

    float shX = str.x - c0.x;
    float shZ = str.z - c0.z;

    float infCylY = rad*rad - shX*shX - shZ*shZ;
    float cyl = intersect_function(c0.y-H/2.0f, infCylY, 0.0, 1.0);

    return cyl;
}

float FRepObj3D::cylinderZ(glm::f32vec3 pos, glm::f32vec3 center, float R, float h)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);
    float H = convertToSTR(h);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;

    float infCylZ = rad*rad - shX*shX - shY*shY;
    float cyl = intersect_function(c0.z - H/2.0f, infCylZ, 0.0, 1.0);

    return cyl;
}

float FRepObj3D::ellipticCylinderX(glm::f32vec3 pos, glm::f32vec3 center, float a, float b)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return 1.0f/resolutionX - shY*shY/(a*a) + shZ*shZ/(b*b);
}

float FRepObj3D::ellipticCylinderY(glm::f32vec3 pos, glm::f32vec3 center, float a, float b)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shZ = str.z - c0.z;

    return 1.0f/resolutionX - shX*shX/(a*a) - shZ*shZ/(b*b);
}

float FRepObj3D::ellipticCylinderZ(glm::f32vec3 pos, glm::f32vec3 center, float a, float b)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;

    return 1.0f/resolutionX - shX*shX/(a*a) - shY*shY/(b*b);
}

float FRepObj3D::coneX(glm::f32vec3 pos, glm::f32vec3 center, float R, float h)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float height = convertToSTR(h);
    float rad    = convertToSTR(R);

    float shX = scale * (str.x - c0.x);
    float shY = scale * (str.y - c0.y);
    float shZ = scale * (str.z - c0.z);
    return shX*shX - (shY*shY + shZ*shZ)/(rad*rad);
}

float FRepObj3D::coneY(glm::f32vec3 pos, glm::f32vec3 center, float R)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return shY*shY - shX*shX/(rad*rad) - shZ*shZ/(rad*rad);
}

float FRepObj3D::coneZ(glm::f32vec3 pos, glm::f32vec3 center, float R)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rad   = convertToSTR(R);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return shZ*shZ - shX*shX/(rad*rad) - shY*shY/(rad*rad);
}

float FRepObj3D::ellipticConeX(glm::f32vec3 pos, glm::f32vec3 center, float a, float b)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return shX*shX - shY*shY/(a*a) - shZ*shZ/(b*b);
}

float FRepObj3D::ellipticConeY(glm::f32vec3 pos, glm::f32vec3 center, float h, float a, float b)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return shY*shY - shX*shX/(a*a) - shZ*shZ/(b*b);
}

float FRepObj3D::ellipticConeZ(glm::f32vec3 pos, glm::f32vec3 center, float h, float a, float b)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return shZ*shZ - shX*shX/(a*a) - shY*shY/(b*b);
}

float FRepObj3D::torusX(glm::f32vec3 pos, glm::f32vec3 center, float R, float r)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rev   = convertToSTR(R);
    float rad   = convertToSTR(r);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return rad*rad - shX*shX - shY*shY - shZ*shZ - rev*rev + 2.0f*rev*std::sqrt(shZ*shZ + shY*shY);
}

float FRepObj3D::torusY(glm::f32vec3 pos, glm::f32vec3 center, float R, float r)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rev   = convertToSTR(R);
    float rad   = convertToSTR(r);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return rad*rad - shX*shX - shY*shY - shZ*shZ - rev*rev + 2.0f*rev*std::sqrt(shZ*shZ + shX*shX);
}

float FRepObj3D::torusZ(glm::f32vec3 pos, glm::f32vec3 center, float R, float r)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);
    float rev   = convertToSTR(R);
    float rad   = convertToSTR(r);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return rad*rad - shX*shX - shY*shY - shZ*shZ - rev*rev + 2.0f*rev*std::sqrt(shX*shX + shY*shY);
}

float FRepObj3D::box(glm::f32vec3 pos, glm::f32vec3 center, float w, float h, float d)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0   = convertToSTR(center);
    float width  = convertToSTR(w);
    float height = convertToSTR(h);
    float depth  = convertToSTR(d);

    float planeXYfr  = c0.z + depth / 2.0f;
    float planeXYbc  = c0.z - depth / 2.0f;
    float planeXZbt  = c0.x + height / 2.0f;
    float planeXZtop = c0.x - height / 2.0f;
    float planeYZr   = c0.y + width / 2.0f;
    float planeYZl   = c0.y - width / 2.0f;

    float halfBox1 = intersect_function(planeXYfr, intersect_function(planeYZr, planeXZbt, 0.0, 1.0), 0.0f, 1.0f);
    float halfBox2 = intersect_function(planeXYbc, intersect_function(planeYZl, planeXZtop, 0.0, 1.0), 0.0f, 1.0f);
    return intersect_function(halfBox1, halfBox2, 0.0f, 1.0f);
}

float FRepObj3D::ellipsoid(glm::f32vec3 pos, glm::f32vec3 center, float a, float b, float c)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return 1.0f/resolutionX - shX*shX/(a*a) - shY*shY/(b*b) - shZ*shZ/(c*c);
}

float FRepObj3D::superEllipsoid(glm::f32vec3 pos, glm::f32vec3 center, float a, float b, float c, float s1, float s2)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    float val0 = shX*shX/(s2*a*a) + shY*shY/(s2*b*b);
    return 1.0f/resolutionX - std::pow(val0, s2/s1) - std::pow(shZ/c, 2.0f/s1);
}

float FRepObj3D::heart3D(glm::f32vec3 pos, glm::f32vec3 center)
{
    glm::f32vec3 str = convertToSTR(pos);
    glm::f32vec3 c0  = convertToSTR(center);

    float shX = str.x - c0.x;
    float shY = str.y - c0.y;
    float shZ = str.z - c0.z;

    return - std::pow((2.0f*scale*scale*shZ*shZ + scale*scale*shY*shY + scale*scale*shX*shX - 1.0f), 3.0f) +
             scale*scale*shY*shY * scale*scale*scale*shZ*shZ*shZ/20.0f +
            scale*scale*shX*shX * scale*scale*scale*shY*shY*shY;
}

float FRepObj3D::blending_union(float f1, float f2, float a0, float a1, float a2)
{
    float f1f2  = union_function(f1, f2, 0.0f, 1.0f);
    float blFun = a0 / (1.0f + (f1/a1)*(f1/a1) + (f2/a2)*(f1/a1));
    return f1f2 + blFun;
}

float FRepObj3D::blending_intersection(float f1, float f2, float a0, float a1, float a2)
{
    float f1f2  = intersect_function(f1, f2, 0.0f, 1.0f);
    float blFun = a0 / (1.0f + (f1/a1)*(f1/a1) + (f2/a2)*(f1/a1));
    return f1f2 + blFun;
}

float FRepObj3D::blending_subtraction(float f1, float f2, float a0, float a1, float a2)
{
    float f1f2  = subtract_function(f1, f2, 0.0f, 1.0f);
    float blFun = a0 / (1.0f + (f1/a1)*(f1/a1) + (f2/a2)*(f1/a1));
    return f1f2 + blFun;
}

glm::vec2 FRepObj3D::findZeroLevelSetInterval(std::vector<float> field, int numElemToAverage)
{
    std::vector<float> posVals, negVals;

    for(size_t i = 0; i < field.size(); i++)
    {
        if(field[i] >= 0)
            posVals.push_back(field[i]);
    }

    std::vector<float> minPosV;
    for(int j = 0; j < numElemToAverage; j++)
    {
        std::vector<float>::iterator minPos = std::min_element(posVals.begin(), posVals.end());
        minPosV.push_back(*minPos);
        posVals.erase(std::remove(posVals.begin(), posVals.end(), *minPos), posVals.end());
    }

    float minPos, minPsum = 0;
    for(size_t i = 0; i < minPosV.size(); i++)
        minPsum += minPosV[i];

    minPos = minPsum / static_cast<float>(minPosV.size());

    glm::vec2 result = glm::vec2(0.0f, minPos);
    return result;
}

std::vector<float> FRepObj3D::getFRep3D(glm::f32vec3 cent, std::function<float (glm::f32vec3, glm::f32vec3)> fun)
{
    frep.clear();
    for (int z = 0; z < resolutionZ; z++)
        for (int y = 0; y < resolutionY; y++)
            for (int x = 0; x < resolutionX; x++)
                frep.push_back(fun(glm::f32vec3(x, y, z), cent));

    return frep;
}

std::vector<float> FRepObj3D::getFRep3D(glm::f32vec3 cent, float a, std::function<float (glm::f32vec3, glm::f32vec3, float)> fun)
{
    frep.clear();
    for (int z = 0; z < resolutionZ; z++)
        for (int y = 0; y < resolutionY; y++)
            for (int x = 0; x < resolutionX; x++)
                frep.push_back(fun(glm::f32vec3(x, y, z), cent, a));

    return frep;
}

std::vector<float> FRepObj3D::getFRep3D(glm::f32vec3 cent, float a, float b, std::function<float (glm::f32vec3, glm::f32vec3, float, float)> fun)
{
    frep.clear();
    for (int z = 0; z < resolutionZ; z++)
        for (int y = 0; y < resolutionY; y++)
            for (int x = 0; x < resolutionX; x++)
                frep.push_back(fun(glm::f32vec3(x, y, z), cent, a, b));
    return frep;
}

std::vector<float> FRepObj3D::getFRep3D(glm::f32vec3 cent, float a, float b, float c, std::function<float(glm::f32vec3, glm::f32vec3, float, float, float)> fun)
{
    frep.clear();
    for (int z = 0; z < resolutionZ; z++)
        for (int y = 0; y < resolutionY; y++)
            for (int x = 0; x < resolutionX; x++)
                frep.push_back(fun(glm::f32vec3(x, y, z), cent, a, b, c));
    return frep;
}

std::vector<float> FRepObj3D::getFRep3D(glm::f32vec3 cent, float a, float b, float c, float d, std::function<float (glm::f32vec3, glm::f32vec3, float, float, float, float)> fun)
{
    frep.clear();
    for (int z = 0; z < resolutionZ; z++)
        for (int y = 0; y < resolutionY; y++)
            for (int x = 0; x < resolutionX; x++)
                frep.push_back(fun(glm::f32vec3(x, y, z), cent, a, b, c, d));
    return frep;

}

std::vector<float> FRepObj3D::getFRep3D(glm::f32vec3 cent, float a, float b, float c, float d, float e, std::function<float (glm::f32vec3, glm::f32vec3, float, float, float, float, float)> fun)
{
    frep.clear();
    for (int z = 0; z < resolutionZ; z++)
        for (int y = 0; y < resolutionY; y++)
            for (int x = 0; x < resolutionX; x++)
                frep.push_back(fun(glm::f32vec3(x, y, z), cent, a, b, c, d, e));
    return frep;
}

std::vector<float> FRepObj3D::getFRep3D(std::vector<float> f1, std::vector<float> f2, float alpha, float m, std::function<float (float, float, float, float)> fun)
{
    assert(f1.size() == f2.size());
    frep.clear();
    for (size_t i = 0; i < f1.size(); i++)
        frep.push_back(fun(f1[i], f2[i], alpha, m));
    return frep;
}

} //namespace frep3D_object

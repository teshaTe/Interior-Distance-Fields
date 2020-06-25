#ifndef H_FREP_3D_CLASS
#define H_FREP_3D_CLASS

#include <cmath>
#include <vector>
#include <functional>

#include <glm/glm.hpp>

namespace frep {

class FRepObj3D
{
//primitives
public:
    FRepObj3D(int resX, int resY, int resZ, float scaleF);
    ~FRepObj3D(){}

    float sphere ( glm::f32vec3 pos, glm::f32vec3 center, float R );
    float halfSphere( glm::f32vec3 pos, glm::f32vec3 center, float R );
    float blobby ( glm::f32vec3 pos, glm::f32vec3 center, float R );

    float cylinderX( glm::f32vec3 pos, glm::f32vec3 center, float R, float h );
    float cylinderY( glm::f32vec3 pos, glm::f32vec3 center, float R, float h );
    float cylinderZ( glm::f32vec3 pos, glm::f32vec3 center, float R, float h );

    float ellipticCylinderX( glm::f32vec3 pos, glm::f32vec3 center, float a, float b );
    float ellipticCylinderY( glm::f32vec3 pos, glm::f32vec3 center, float a, float b );
    float ellipticCylinderZ( glm::f32vec3 pos, glm::f32vec3 center, float a, float b );

    float coneX ( glm::f32vec3 pos, glm::f32vec3 center, float R, float h );
    float coneY ( glm::f32vec3 pos, glm::f32vec3 center, float R );
    float coneZ ( glm::f32vec3 pos, glm::f32vec3 center, float R );

    float ellipticConeX(glm::f32vec3 pos, glm::f32vec3 center, float a, float b );
    float ellipticConeY( glm::f32vec3 pos, glm::f32vec3 center, float h, float a, float b );
    float ellipticConeZ( glm::f32vec3 pos, glm::f32vec3 center, float h, float a, float b );

    float torusX ( glm::f32vec3 pos, glm::f32vec3 center, float R, float r );
    float torusY ( glm::f32vec3 pos, glm::f32vec3 center, float R, float r );
    float torusZ ( glm::f32vec3 pos, glm::f32vec3 center, float R, float r );

    float box       (glm::f32vec3 pos, glm::f32vec3 center, float w, float h, float d );
    float ellipsoid ( glm::f32vec3 pos, glm::f32vec3 center, float a, float b, float c );
    float plane( glm::f32vec3 pos, glm::f32vec3 center );

    float superEllipsoid( glm::f32vec3 pos, glm::f32vec3 center, float a, float b, float c, float s1, float s2 );
    float heart3D ( glm::f32vec3 pos, glm::f32vec3 center );

    float blending_union(float f1, float f2, float a0, float a1, float a2);
    float blending_intersection(float f1, float f2, float a0, float a1, float a2);
    float blending_subtraction(float f1, float f2, float a0, float a1, float a2);
    glm::vec2 findZeroLevelSetInterval(std::vector<float> field , int numElemToAverage = 20);

//operations over primitives and some other operations
public:
    std::vector<float> getFRep3D( glm::f32vec3 cent, std::function<float( glm::f32vec3, glm::f32vec3 )> fun);
    std::vector<float> getFRep3D( glm::f32vec3 cent, float a, std::function<float( glm::f32vec3, glm::f32vec3, float )> fun );
    std::vector<float> getFRep3D( glm::f32vec3 cent, float a, float b, std::function<float( glm::f32vec3, glm::f32vec3, float, float )> fun );
    std::vector<float> getFRep3D( glm::f32vec3 cent, float a, float b, float c, std::function<float(glm::f32vec3, glm::f32vec3, float, float, float )> fun );
    std::vector<float> getFRep3D( glm::f32vec3 cent, float a, float b, float c, float d, std::function<float( glm::f32vec3, glm::f32vec3, float, float, float, float )> fun );
    std::vector<float> getFRep3D( glm::f32vec3 cent, float a, float b, float c, float d, float e, std::function<float( glm::f32vec3, glm::f32vec3, float, float, float, float, float )> fun );
    std::vector<float> getFRep3D(std::vector<float> f1, std::vector<float> f2, float alpha, float m,
                                  std::function<float(float, float, float, float)> fun );

    inline void setNewResolution( int resX, int resY, int resZ ) { resolutionX = resX;
                                                                   resolutionY = resY;
                                                                   resolutionZ = resZ; }
    inline void setNewScale( float scaleF) { scale = scaleF; }

    inline float intersect_function(float f1, float f2, float alpha = 0.0f, float m = 0.0f)
    { return (1.0f/(1.0f+alpha))*( f1 + f2 - std::sqrt(f1 * f1 + f2 * f2 - 2*alpha*f1*f2)*std::pow(f1*f1+f2*f2, m/2.0f)); }

    inline float union_function(float f1, float f2, float alpha = 0.0f, float m = 0.0f)
    { return (1.0f/(1.0f+alpha))*( f1 + f2 + std::sqrt(f1 * f1 + f2 * f2 - 2*alpha*f1*f2)*std::pow(f1*f1+f2*f2, m/2.0f)); }

    inline float subtract_function(float f1, float f2, float alpha = 0.0f, float m = 0.0f) { return intersect_function(f1, -f2, alpha, m); }

    std::vector<float> rotateObj( std::vector<float> obj, float angle );
    std::vector<float> scaleObj ( std::vector<float> obj, float scale );

private:
    int resolutionX, resolutionY, resolutionZ;
    float scale;

    std::vector<float> frep;

    inline float convertToSTR( float val ) { return val / resolutionX; }
    inline glm::f32vec3 convertToSTR( glm::f32vec3 val ) { return glm::f32vec3( val.x/resolutionX, val.y/resolutionY, val.z/resolutionZ ); }

};

} //namespace frep3D_object
#endif

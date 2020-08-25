#include "include/frep2D.h"
#include <iostream>
#include <cmath>

namespace frep {

FRepObj2D::FRepObj2D(int resX, int resY, float scaleF) : resolutionX(resX), resolutionY(resY), scale(scaleF) {}

float FRepObj2D::triangle(glm::vec2 pos, glm::vec2 cent, float a, float b, float c)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    a = convertToUV(a); // right side
    b = convertToUV(b); // left side
    c = convertToUV(c); // the base of the triangle

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);

    float lb1 = shY + a - ( shX + b);
    float rb1 = shY + a - (-shX + b);
    return intersect_function(intersect_function(lb1, rb1), -shY+c);
}

float FRepObj2D::triangle2(glm::vec2 pos, glm::vec2 cent, const float a, const float b)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    glm::vec2 sides = convertToUV(glm::vec2(a, b));

    float rec = rectangle(pos, cent, a, b);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    float line = shY - (sides.y/sides.x)*shX;

    return intersect_function(rec, line);
}

float FRepObj2D::circle(glm::vec2 pos, glm::vec2 cent, float R)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    R = convertToUV(R);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return R*R - shX*shX - shY*shY;
}

float FRepObj2D::blobby2D(glm::vec2 pos, glm::vec2 cent, float R)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    R = convertToUV(R);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return R*R/(shX*shX + shY*shY);
}

float FRepObj2D::ellipticCylZ2D(glm::vec2 pos, glm::vec2 cent, float a, float b)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return 1.0f/resolutionX - std::pow(shX, 2.0f)/(a*a) - std::pow(shY, 2.0f)/(b*b);
}

float FRepObj2D::ellipsoid2D(glm::vec2 pos, glm::vec2 cent, float a, float b)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return 1.0f/resolutionX - std::pow(shX/a, 2.0f) - std::pow(shY/b, 2.0f);
}

float FRepObj2D::torusY2D(glm::vec2 pos, glm::vec2 cent, float R, float rev)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    R = convertToUV(R);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return R*R - std::pow(shY, 2.0f) - std::pow(shY, 2.0f) - rev*rev +
           2.0f*rev*std::sqrt(std::pow(shX, 2.0f) + std::pow(shY, 2.0f));
}

float FRepObj2D::torusZ2D(glm::vec2 pos, glm::vec2 cent, float R, float rev)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    R = convertToUV(R);

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);

    return R*R - std::pow(shX, 2.0f) - std::pow(shY, 2.0f) - rev*rev +
           2.0f*rev*std::sqrt(std::pow(shX, 2.0f) + std::pow(shY, 2.0f));
}

float FRepObj2D::rectangle(glm::vec2 pos, glm::vec2 cent, float w, float h)
{
    glm::vec2 uv  = convertToUV(pos);
    glm::vec2 c0  = convertToUV(cent);
    glm::vec2 rec = convertToUV(glm::vec2(w, h));

    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return intersect_function(rec.x - std::abs(shX), rec.y - std::abs(shY));
}

float FRepObj2D::heart2D(glm::vec2 pos, glm::vec2 cent)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return - std::pow(shX*shX + shY*shY - 1.0f, 3.0f) + shX*shX*std::pow(shY, 3.0f);
}

float FRepObj2D::decocube2D(glm::vec2 pos, glm::vec2 cent)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);
    return -(std::pow(shX*shX + shY*shY - 2.0f, 2.0f) + 1.0f) *
            (std::pow(shY*shY - 2.0f, 2.0f) + std::pow(shX*shX-1.0f, 2.0f)) *
            (std::pow(shX*shX - 2.0f, 2.0f) + std::pow(shY*shY-1.0f, 2.0f)) + 50.0f/resolutionX;
}

float FRepObj2D::suriken(glm::vec2 pos, glm::vec2 cent)
{
    glm::vec2 uv = convertToUV(pos);
    glm::vec2 c0 = convertToUV(cent);
    float shX = scale*(uv.x - c0.x);
    float shY = scale*(uv.y - c0.y);

    float lb1 = shY + scaleToNewRange(358.4)/resolutionY - ( shX*2.0f + scaleToNewRange(153.6)/resolutionX);
    float rb1 = shY + scaleToNewRange(358.4)/resolutionY - (-shX*2.0f + scaleToNewRange(153.6)/resolutionX);
    float lb2 = shY - scaleToNewRange(51.2)/resolutionY  - ( shX*2.0f + scaleToNewRange(153.6)/resolutionY);
    float rb2 = shY - scaleToNewRange(51.2)/resolutionY  - (-shX*2.0f + scaleToNewRange(153.6)/resolutionY);
    float lb3 = shX + scaleToNewRange(358.4)/resolutionX - ( shY*2.0f + scaleToNewRange(153.6)/resolutionY);
    float rb3 = shX + scaleToNewRange(358.4)/resolutionX - (-shY*2.0f + scaleToNewRange(153.6)/resolutionY);
    float lb4 = shX - scaleToNewRange(51.2)/resolutionX  - ( shY*2.0f + scaleToNewRange(153.6)/resolutionY);
    float rb4 = shX - scaleToNewRange(51.2)/resolutionX  - (-shY*2.0f + scaleToNewRange(153.6)/resolutionY);

    float trian1 = intersect_function(intersect_function( lb1,  rb1), -shY+scaleToNewRange(102.4)/resolutionY);
    float trian2 = intersect_function(intersect_function(-lb2, -rb2),  shY+scaleToNewRange(102.4)/resolutionY);
    float trian3 = intersect_function(intersect_function( lb3,  rb3), -shX+scaleToNewRange(102.4)/resolutionX);
    float trian4 = intersect_function(intersect_function(-lb4, -rb4),  shX+scaleToNewRange(102.4)/resolutionX);

    return union_function(union_function(union_function(trian1, trian2), trian3), trian4);
}

float FRepObj2D::bat(glm::vec2 pos, glm::vec2 center)
{
    //glm::vec2 center = glm::vec2(resolutionX/2.0f, resolutionY/2.0f);

    float rec0 = rectangle(pos, center, scaleToNewRange(50), scaleToNewRange(30)); //main body
    float rec1 = rectangle(pos, glm::vec2(center.x + scaleToNewRange(80)/scale, center.y - scaleToNewRange(5.5)/scale),
                                                                 scaleToNewRange(40), scaleToNewRange(36));
    float rec2 = rectangle(pos, glm::vec2(center.x - scaleToNewRange(80)/scale, center.y - scaleToNewRange(5.5)/scale),
                                                                 scaleToNewRange(40), scaleToNewRange(36));
    float rec3 = rectangle(pos, glm::vec2(center.x + scaleToNewRange(85)/scale, center.y - scaleToNewRange(45)/scale),
                                                                 scaleToNewRange(35), scaleToNewRange(9)); //right part of the wing
    float rec4 = rectangle(pos, glm::vec2(center.x - scaleToNewRange(85)/scale, center.y - scaleToNewRange(45)/scale),
                                                                 scaleToNewRange(35), scaleToNewRange(9)); //left part of the wing
    float rec5 = rectangle(pos, glm::vec2(center.x, center.y - scaleToNewRange(35)/scale),
                                           scaleToNewRange(12), scaleToNewRange(6));  //head

    float tri0 = triangle2(getRotatedCoords(pos, 45), getRotatedCoords(glm::vec2(center.x, center.y + scaleToNewRange(20)/scale), 45),
                                                     scaleToNewRange(40), scaleToNewRange(40)); //bat tale, bottom triangle
    float tri1 = triangle2(pos, glm::vec2(center.x - scaleToNewRange(34)/scale, center.y - scaleToNewRange(33)/scale),
                                                                  scaleToNewRange(7), scaleToNewRange(7));                                              //supportive triangle (bot) for left wing, internal
    float tri2 = triangle2(getRotatedCoords(pos, 90), getRotatedCoords(glm::vec2(center.x + scaleToNewRange(34)/scale, center.y - scaleToNewRange(34)/scale), 90),
                           scaleToNewRange(7), scaleToNewRange(7));   //supportive triangle (bot) for right wing, internal
    float tri3 = triangle2(pos, glm::vec2(center.x - scaleToNewRange(44)/scale, center.y - scaleToNewRange(45)/scale),
                           scaleToNewRange(7), scaleToNewRange(9));   //supportive triangle (up) for left wing, internal
    float tri4 = triangle2(getRotatedCoords(pos, 90), getRotatedCoords(glm::vec2(center.x + scaleToNewRange(44)/scale, center.y - scaleToNewRange(45)/scale), 90),
                           scaleToNewRange(9), scaleToNewRange(7));   //supportive triangle (up) for right wing, internal
    float tri5 = triangle2(pos, glm::vec2(center.x + scaleToNewRange(8)/scale, center.y - scaleToNewRange(50)/scale),
                           scaleToNewRange(4), scaleToNewRange(10));  //right ear
    float tri6 = triangle2(getRotatedCoords(pos, 90), getRotatedCoords(glm::vec2(center.x - scaleToNewRange(8)/scale, center.y - scaleToNewRange(50)/scale), 90),
                           scaleToNewRange(10), scaleToNewRange(4));  //left ear
    float tri7 = triangle2(getRotatedCoords(pos, 180), getRotatedCoords(glm::vec2(center.x - scaleToNewRange(152)/scale, center.y - scaleToNewRange(19.5)/scale), 180),
                           scaleToNewRange(41), scaleToNewRange(35)); //left wing
    float tri8 = triangle2(getRotatedCoords(pos, 270), getRotatedCoords(glm::vec2(center.x + scaleToNewRange(150)/scale, center.x - scaleToNewRange(19.5)/scale), 270),
                           scaleToNewRange(35), scaleToNewRange(41)); //right wing

    float union1  = union_function(rec0, rec1);
    float union2  = union_function(rec2, union1);
    float union3  = union_function(rec3, union2);
    float union4  = union_function(rec4, union3);
    float union5  = union_function(rec5, union4);
    float union6  = union_function(tri0, union5);
    float union7  = union_function(tri1, union6);
    float union8  = union_function(tri2, union7);
    float union9  = union_function(tri3, union8);
    float union10 = union_function(tri4, union9);
    float union11 = union_function(tri5, union10);
    float union12 = union_function(tri6, union11);
    float union13 = union_function(tri7, union12);
    float bat     = union_function(tri8, union13);

    return bat;
}

float FRepObj2D::trebleClef(glm::vec2 pos, glm::vec2 center)
{
    //glm::vec2 center = glm::vec2(resolutionX/2.0f, resolutionY/2.0f);

    //constructing bottom part
    float rec_base = rectangle(pos, glm::vec2(center.x, center.y - scaleToNewRange(20)), scaleToNewRange(15), scaleToNewRange(120));
    float ark1 = circle(pos, glm::vec2(center.x - scaleToNewRange(20)/scale, center.y + scaleToNewRange(95)/scale), scaleToNewRange(35));
    float ark2 = circle(pos, glm::vec2(center.x - scaleToNewRange(50)/scale, center.y + scaleToNewRange(75)/scale), scaleToNewRange(40));
    float ark3 = circle(pos, glm::vec2(center.x - scaleToNewRange(30)/scale, center.y + scaleToNewRange(95)/scale), scaleToNewRange(35));
    float ark4 = circle(pos, glm::vec2(center.x - scaleToNewRange(20)/scale, center.y + scaleToNewRange(70)/scale), scaleToNewRange(35));
    float circ_bot = circle(pos, glm::vec2(center.x - scaleToNewRange(40)/scale, center.y + scaleToNewRange(85)/scale), scaleToNewRange(20));

    float arkB1 = subtract_function(ark1, ark2);
    float arkB2 = subtract_function(ark3, ark4);

    float union1  = union_function(arkB1, rec_base);
    float union2  = union_function(arkB2, union1);
    float fin_bot = union_function(circ_bot, union2);

    //constructing middle part
    float ark5 = circle(pos, glm::vec2(center.x + scaleToNewRange(12)/scale, center.y - scaleToNewRange(10)/scale), scaleToNewRange(60));
    float ark6 = circle(pos, glm::vec2(center.x - scaleToNewRange(5)/scale, center.y - scaleToNewRange(5)/scale), scaleToNewRange(49));
    float ark7 = circle(pos, glm::vec2(center.x, center.y - scaleToNewRange(28)/scale), scaleToNewRange(80)/scale);
    float ark8 = circle(pos, glm::vec2(center.x + scaleToNewRange(15)/scale, center.y - scaleToNewRange(28)/scale), scaleToNewRange(70));

    float arkM1 = subtract_function(ark5, ark6);
    float arkM2 = subtract_function(ark7, ark8);

    float union3     = union_function(arkM1, fin_bot);
    float fin_middle = union_function(arkM2, union3);

    //constructing upper part
    float ark9  = circle(pos, glm::vec2(center.x, center.y - scaleToNewRange(155)/scale), scaleToNewRange(57));
    float ark10 = circle(pos, glm::vec2(center.x - scaleToNewRange(30)/scale, center.y - scaleToNewRange(169)/scale), scaleToNewRange(65));
    float ark11 = circle(pos, glm::vec2(center.x + scaleToNewRange(44)/scale, center.y - scaleToNewRange(150)/scale), scaleToNewRange(60));
    float ark12 = circle(pos, glm::vec2(center.x + scaleToNewRange(87)/scale, center.y - scaleToNewRange(160)/scale), scaleToNewRange(75));
    float rec_helper1 = rectangle(pos, glm::vec2(center.x, center.y - scaleToNewRange(87)/scale), scaleToNewRange(100), scaleToNewRange(50));
    float circ_top    = circle(pos, glm::vec2(center.x + scaleToNewRange(27)/scale, center.y - scaleToNewRange(201)/scale), scaleToNewRange(5));

    float arkU1  = subtract_function(ark9, ark10);
    float arkU2  = subtract_function(ark11, ark12);
    float arkU2f = subtract_function(arkU2, rec_helper1);

    float union4 = union_function(arkU1, fin_middle);
    float union5 = union_function(circ_top, union4);
    float trebleclef = union_function(arkU2f, union5);

    return trebleclef;
}

float FRepObj2D::bounded_blending(float f1, float f2, float a0, float a1, float a2, float a3,
                                  float time, float alpha, float m)
{
      f1 = intersect_function(f1, -time, alpha, m);
      f2 = intersect_function(f2, time - 1.0f, alpha, m);
      float f3 = intersect_function(time + 10.0f, 10.0f - time, alpha, m);

      float r1 = (f1/a1)*(f1/a1) + (f2/a2)*(f2/a2);
      float r2 = 0.0f;

      if(f3 > 0.0f)
        r2 = (f3/a3) * (f3/a3);

      float rr = 0.0f;
      if(r1 > 0.0f)
        rr = r1 / (r1 + r2);

      float d = 0.0f;
      if(rr < 1.0f)
        d = a0 * (1.0f - rr)*(1.0f - rr)*(1.0f - rr) / (1.0f + rr);

      return union_function(f1, f2, alpha, m) + d;
}


float FRepObj2D::union_function_R0(float f1, float f2, float n)
{
    float result = 0.0f;

    if(f1 > 0.0f && f2 > 0.0f)
       result = std::pow(std::pow(f1, n) + std::pow(f2, n), 1.0f/n);
    else if(f1 <= 0.0f && f2 >= 0.0f)
       result = f2;
    else if(f1 >= 0.0f && f2 <= 0.0f)
       result = f1;
    else if(f1 < 0.0f && f2 < 0.0f)
       result = std::pow(-1.0f, n+1.0f) * f1*f2 * std::pow(std::pow(f1, n) + std::pow(f2, n), -1.0f/n);

    return result;
}

glm::vec2 FRepObj2D::findZeroLevelSetInterval(std::vector<float> field, int numElemToAverage)
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

std::vector<float> FRepObj2D::getFRep2D(std::function<float(glm::vec2)> fun)
{
    frep.clear();
    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(glm::vec2(x, y)));

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, std::function<float(glm::vec2, glm::vec2)> fun)
{
    frep.clear();
    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(glm::vec2(x, y), cent));

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, float R,
                                        std::function<float(glm::vec2, glm::vec2, float)> fun)
{
    frep.clear();
    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(glm::vec2(x, y), cent, R));

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, float p1, float p2,
                                        std::function<float(glm::vec2, glm::vec2, float, float)> fun)
{
    frep.clear();
    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(glm::vec2(x, y), cent, p1, p2));

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, float p1, float p2, float p3,
                                        std::function<float(glm::vec2, glm::vec2, float, float, float)> fun)
{
    frep.clear();
    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(glm::vec2(x, y), cent, p1, p2, p3));

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(std::vector<float> f1, std::vector<float> f2, float alpha, float m,
                                         std::function<float (float, float, float, float)> fun)
{
    frep.clear();
    for(int i = 0; i < resolutionX*resolutionY; i++)
        frep.push_back(fun(f1[i], f2[i], alpha, m));

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(std::vector<float> f1, std::vector<float> f2, float n, std::function<float (float, float, float)> fun)
{
    frep.clear();
    for(int i = 0; i < resolutionX*resolutionY; i++)
        frep.push_back(fun(f1[i], f2[i], n));

    return frep;
}

glm::vec2 FRepObj2D::getRotatedCoords(glm::vec2 inCoords, const float angle)
{
    float angle0 = static_cast<float>((angle/180.0f) * M_PI);

    float rotX = inCoords.x*cosf(angle0) - inCoords.y*sinf(angle0);
    float rotY = inCoords.x*sinf(angle0) + inCoords.y*cosf(angle0);
    return glm::vec2(rotX, rotY);
}

std::vector<float> FRepObj2D::getRotatedFrep2D(glm::vec2 cent, float w, float h,
                                                float angle, std::function<float (glm::vec2, glm::vec2, float, float)> fun)
{
    frep.clear();
    glm::vec2 c0 = getRotatedCoords(cent, angle);

    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(getRotatedCoords(glm::vec2(x,y), angle), c0, w, h));

    return frep;
}

std::vector<float> FRepObj2D::getRotatedFrep2D(glm::vec2 cent, float a, float b, float c, float angle,
                                                std::function<float (glm::vec2, glm::vec2, float, float, float, float)> fun)
{
    frep.clear();
    glm::vec2 c0 = getRotatedCoords(cent, angle);

    for(int y = 0; y < resolutionY; y++)
        for(int x = 0; x < resolutionX; x++)
            frep.push_back(fun(getRotatedCoords(glm::vec2(x, y), angle), c0, a, b, c, angle));

    return frep;
}

} //namespace frep

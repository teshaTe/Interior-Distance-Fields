#include "include/frep2D.h"
#include <iostream>
#include <cmath>

namespace frep {

FRepObj2D::FRepObj2D(int resX, int resY, float scaleF )
{
    resolutionX = resX;
    resolutionY = resY;
    scale       = scaleF;
}

float FRepObj2D::triangle(glm::vec2 pos, glm::vec2 cent, float a, float b, float c)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    a = convertToUV( a ); // right side
    b = convertToUV( b ); // left side
    c = convertToUV( c ); // the base of the triangle

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;

    float lb1 = shY + a - ( shX + b);
    float rb1 = shY + a - (-shX + b);
    return intersect_function( intersect_function( lb1,  rb1, 0.0f, 0.0f), -shY+c, 0.0f, 0.0f );
}

float FRepObj2D::triangle2(glm::vec2 pos, glm::vec2 cent, const float a, const float b)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    glm::vec2 sides = convertToUV( glm::vec2(a, b) );

    float rec = rectangle( pos, cent, a, b );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    float line = shY - (sides.y/sides.x)*shX;

    return intersect_function( rec, line, 0.0f, 0.0f );
}

float FRepObj2D::circle(glm::vec2 pos, glm::vec2 cent, float R)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    R = convertToUV( R );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return R*R - shX*shX - shY*shY;
}

float FRepObj2D::blobby2D(glm::vec2 pos, glm::vec2 cent, float R)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    R = convertToUV( R );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return R*R/( shX*shX + shY*shY);
}

float FRepObj2D::ellipticCylZ2D(glm::vec2 pos, glm::vec2 cent, float a, float b)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return 1.0f/resolutionX - std::pow( shX, 2.0f )/( a*a ) - std::pow( shY, 2.0f )/( b*b );
}

float FRepObj2D::ellipsoid2D(glm::vec2 pos, glm::vec2 cent, float a, float b)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return 1.0f/resolutionX - std::pow( ( shX )/a, 2.0f ) - std::pow( ( shY )/b, 2.0f );
}

float FRepObj2D::torusY2D(glm::vec2 pos, glm::vec2 cent, float R, float rev)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    R = convertToUV( R );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return R*R - std::pow( shY, 2.0f ) - std::pow( shY, 2.0f ) - rev*rev +
            2.0f*rev*std::sqrt( std::pow( shX, 2.0f ) + std::pow( shY, 2.0f ));
}

float FRepObj2D::torusZ2D(glm::vec2 pos, glm::vec2 cent, float R, float rev)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    R = convertToUV( R );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;

    return R*R - std::pow( shX, 2.0f ) - std::pow( shY, 2.0f ) -
            rev*rev + 2.0f*rev*std::sqrt( std::pow( shX, 2.0f ) + std::pow( shY, 2.0f ) );
}

float FRepObj2D::rectangle(glm::vec2 pos, glm::vec2 cent, float w, float h )
{
    glm::vec2 uv  = convertToUV( pos );
    glm::vec2 c0  = convertToUV( cent );
    glm::vec2 rec = convertToUV( glm::vec2(w, h) );

    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return intersect_function( rec.x - std::abs(shX), rec.y - std::abs(shY), 0.0f, 0.0f );
}

float FRepObj2D::heart2D(glm::vec2 pos, glm::vec2 cent)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;
    return - std::pow( scale*scale*shX*shX + scale*scale*shY*shY - 1.0f, 3.0f ) +
            scale*scale*shX*shX*std::pow( scale*shY, 3.0f );
}

float FRepObj2D::decocube2D(glm::vec2 pos, glm::vec2 cent)
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;

    return -( std::pow( scale*scale*shX*shX + scale*scale*shY*shY - 2.0f, 2.0f) + 1.0f ) *
           ( std::pow( scale*scale*shY*shY - 2.0f, 2.0f ) + std::pow(scale*scale*shX*shX-1.0f, 2.0f)) *
           ( std::pow( scale*scale*shX*shX - 2.0f, 2.0f ) + std::pow(scale*scale*shY*shY-1.0f, 2.0f)) + 50.0f;
}

float FRepObj2D::suriken( glm::vec2 pos, glm::vec2 cent )
{
    glm::vec2 uv = convertToUV( pos );
    glm::vec2 c0 = convertToUV( cent );
    float shX = uv.x - c0.x;
    float shY = uv.y - c0.y;

    float lb1 = shY + 0.7f - ( shX*2.0f + 0.3f);
    float rb1 = shY + 0.7f - (-shX*2.0f + 0.3f);
    float lb2 = shY - 0.1f - ( shX*2.0f + 0.3f);
    float rb2 = shY - 0.1f - (-shX*2.0f + 0.3f);
    float lb3 = shX + 0.7f - ( shY*2.0f + 0.3f);
    float rb3 = shX + 0.7f - (-shY*2.0f + 0.3f);
    float lb4 = shX - 0.1f - ( shY*2.0f + 0.3f);
    float rb4 = shX - 0.1f - (-shY*2.0f + 0.3f);

    float trian1 = intersect_function( intersect_function( lb1,  rb1, 0.0f, 0.0f), -shY+0.2f, 0.0f, 0.0f );
    float trian2 = intersect_function( intersect_function(-lb2, -rb2, 0.0f, 0.0f),  shY+0.2f, 0.0f, 0.0f );
    float trian3 = intersect_function( intersect_function( lb3,  rb3, 0.0f, 0.0f), -shX+0.2f, 0.0f, 0.0f );
    float trian4 = intersect_function( intersect_function(-lb4, -rb4, 0.0f, 0.0f),  shX+0.2f, 0.0f, 0.0f );

    return union_function( union_function( union_function( trian1, trian2, 0.0f, 0.0f ),
                                             trian3, 0.0f, 0.0f ), trian4, 0.0f, 0.0f );
}

float FRepObj2D::elf( glm::vec2 pos )
{
    glm::vec2 uv = convertToUV( pos );
    float shX  = resolutionX/2.0f,
          shY  = resolutionY/3.0f;
    //head
    glm::vec2 headC = glm::vec2( (0.0f+shX), (4.5f+shY) );
    float head    = ellipsoid2D( pos, headC, 6.0f, 5.5f );

    //nose
    glm::vec2 noseC = glm::vec2( 0.0f+shX, 25.5f+shY );
    float nose    = ellipsoid2D( pos, noseC, 0.5f, 0.5f );

    //face
    glm::vec2 faceC1 = glm::vec2( 0.0f+shX, 4.5f+shY );
    float xFace1   = ellipsoid2D( pos, faceC1, 6.0f, 5.5f );
    glm::vec2 faceC2 = glm::vec2( 0.0f+shX, 3.7f+shY );
    float xFace2   = ellipticCylZ2D( pos, faceC2, 4.7f, 4.5f );
    float face     = intersect_function( xFace1, xFace2, 0.0f, 0.0f );

    //mouth
    glm::vec2 mouthC = glm::vec2( 0.0f+shX, 10.5f+shY );
    float mouth1   = ellipsoid2D( pos, mouthC, 3.0f, 2.3f );
    float mouth2   = -(60.5f+shY)/resolutionY + uv.y;
    float mouth    = intersect_function( mouth1, mouth2, 0.0f, 0.0f );

    //eyes
    glm::vec2 eye1C = glm::vec2( 0.0f+shX, 2.5f+shY );
    float eye1    = ellipsoid2D( pos, eye1C, 6.0f, 5.5f );

    glm::vec2 eye2C = glm::vec2( 36.5f+shX, -10.5f+shY );
    float eye2    = ellipticCylZ2D( pos, eye2C, 1.5f, 1.4f );
    float eye3    = intersect_function( eye1, eye2, 0.0f, 0.0f );

    glm::vec2 eye4C = glm::vec2( -36.5f+shX, -10.5f+shY );
    float eye4    = ellipticCylZ2D( pos, eye4C, 1.5f, 1.4f );
    float eye5    = intersect_function( eye1, eye4, 0.0f, 0.0f);

    glm::vec2 eye6C = glm::vec2( 20.3f+shX, -10.5f+shY );
    float eye6    = ellipsoid2D( pos, eye6C, 0.3f, 0.3f );
    float eye7    = subtract_function( eye3, eye6, 0.0f, 0.0f );

    glm::vec2 eye8C = glm::vec2( -20.3f+shX, -10.5f+shY );
    float eye8    = ellipsoid2D( pos, eye8C, 0.3f, 0.3f );
    float eye9    = subtract_function( eye5, eye8, 0.0f, 0.0f );

    //head final
    float fin_head0 = union_function( eye9, eye7, 0.0f, 0.0f );
    float fin_head1 = union_function( face, head, 0.0f, 0.0f );
    float fin_head2 = subtract_function( fin_head1, mouth, 0.0f, 0.0f );
    float fin_head3 = subtract_function( fin_head2, nose, 0.0f, 0.0f );
    float fin_head  = subtract_function( fin_head3, fin_head0, 0.0f, 0.0f );

    //body
    glm::vec2 bodyC = glm::vec2( 0.0f+shX, 122.0f+shY );
    float body1   = ellipsoid2D( pos, bodyC, 4.5f, 6.0f );
    float body2   = intersect_function( intersect_function( body1, (126.0f+shY)/resolutionY+uv.y, 0.0f, 0.0f ),
                                                                  -(121.0f+shY)/resolutionY+uv.y, 0.0f, 0.0f );

    //neck
    glm::vec2 neckC = glm::vec2( 0.0f+shX, 120.0f+shY );
    float neck    = ellipsoid2D( pos, neckC, 1.8f, 0.8f );

    //bell
    glm::vec2 bellC = glm::vec2( 0.0f+shX, 180.0f+shY );
    float bell    = circle( pos, bellC, 10.8f);

    float body3 = union_function( body2, neck, 0.0f, 0.0f );

    //belly
    glm::vec2 bellyC1 = glm::vec2( 0.0f+shX, 32.0f+shY );
    float belly1    = ellipsoid2D( pos, bellyC1, 4.5, 6.0 );

    glm::vec2 bellyC2 = glm::vec2( 0.0f+shX, 200.0f+shY );
    float belly2    = ellipticCylZ2D( pos, bellyC2, 0.5f, 0.5f );
    float belly3    = intersect_function( belly1, belly2, 0.0f, 0.0f );

    glm::vec2 bellyC4 = glm::vec2( 0.0f+shX, 138.0f+shY );
    float belly4    = ellipticCylZ2D( pos, bellyC4, 2.0f, 1.6f );

    glm::vec2 bellyC5 = glm::vec2( 0.0f+shX, 200.0f+shY );
    float belly5    = ellipsoid2D( pos, bellyC5, 1.1f, 1.1f );
    float belly6    = -(103.0f+shY)/resolutionY + uv.y;

    float fin_body0 = intersect_function( belly4, belly6, 0.0f, 0.0f);
    float fin_body1 = union_function( fin_body0, belly5, 0.0f, 0.0f );
    float fin_body2 = subtract_function( body3, fin_body1, 0.0f, 0.0f );
    float fin_body3 = union_function( fin_body2, belly3, 0.0f, 0.0f );
    float fin_body  = union_function( fin_body3, bell, 0.0f, 0.0f );

    //legs and feet
    glm::vec2 l_legC = glm::vec2( 93.0f+shX, 247.5f+shY );
    float l_leg    = ellipsoid2D( pos, l_legC, 3.0f, 1.5f );
    glm::vec2 r_legC = glm::vec2( -93.0f+shX, 247.5f+shY );
    float r_leg    = ellipsoid2D( pos, r_legC, 3.0, 1.5 );

    glm::vec2 l_footC = glm::vec2( 92.0f+shX, 246.0f+shY );
    float l_foot    = ellipsoid2D( pos, l_footC, 1.0f, 1.0f );
    glm::vec2 r_footC = glm::vec2(-92.0f+shX, 246.0f+shY);
    float r_foot    = ellipsoid2D( pos, r_footC, 1.0f, 1.0f );

    float legs0 = union_function( l_leg, l_foot, 0.0f, 0.0f);
    float legs1 = union_function( r_leg, r_foot, 0.0f, 0.0f );
    float legs  = union_function( legs0, legs1, 0.0f, 0.0f );

    //right mustache
    float rh = 0.15f;
    glm::vec2 r_mustC1 = glm::vec2( 112.0f+shX, 100.0f+shY );
    float r_must1    = torusZ2D( pos, r_mustC1, 7.0f, rh );
    glm::vec2 r_mustC2 = glm::vec2( 122.0f+shX, 102.0f+shY );
    float r_must2    = torusZ2D( pos, r_mustC2, 7.0f, rh );
    glm::vec2 r_mustC3 = glm::vec2( 112.0f+shX, 104.0f+shY );
    float r_must3    = torusZ2D( pos, r_mustC3, 7.0f, rh );

    //left mustache
    glm::vec2 l_mustC1 = glm::vec2( -112.0f+shX, 100.0f+shY );
    float l_must1    = torusZ2D( pos, l_mustC1, 7.0f, rh );
    glm::vec2 l_mustC2 = glm::vec2( -112.0f+shX, 102.0f+shY );
    float l_must2    = torusZ2D( pos, l_mustC2, 7.0f, rh );
    glm::vec2 l_mustC3 = glm::vec2( -112.0f+shX, 104.0f+shY );
    float l_must3    = torusZ2D( pos, l_mustC3, 7.0f, rh );

    float r_must11 = intersect_function( uv.y, (1.4f*shY)/resolutionY - uv.y, 0.0f, 0.0f );
    float r_must22 = union_function( union_function( r_must1, r_must2, 0.0f, 0.0f ), r_must3, 0.0f, 0.0f );
    float r_mustF  = intersect_function( intersect_function( r_must22, r_must11, 0.0f, 0.0f ), uv.y, 0.0f, 0.0f );

    float l_must11 = intersect_function( -uv.y + (1.4f*shY)/resolutionY, uv.y, 0.0f, 0.0f );
    float l_must22 = union_function( union_function( l_must1, l_must2, 0.0f, 0.0f ), l_must3, 0.0f, 0.0f );
    float l_mustF  = intersect_function( intersect_function( l_must22, l_must11, 0.0f, 0.0f), uv.y, 0.0f, 0.0f );
    float mustache = union_function( r_mustF, l_mustF, 0.0f, 0.0f );

    //hands
    glm::vec2 r_handC = glm::vec2( 158.0f+shX, 152.5f+shY );
    float r_hand    = circle( pos, r_handC, 30.0f );
    glm::vec2 l_handC = glm::vec2( -158.0f+shX, 152.5f+shY );
    float l_hand    = circle( pos, l_handC, 30.0f );

    float ude2   = uv.y + (150.0f+shY)/resolutionY;
    float ude3   = (150.0f+shY)/resolutionY - uv.y;
    glm::vec2 udeC = glm::vec2( 0.0f+shX, 251.0f+shY );
    float ude1   = torusZ2D( pos, udeC, 30.0f, 0.3f );
    float te1    = intersect_function( intersect_function( ude1, ude2, 0.0f, 0.0f ), ude3, 0.0f, 0.0f );
    float te2    = union_function( union_function( te1, r_hand, 0.0f, 0.0f ), l_hand, 0.0f, 0.0f );

    //final model
    float elf0 = union_function( fin_head, fin_body, 0.0f, 0.0f );
    float elf1 = union_function( elf0, legs, 0.0f, 0.0f );
    float elf2 = union_function( elf1, te2, 0.0f, 0.0f );
    float elf  = union_function( elf2, mustache, 0.0f, 0.0f );

    return elf;
}

float FRepObj2D::bat(glm::vec2 pos)
{
    float rec0 = rectangle( pos, glm::vec2(250, 250), 50, 30 ); //main body
    float rec1 = rectangle( pos, glm::vec2(330, 244.5), 40, 36 );
    float rec2 = rectangle( pos, glm::vec2(170, 244.5), 40, 36 );
    float rec3 = rectangle( pos, glm::vec2(335, 205), 35, 9 ); //right part of the wing
    float rec4 = rectangle( pos, glm::vec2(165, 205), 35, 9 ); //left part of the wing
    float rec5 = rectangle( pos, glm::vec2(250, 215), 12, 6 );  //head

    float tri0 = triangle2( getRotatedCoords(pos, 45), getRotatedCoords(glm::vec2(250, 270),45), 40, 40 ); //bat tale, bottom triangle
    float tri1 = triangle2( pos, glm::vec2(216, 217), 7, 7 );                                              //supportive triangle (bot) for left wing, internal
    float tri2 = triangle2( getRotatedCoords(pos, 90), getRotatedCoords(glm::vec2(284, 216),90), 7, 7 );   //supportive triangle (bot) for right wing, internal
    float tri3 = triangle2( pos, glm::vec2(206, 205), 7, 9 );                                            //supportive triangle (up) for left wing, internal
    float tri4 = triangle2( getRotatedCoords(pos, 90), getRotatedCoords(glm::vec2(294, 205),90), 9, 7 ); //supportive triangle (up) for right wing, internal
    float tri5 = triangle2( pos, glm::vec2(258, 200), 4, 10 );                                             //right ear
    float tri6 = triangle2( getRotatedCoords(pos, 90), getRotatedCoords(glm::vec2(242, 200),90), 10, 4 );  //left ear
    float tri7 = triangle2( getRotatedCoords(pos, 180), getRotatedCoords(glm::vec2(98, 230.5),180),  41, 35 ); //left wing
    float tri8 = triangle2( getRotatedCoords(pos, 270), getRotatedCoords(glm::vec2(400, 230.5),270), 35, 41 ); //right wing

    float union1  = union_function( rec0, rec1, 0, 0 );
    float union2  = union_function( rec2, union1, 0, 0 );
    float union3  = union_function( rec3, union2, 0, 0 );
    float union4  = union_function( rec4, union3, 0, 0 );
    float union5  = union_function( rec5, union4, 0, 0 );
    float union6  = union_function( tri0, union5, 0, 0 );
    float union7  = union_function( tri1, union6, 0, 0 );
    float union8  = union_function( tri2, union7, 0, 0 );
    float union9  = union_function( tri3, union8, 0, 0 );
    float union10 = union_function( tri4, union9, 0, 0 );
    float union11 = union_function( tri5, union10, 0, 0 );
    float union12 = union_function( tri6, union11, 0, 0 );
    float union13 = union_function( tri7, union12, 0, 0 );
    float bat     = union_function( tri8, union13, 0, 0 );

    return bat;
}

float FRepObj2D::trebleClef(glm::vec2 pos)
{
    //constructing bottom part
    float rec_base = rectangle( pos, glm::vec2(250, 230), 15, 120 );
    float ark1     = circle( pos, glm::vec2(230, 345), 35 );
    float ark2     = circle( pos, glm::vec2(200, 325), 40 );
    float ark3     = circle( pos, glm::vec2(220, 345), 35 );
    float ark4     = circle( pos, glm::vec2(230, 320), 35 );
    float circ_bot = circle( pos, glm::vec2(210, 335), 20 );

    float arkB1 = subtract_function( ark1, ark2, 0.0f, 0.0f );
    float arkB2 = subtract_function( ark3, ark4, 0.0f, 0.0f );

    float union1  = union_function( arkB1, rec_base, 0.0f, 0.0f );
    float union2  = union_function( arkB2, union1, 0.0f, 0.0f );
    float fin_bot = union_function( circ_bot, union2, 0.0f, 0.0f );

    //constructing middle part
    float ark5 = circle( pos, glm::vec2(262, 240), 60 );
    float ark6 = circle( pos, glm::vec2(245, 245), 49 );
    float ark7 = circle( pos, glm::vec2(250, 222), 80 );
    float ark8 = circle( pos, glm::vec2(265, 222), 70 );

    float arkM1 = subtract_function( ark5, ark6, 0.0f, 0.0f );
    float arkM2 = subtract_function( ark7, ark8, 0.0f, 0.0f );

    float union3     = union_function( arkM1, fin_bot, 0.0f, 0.0f );
    float fin_middle = union_function( arkM2, union3, 0.0f, 0.0f );

    //constructing upper part
    float ark9  = circle( pos, glm::vec2(250, 95), 57 );
    float ark10 = circle( pos, glm::vec2(220, 81), 65 );
    float ark11 = circle( pos, glm::vec2(294, 100), 60 );
    float ark12 = circle( pos, glm::vec2(337, 90), 75 );
    float rec_helper1 = rectangle( pos, glm::vec2(250, 163), 100, 50 );
    float circ_top    = circle( pos, glm::vec2(277, 49), 5 );

    float arkU1  = subtract_function( ark9, ark10, 0.0f, 0.0f );
    float arkU2  = subtract_function( ark11, ark12, 0.0f, 0.0f );
    float arkU2f = subtract_function( arkU2, rec_helper1, 0.0f, 0.0f );

    float union4 = union_function( arkU1, fin_middle, 0.0f, 0.0f );
    float union5 = union_function( circ_top, union4, 0.0f, 0.0f );
    float trebleclef = union_function( arkU2f, union5, 0.0f, 0.0f );

    return trebleclef;
}

float FRepObj2D::bounded_blending(float f1, float f2, float a0, float a1, float a2, float a3,
                                  float time, float alpha, float m)
{
      f1 = intersect_function( f1, -time, alpha, m);
      f2 = intersect_function( f2, time - 1.0f, alpha, m);
      float f3 = intersect_function(time + 10.0f, 10.0f - time, alpha, m);

      float r1 = (f1/a1)*(f1/a1) + (f2/a2)*(f2/a2);
      float r2 = 0.0f;

      if( f3 > 0.0f )
      {
        r2 = (f3/a3) * (f3/a3);
      }

      float rr = 0.0f;
      if( r1 > 0.0f )
      {
        rr = r1 / (r1 + r2);
      }

      float d = 0.0f;
      if( rr < 1.0f )
      {
        d = a0 * (1.0f - rr)*(1.0f - rr)*(1.0f - rr) / (1.0f + rr);
      }

      return union_function(f1, f2, alpha, m) + d;
}

float FRepObj2D::intersect_function_R0(float f1, float f2, glm::vec2 gradf1, glm::vec2 gradf2, float n)
{
    float result = 0.0f;

    float k = 2.0f * ( gradf1.x * gradf2.x + gradf1.y * gradf1.y );

    float f12  = f1*f1;
    float f22  = f2*f2;
    float f13  = f12*f1;
    float f23  = f22*f2;

    float znam  = f1*f1 + f2*f2;
    float znam2 = znam  * znam;
    float znam3 = znam2 * znam;

    float A1 = ( f13*f1*f22 + f12*f23*f2 + f1*f2*2.0f ) / znam3;
    float B1 = ( - f12*f22*4.0f - f13*f2*2.0f - f1*f23*2.0f ) / znam2;
    float C1 = ( f12 + f22 + f1*f2*2.0f ) / znam;

    float alpha1 = 1.0f/std::sqrt( A1+B1+C1 );
    float alpha2 = 1.0f/std::sqrt( (f12+f22+f1*f2*2.0f) / znam );

    if( f1 > 0.0f && f2 > 0.0f)
       result = alpha1 * f1*f2 * std::pow( std::pow(f1, n) + std::pow(f2, n), -1.0f/n );
    else if( f1 <= 0.0f && f2 >= 0.0f )
       result = f1;
    else if( f1 >= 0.0f && f2 <= 0.0f )
       result = f2;
    else if( f1 < 0.0f && f2 < 0.0f )
       result = alpha2 * std::pow( -1.0f, n+1.0f ) * std::pow( std::pow(f1, n) + std::pow(f2, n), 1.0f/n );

    return result;
}


float FRepObj2D::union_function_R0(float f1, float f2, float n)
{
    float result = 0.0f;

    if( f1 > 0.0f && f2 > 0.0f)
       result = std::pow( std::pow(f1, n) + std::pow(f2, n), 1.0f/n );
    else if( f1 <= 0.0f && f2 >= 0.0f )
       result = f2;
    else if( f1 >= 0.0f && f2 <= 0.0f )
       result = f1;
    else if( f1 < 0.0f && f2 < 0.0f )
       result = std::pow( -1.0f, n+1.0f ) * f1*f2 * std::pow( std::pow(f1, n) + std::pow(f2, n), -1.0f/n );

    return result;
}

glm::vec2 FRepObj2D::findZeroLevelSetInterval(std::vector<float> field, int numElemToAverage)
{
    std::vector<float> posVals, negVals;

    for(size_t i = 0; i < field.size(); i++)
    {
        if( field[i] >= 0 )
            posVals.push_back(field[i]);
    }

    std::vector<float> minPosV;
    for(int j = 0; j < numElemToAverage; j++)
    {
        std::vector<float>::iterator minPos = std::min_element( posVals.begin(), posVals.end() );
        minPosV.push_back(*minPos);
        posVals.erase(std::remove( posVals.begin(), posVals.end(), *minPos), posVals.end());
    }

    float minPos, minPsum = 0;
    for(size_t i = 0; i < minPosV.size(); i++)
        minPsum += minPosV[i];

    minPos = minPsum / static_cast<float>(minPosV.size());

    glm::vec2 result = glm::vec2( 0.0f, minPos );
    return result;
}

float FRepObj2D::union_function_R0(float f1, float f2, glm::vec2 gradf1, glm::vec2 gradf2, float n)
{
    float result = 0.0f;

    float k = 2.0f * ( gradf1.x * gradf2.x + gradf1.y * gradf2.y );

    float f12  = f1*f1;
    float f22  = f2*f2;
    float f13  = f12*f1;
    float f23  = f22*f2;

    float znam  = f12 + f22;
    float znam2 = znam  * znam;
    float znam3 = znam2 * znam;

    float A1 = (  f13*f1*f22   + f12*f23*f2 + k*f13*f23 ) / znam3;
    float B1 = ( -4.0f*f12*f22 - k*f13*f2   - k*f1*f23  ) / znam2;
    float C1 = ( f12 + f22 + k*f1*f2 ) / znam;

    float alpha2 = 1.0f / std::sqrt( A1 + B1 + C1 );

    if( f1 > 0.0f && f2 > 0.0f)
    {
       result = std::pow( std::pow(f1, n) + std::pow(f2, n), 1.0f/n );
    }
    else if( f1 <= 0.0f && f2 >= 0.0f )
       result = f2;
    else if( f1 >= 0.0f && f2 <= 0.0f )
       result = f1;
    else if( f1 < 0.0f && f2 < 0.0f )
    {
        if( !std::isnan( alpha2 ) && !std::isinf( alpha2 ))
            result = alpha2 * std::pow( -1.0f, n+1.0f ) * f1*f2 * std::pow( std::pow(f1, n) + std::pow(f2, n), -1.0f/n );
        else
            result = std::pow( -1.0f, n+1.0f ) * f1*f2 * std::pow( std::pow(f1, n) + std::pow(f2, n), -1.0f/n );
    }

    return result;
}

std::vector<float> FRepObj2D::getFRep2D(std::function<float(glm::vec2)> fun)
{
    frep.clear();
    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( glm::vec2( x, y )) );
        }
    }
    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, std::function<float( glm::vec2, glm::vec2 )> fun)
{
    frep.clear();
    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( glm::vec2( x, y ), cent ));
        }
    }
    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, float R,
                                        std::function<float(glm::vec2, glm::vec2, float)> fun )
{
    frep.clear();
    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( glm::vec2( x, y ), cent, R ));
        }
    }
    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, float p1, float p2,
                                        std::function<float(glm::vec2, glm::vec2, float, float )> fun )
{
    frep.clear();
    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( glm::vec2( x, y ), cent, p1, p2 ));
        }
    }
    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(glm::vec2 cent, float p1, float p2, float p3,
                                        std::function<float(glm::vec2, glm::vec2, float, float, float)> fun)
{
    frep.clear();
    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( glm::vec2( x, y ), cent, p1, p2, p3 ));
        }
    }
    return frep;
}

std::vector<float> FRepObj2D::getFRep2D( std::vector<float> f1, std::vector<float> f2, float alpha, float m,
                                         std::function<float (float, float, float, float)> fun)
{
    frep.clear();
    for(int i = 0; i < resolutionX*resolutionY; i++)
    {
        frep.push_back( fun( f1[i], f2[i], alpha, m ));
    }

    return frep;
}

std::vector<float> FRepObj2D::getFRep2D(std::vector<float> f1, std::vector<float> f2, float n, std::function<float (float, float, float)> fun)
{
    frep.clear();
    for(int i = 0; i < resolutionX*resolutionY; i++)
    {
        frep.push_back( fun( f1[i], f2[i], n ));
    }

    return frep;
}

glm::vec2 FRepObj2D::getRotatedCoords(glm::vec2 inCoords, const float angle )
{
    float angle0 = static_cast<float>( (angle/180.0f) * M_PI );

    float rotX = inCoords.x*cosf( angle0 ) - inCoords.y*sinf( angle0);
    float rotY = inCoords.x*sinf( angle0 ) + inCoords.y*cosf( angle0 );
    return glm::vec2( rotX, rotY );
}

std::vector<float> FRepObj2D::getRotatedFrep2D( glm::vec2 cent, float w, float h,
                                                float angle, std::function<float (glm::vec2, glm::vec2, float, float)> fun)
{
    frep.clear();
    glm::vec2 c0 = getRotatedCoords( cent, angle );

    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( getRotatedCoords( glm::vec2(x,y), angle ), c0, w, h ));
        }
    }
    return frep;
}

std::vector<float> FRepObj2D::getRotatedFrep2D( glm::vec2 cent, float a, float b, float c, float angle,
                                                std::function<float (glm::vec2, glm::vec2, float, float, float, float)> fun )
{
    frep.clear();
    glm::vec2 c0 = getRotatedCoords( cent, angle );

    for( int y = 0; y < resolutionY; y++ )
    {
        for( int x = 0; x < resolutionX; x++ )
        {
            frep.push_back( fun( getRotatedCoords( glm::vec2(x, y), angle ), c0, a, b, c, angle ));
        }
    }
    return frep;
}

} //namespace frep2D_object

//headers for triangulation of the zero-level set
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/draw_triangulation_3.h>
//***********************************************
#include <fstream>

//defining new types (implicit function triangulation)
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef CGAL::Delaunay_triangulation_3<CGAL::Exact_predicates_inexact_constructions_kernel> DT3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

FT sphereFunc(Point_3 p)
{
    /*const FT x2 = (p.x() - c.x())*(p.x() - c.x()),
             y2 = (p.y() - c.y())*(p.y() - c.y()),
             z2 = (p.z() - c.z())*(p.z() - c.z());*/
    const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
    return x2+y2+z2 - 1;
}

int main(int argc, char *argv[])
{
    //triangulation of the zero-level set of the implicitly defined object
    Tr tr;
    C2t3 c2t3(tr);

    float R = 1.5f;
    Surface_3 surface(sphereFunc, Sphere_3(CGAL::ORIGIN, R*R));
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.0, 0.1, 0.1);
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

    DT3 dt3(tr.points_begin(), tr.points_end());
    CGAL::draw(dt3);
    return 0;
}

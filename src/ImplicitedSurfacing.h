#ifndef IMPLICITEDSURFACING_H
#define IMPLICITEDSURFACING_H

#include "Polygonizer.h"
#include "mymesh/readers.h"

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;



class Surfacer{

public:

    R3Pt							s_ptMin, s_ptMax;
    Array<R3Pt>					s_aptSurface;
    Array<R3Vec>					s_avecSurface;
    Array<R3Pt_i>				s_afaceSurface;


    vector<double>all_v;
    vector<uint>all_fv;

    R3Pt st;
    double dSize;
    int iBound;


    Surfacer(){}

    void CalSurfacingPara(vector<double>&Vs, vector<int>&labels);

    double Surfacing_Implicit(vector<double>&Vs, vector<int>&labels, bool ischeckall,
                   double (*function)(const R3Pt &in_pt));



    double Surfacing_CGALImplicit(FT (*function)(const Point_3 in_pt));

    double Surfacing_Spectral(vector<double>&v, vector<double> &vn);

    double Surfacing_Poisson(vector<double>&v,vector<double>&vn);

    double Surfacing_PowerCrust(vector<double>&v);

    double Surfacing_ScreenedPoisson(vector<double>&v,vector<double>&vn);

    double Surfacing_NASR(vector<double>&v);


    void WriteSurface(string fname);
    void WriteSurface(vector<double> &v, vector<uint>&fv);
    void WriteSurface(vector<double> **v, vector<uint>**fv);

    void ClearBuffer();

    void ClearSingleComponentBuffer();

private:
    void GetCurSurface(vector<double> &v, vector<uint>&fv);
    void InsertToCurSurface(vector<double>&v,vector<uint>&fv);





};


















#endif // IMPLICITEDSURFACING_H

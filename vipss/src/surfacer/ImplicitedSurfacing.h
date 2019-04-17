#ifndef IMPLICITEDSURFACING_H
#define IMPLICITEDSURFACING_H

#include "Polygonizer.h"
#include "../readers.h"


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

    void CalSurfacingPara(vector<double>&Vs, int nvoxels);

    double Surfacing_Implicit(vector<double>&Vs, int n_voxels, bool ischeckall,
                   double (*function)(const R3Pt &in_pt));



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

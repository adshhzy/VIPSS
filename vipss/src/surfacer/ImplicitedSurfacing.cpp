#include "ImplicitedSurfacing.h"


typedef std::chrono::high_resolution_clock Clock;


static Surfacer *p_ImplicitSurfacer;

static int TriProc(int in_i1, int in_i2, int in_i3, VERTICES vs) {
    const R3Pt pt = vs.ptr[in_i1].position;
    p_ImplicitSurfacer->s_afaceSurface.addItem( R3Pt_i( in_i3, in_i2, in_i1 ) );
    return 1;
}

static void VertProc(VERTICES vs) {
    p_ImplicitSurfacer->s_aptSurface.need( vs.count );
    p_ImplicitSurfacer->s_avecSurface.need( vs.count );
    for ( int i = 0; i < vs.count; i++ ) {
        p_ImplicitSurfacer->s_aptSurface[i] = vs.ptr[i].position;
        p_ImplicitSurfacer->s_avecSurface[i] = vs.ptr[i].normal;
    }
}



void Surfacer::CalSurfacingPara(vector<double>&Vs, int nvoxels){

    vector<double>leftcorner(3, DBL_MAX);
    vector<double>rightcorner(3, DBL_MIN);
    vector<double>midpoint(3,0);

    int nv = Vs.size()/3;
    for(int i=0;i<nv;++i){
        auto p_v = Vs.data()+i*3;
        for(int j=0;j<3;++j){
            leftcorner[j] = min(leftcorner[j],p_v[j]);
            rightcorner[j] = max(rightcorner[j],p_v[j]);
            //midpoint[j] += p_v[j];
        }
    }
    for(int j=0;j<3;++j){
        //midpoint[j] /= nv;
        midpoint[j] = (leftcorner[j] + rightcorner[j]) / 2.;
    }

    //double width = MyUtility::_VerticesDistance(leftcorner.data(),rightcorner.data());
    double width = -1;
    for(int j=0;j<3;++j)width = max(width, fabs(rightcorner[j]-leftcorner[j]));
    //dSize = width * 0.02;
    dSize = width * (1./nvoxels);

    if(0){

    }else{
        for(int j=0;j<3;++j){
            st[j] = midpoint[j];
        }
        iBound = (int) (width / dSize / 2. * 1.75);

    }


}


double Surfacer::Surfacing_Implicit(vector<double>&Vs,int n_voxels, bool ischeckall,
                                    double (*function)(const R3Pt &in_pt)){

    p_ImplicitSurfacer = this;
    ClearBuffer();

    CalSurfacingPara(Vs, n_voxels);


    //polygonize(function, size, bounds, st, triproc, vertproc);
    double thresDist = 1e-3;

    double re_time;
    cout<<"Implicit Surfacing: "<<endl;

    auto t1 = Clock::now();




    if(!ischeckall){
        polygonize(function, dSize, iBound, st, TriProc, VertProc);
        GetCurSurface(all_v,all_fv);
    }else{
        cout << "Deprecated, Please seek to CGAL implicit surfacer for surfaces with multiple connected components" <<endl;
    }

    cout<<"Implicit Surfacing Done."<<endl;
    auto t2 = Clock::now();
    cout << "Total Surfacing time: " <<  (re_time = std::chrono::nanoseconds(t2 - t1).count()/1e9) <<endl;

    return re_time;

}


void Surfacer::WriteSurface(string fname){

    writeObjFile(fname,all_v,all_fv);

}

void Surfacer::WriteSurface(vector<double> &v, vector<uint>&fv){

    v = all_v;
    fv = all_fv;
}

void Surfacer::WriteSurface(vector<double> **v, vector<uint> **fv){

    *v = &all_v;
    *fv = &all_fv;
}

void Surfacer::ClearBuffer(){

    ClearSingleComponentBuffer();
    all_v.clear();
    all_fv.clear();
}

void Surfacer::ClearSingleComponentBuffer(){
    s_aptSurface.clearcompletely();
    s_avecSurface.clearcompletely();
    s_afaceSurface.clearcompletely();
}

void Surfacer::GetCurSurface(vector<double>&v,vector<uint>&fv){

    int beInd = v.size()/3;
    for(int i=0;i<s_aptSurface.num();++i){
        for(int j=0;j<3;++j)v.push_back( s_aptSurface[i][j] );
    }

    for(int i=0;i<s_afaceSurface.num();++i){
        for(int j=0;j<3;++j)fv.push_back( beInd + s_afaceSurface[i][2-j]);
    }

}


void Surfacer::InsertToCurSurface(vector<double>&v,vector<uint>&fv){

    int beInd = all_v.size()/3;
    all_v.insert(all_v.end(),v.begin(),v.end());

    for(auto a:fv)all_fv.push_back(beInd+a);

}

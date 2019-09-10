#include "rbfcore.h"
#include "utility.h"
#include "Solver.h"
#include <armadillo>
#include <fstream>
#include <limits>
#include <iomanip>
#include <ctime>
#include <chrono>
#include<algorithm>



double sigma = 2.0;
double inv_sigma_squarex2 = 1/(2 * pow(sigma, 2));
double Gaussian_Kernel(const double x_square){

    return exp(-x_square*inv_sigma_squarex2);

}

double Gaussian_Kernel_2p(const double *p1, const double *p2){



    return Gaussian_Kernel(MyUtility::vecSquareDist(p1,p2));


}

double Gaussian_PKernel_Dirichlet_2p(const double *p1, const double *p2){


    double d2 = MyUtility::vecSquareDist(p1,p2);
    return (6*sigma*sigma-d2)*sqrt(Gaussian_Kernel(d2));


}

double Gaussian_PKernel_Bending_2p(const double *p1, const double *p2){


    double d2 = MyUtility::vecSquareDist(p1,p2);
    double d4 = d2*d2;
    double sigma2 = sigma * sigma;
    double sigma4 = sigma2 * sigma2;
    return (60*sigma4-20*sigma2*d2+d4)*sqrt(Gaussian_Kernel(d2));


}


double XCube_Kernel(const double x){

    return pow(x,3);
}

double XCube_Kernel_2p(const double *p1, const double *p2){


    return XCube_Kernel(MyUtility::_VerticesDistance(p1,p2));

}

void XCube_Gradient_Kernel_2p(const double *p1, const double *p2, double *G){


    double len_dist  = MyUtility::_VerticesDistance(p1,p2);
    for(int i=0;i<3;++i)G[i] = 3*len_dist*(p1[i]-p2[i]);
    return;

}

double XCube_GradientDot_Kernel_2p(const double *p1, const double *p2, const double *p3){


    double G[3];
    XCube_Gradient_Kernel_2p(p1,p2,G);
    return MyUtility::dot(p3,G);

}

void XCube_Hessian_Kernel_2p(const double *p1, const double *p2, double *H){


    double diff[3];
    for(int i=0;i<3;++i)diff[i] = p1[i] - p2[i];
    double len_dist  = sqrt(MyUtility::len(diff));

    if(len_dist<1e-8){
        for(int i=0;i<9;++i)H[i] = 0;
    }else{
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            if(i==j)H[i*3+j] = 3 * pow(diff[i],2) / len_dist + 3 * len_dist;
            else H[i*3+j] = 3 * diff[i] * diff[j] / len_dist;
    }


    return;

}

void XCube_HessianDot_Kernel_2p(const double *p1, const double *p2, const double *p3, vector<double>&dotout){


    double H[9];
    XCube_Gradient_Kernel_2p(p1,p2,H);
    dotout.resize(3);
    for(int i=0;i<3;++i){
        dotout[i] = 0;
        for(int j=0;j<3;++j){
            dotout[i] += H[i*3+j] * p3[j];
        }
    }

}

RBF_Core::RBF_Core(){

    Kernal_Function = Gaussian_Kernel;
    Kernal_Function_2p = Gaussian_Kernel_2p;
    P_Function_2p = Gaussian_PKernel_Dirichlet_2p;

    isHermite = false;

    mp_RBF_INITMETHOD.insert(make_pair(GT_NORMAL,"GT_NORMAL"));
    mp_RBF_INITMETHOD.insert(make_pair(GlobalEigen,"GlobalEigen"));
    mp_RBF_INITMETHOD.insert(make_pair(GlobalEigenWithMST,"GlobalEigenWithMST"));
    mp_RBF_INITMETHOD.insert(make_pair(GlobalEigenWithGT,"GlobalEigenWithGT"));
    mp_RBF_INITMETHOD.insert(make_pair(LocalEigen,"LocalEigen"));
    mp_RBF_INITMETHOD.insert(make_pair(IterativeEigen,"IterativeEigen"));
    mp_RBF_INITMETHOD.insert(make_pair(ClusterEigen,"ClusterEigen"));


    mp_RBF_METHOD.insert(make_pair(Variational,"Variational"));
    mp_RBF_METHOD.insert(make_pair(Variational_P,"Variation_P"));
    mp_RBF_METHOD.insert(make_pair(LS,"LS"));
    mp_RBF_METHOD.insert(make_pair(LSinterp,"LSinterp"));
    mp_RBF_METHOD.insert(make_pair(Interp,"Interp"));
    mp_RBF_METHOD.insert(make_pair(RayleighQuotients,"Rayleigh"));
    mp_RBF_METHOD.insert(make_pair(RayleighQuotients_P,"Rayleigh_P"));
    mp_RBF_METHOD.insert(make_pair(RayleighQuotients_I,"Rayleigh_I"));
    mp_RBF_METHOD.insert(make_pair(Hermite,"Hermite"));
    mp_RBF_METHOD.insert(make_pair(Hermite_UnitNorm,"UnitNorm"));
    mp_RBF_METHOD.insert(make_pair(Hermite_UnitNormal,"UnitNormal"));
    mp_RBF_METHOD.insert(make_pair(Hermite_Tangent_UnitNorm,"T_UnitNorm"));
    mp_RBF_METHOD.insert(make_pair(Hermite_Tangent_UnitNormal,"T_UnitNormal"));

    mp_RBF_Kernal.insert(make_pair(XCube,"TriH"));
    mp_RBF_Kernal.insert(make_pair(ThinSpline,"ThinSpline"));
    mp_RBF_Kernal.insert(make_pair(XLinear,"XLinear"));
    mp_RBF_Kernal.insert(make_pair(Gaussian,"Gaussian"));

}
RBF_Core::RBF_Core(RBF_Kernal kernal){
    isHermite = false;
    Init(kernal);
}

void RBF_Core::Init(RBF_Kernal kernal){

    this->kernal = kernal;
    switch(kernal){
    case Gaussian:
        Kernal_Function = Gaussian_Kernel;
        Kernal_Function_2p = Gaussian_Kernel_2p;
        P_Function_2p = Gaussian_PKernel_Dirichlet_2p;
        break;

    case XCube:
        Kernal_Function = XCube_Kernel;
        Kernal_Function_2p = XCube_Kernel_2p;
        Kernal_Gradient_Function_2p = XCube_Gradient_Kernel_2p;
        Kernal_Hessian_Function_2p = XCube_Hessian_Kernel_2p;
        break;

    default:
        break;

    }

}

void RBF_Core::SetSigma(double x){
    sigma = x;
    inv_sigma_squarex2 = 1/(2 * pow(sigma, 2));
}

double RBF_Core::Dist_Function(const double x, const double y, const double z){


	return -1;

}



inline double RBF_Core::Dist_Function(const double *p){

    n_evacalls++;
    double *p_pts = pts.data();
    static arma::vec kern(npt), kb;
    if(isHermite){
        kern.set_size(npt*4);
        double G[3];
        for(int i=0;i<npt;++i)kern(i) = Kernal_Function_2p(p_pts+i*3, p);
        for(int i=0;i<npt;++i){
            Kernal_Gradient_Function_2p(p,p_pts+i*3,G);
            //for(int j=0;j<3;++j)kern(npt+i*3+j) = -G[j];
            for(int j=0;j<3;++j)kern(npt+i+j*npt) = G[j];
        }
    }else{
        kern.set_size(npt);
        for(int i=0;i<npt;++i)kern(i) = Kernal_Function_2p(p_pts+i*3, p);
    }

    double loc_part = dot(kern,a);

    if(polyDeg==1){
        kb.set_size(4);
        for(int i=0;i<3;++i)kb(i+1) = p[i];
        kb(0) = 1;
    }else if(polyDeg==2){
        vector<double>buf(4,1);
        int ind = 0;
        kb.set_size(10);
        for(int j=0;j<3;++j)buf[j+1] = p[j];
        for(int j=0;j<4;++j)for(int k=j;k<4;++k)kb(ind++) = buf[j] * buf[k];
    }
    double poly_part = dot(kb,b);

    if(0){
        cout<<"dist: "<<p[0]<<' '<<p[1]<<' '<<p[2]<<' '<<p_pts[3]<<' '<<p_pts[4]<<' '<<p_pts[5]<<' '<<
              Kernal_Function_2p(p,p_pts+3)<<endl;
        for(int i=0;i<npt;++i)cout<<kern(i)<<' ';
        for(int i=0;i<bsize;++i)cout<<kb(i)<<' ';
        cout<<endl;
    }

    double re = loc_part + poly_part;
    return re;


}
static RBF_Core * s_hrbf;
double RBF_Core::Dist_Function(const R3Pt &in_pt){
    return s_hrbf->Dist_Function(&(in_pt[0]));
}

//FT RBF_Core::Dist_Function(const Point_3 in_pt){

//    return s_hrbf->Dist_Function(&(in_pt.x()));
//}

void RBF_Core::SetThis(){

    s_hrbf = this;
}

void RBF_Core::Write_Surface(string fname){

    //writeObjFile(fname,finalMesh_v,finalMesh_fv);

    writePLYFile_VF(fname,finalMesh_v,finalMesh_fv);
}

/**********************************************************/


void RBF_Core::Record(RBF_METHOD method, RBF_Kernal kernal, Solution_Struct &rsol, double time){

    npoints.push_back(npt);

    record_initmethod.push_back(mp_RBF_INITMETHOD[curInitMethod]);
    record_method.push_back(mp_RBF_METHOD[method]);
    record_kernal.push_back(mp_RBF_Kernal[kernal]);
    record_initenergy.push_back(rsol.init_energy);
    record_energy.push_back(rsol.energy);
    record_time.push_back(time);


    setup_timev.push_back(setup_time);
    init_timev.push_back(init_time);
    solve_timev.push_back(solve_time);
    callfunc_timev.push_back(callfunc_time);
    invM_timev.push_back(invM_time);
    setK_timev.push_back(setK_time);

}

void RBF_Core::Record(){

    //cout<<"record"<<endl;
    npoints.push_back(npt);

    record_initmethod.push_back(mp_RBF_INITMETHOD[curInitMethod]);
    record_method.push_back(mp_RBF_METHOD[curMethod]);
    //record_kernal.push_back(mp_RBF_Kernal[kernal]);
    record_initenergy.push_back(sol.init_energy);
    record_energy.push_back(sol.energy);
    //cout<<"record"<<endl;


//    setup_timev.push_back(setup_time);
//    init_timev.push_back(init_time);
//    solve_timev.push_back(solve_time);
//    callfunc_timev.push_back(callfunc_time);
//    invM_timev.push_back(invM_time);
//    setK_timev.push_back(setK_time);
   // cout<<"record end"<<endl;
}

void RBF_Core::AddPartition(string pname){

    record_partition.push_back(record_method.size());
    record_partition_name.push_back(pname);
}


void RBF_Core::Print_Record(){

    cout<<"Method\t\t Kernal\t\t Energy\t\t Time"<<endl;
    cout<<std::setprecision(8)<<endl;
    if(record_partition.size()==0){
        for(int i=0;i<record_method.size();++i){
            cout<<record_method[i]<<"\t\t"<<record_kernal[i]<<"\t\t"<<record_energy[i]<<"\t\t"<<record_time[i]<<endl;
        }
        for(int i=0;i<setup_timev.size();++i){
            cout<<setup_timev[i]<<"\t\t"<<init_timev[i]<<"\t\t"<<solve_timev[i]<<"\t\t"<<callfunc_timev[i]<<"\t\t"<<invM_timev[i]<<"\t\t"<<setK_timev[i]<<endl;
        }
    }else{
        for(int j=0;j<record_partition.size();++j){
            cout<<record_partition_name[j]<<endl;
            for(int i=j==0?0:record_partition[j-1];i<record_partition[j];++i){
                cout<<record_method[i]<<"\t\t"<<record_kernal[i]<<"\t\t"<<record_energy[i]<<"\t\t"<<record_time[i]<<endl;
            }
            for(int i=j==0?0:record_partition[j-1];i<record_partition[j];++i){
                cout<<setup_timev[i]<<"\t\t"<<init_timev[i]<<"\t\t"<<solve_timev[i]<<"\t\t"<<callfunc_timev[i]<<endl;
            }
        }
    }

}


void RBF_Core::Print_TimerRecord(string fname){

    ofstream fout(fname);
    fout<<setprecision(5);
    if(!fout.fail()){
        for(int i=0;i<setup_timev.size();++i){
            fout<<npoints[i]<<'\t'<<setup_timev[i]<<"\t"<<init_timev[i]<<"\t"<<solve_timev[i]<<"\t"<<callfunc_timev[i]<<"\t"<<invM_timev[i]<<"\t"<<setK_timev[i]<<endl;
        }
    }
    fout.close();

}

void RBF_Core::Print_TimerRecord_Single(string fname){

    ofstream fout(fname);
    fout<<setprecision(5);
    if(!fout.fail()){
        fout<<"number of points: "<<npt<<endl
           <<"setup_time (Compute H): "<<setup_time<<" s"<<endl
          <<"init_time (Optimize g/Eigen): "<<init_time<<" s"<<endl
         <<"solve_time (Optimize g/LBFGS): "<<solve_time<<" s"<<endl
        <<"surfacing_time: "<<surf_time<<" s"<<endl;
    }
    fout.close();
}

void RBF_Core::Clear_TimerRecord(){
    npoints.clear();
    setup_timev.clear();
    init_timev.clear();
    solve_timev.clear();
    callfunc_timev.clear();
    invM_timev.clear();
    setK_timev.clear();

}

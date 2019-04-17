#include "rbfcore.h"
#include "mymesh/utility.h"
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



double RBF_Core::Dist_Function(const double *p){

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

FT RBF_Core::Dist_Function(const Point_3 in_pt){

    return s_hrbf->Dist_Function(&(in_pt.x()));
}

void RBF_Core::SetThis(){

    s_hrbf = this;
}

void RBF_Core::Set_M(){

    a.randu(npt);
    M.randu(npt,npt);
    double *p_pts = pts.data();
    for(int i=0;i<npt;++i){
        for(int j=i;j<npt;++j){
            M(i,j) = M(j,i) = Kernal_Function_2p(p_pts+i*3, p_pts+j*3);
        }
    }

    //    cout<<"setM "<<p_pts[0]<<' '<<p_pts[1]<<' '<<p_pts[2]<<' '<<p_pts[3]<<' '<<p_pts[4]<<' '<<p_pts[5]<<' '<<
    //          Kernal_Function_2p(p_pts, p_pts+3)<<endl;
    //    for(int i=0;i<npt;++i)cout<<pts[i]<<' ';
    //    cout<<endl;
    //cout<<M.row(0)<<endl;

}

void RBF_Core::Set_N(){


    if(polyDeg==1){
        bsize= 4;
        N.randu(npt,4);
        b.randu(4);

        for(int i=0;i<npt;++i){

            N(i,0) = 1;
            for(int j=0;j<3;++j)N(i,j+1) = pts[i*3+j];

        }
    }else if(polyDeg==2){
        bsize = 10;
        N.randu(npt,10);
        b.randu(10);
        vector<double>buf(4,1);

        for(int i=0;i<npt;++i){

            for(int j=0;j<3;++j)buf[j+1] = pts[i*3+j];
            int ind = 0;
            for(int j=0;j<4;++j)for(int k=j;k<4;++k)N(i,ind++) = buf[j] * buf[k];

        }
    }

}

void RBF_Core::Set_Minv(){

    if(isinv)Minv = inv_sympd(M);
    else {
        arma::mat Eye;
        Eye.eye(npt,npt);
        Minv = solve(M,Eye);
    }
    //cout<<(M*Minv)<<endl;
}
void RBF_Core::Set_Bprey(){

    if(isinv)bprey = inv_sympd(N.t() * Minv * N) * N.t() * Minv;
    else {
        arma::mat Eye2;
        Eye2.eye(bsize,bsize);
        bprey = solve(N.t() * Minv * N, Eye2) * N.t() * Minv;
    }
}
void RBF_Core::Set_VariationalRBF(vector<double>&pts, vector<int>&labels){

    //isinv = false;


    cout<<npt<<endl;
    cout<<"Building up matrix"<<endl;

    Set_M();
    Set_N();

    Set_Minv();
    Set_Bprey();



    arma::mat Eye;
    Eye.eye(npt,npt);
    arma::mat c = N * bprey;
    arma::mat IminusC =  Eye - c;



    if(kernal == Gaussian){
        P.randu(npt,npt);

        double *p_pts = pts.data();
        for(int i=0;i<npt;++i){
            for(int j=i;j<npt;++j){
                P(i,j) = P(j,i) = P_Function_2p(p_pts+i*3, p_pts+j*3);
            }
        }
        arma::mat MinvIminusC = Minv * IminusC;
        K = MinvIminusC.t() * P * MinvIminusC;
    }else{
        cout<<"appro"<<endl;

        P = M;
        K = IminusC.t() * Minv.t() * IminusC;

    }

    arma::vec eigval = eig_sym( K ) ;
    K += (0.1-eigval(0))*Eye;
    //cout<<eigval<<endl;


    cout<<"Done"<<endl;

}



int RBF_Core::Solve_VariationalRBF(){



    vector<LinearVec>LVs(npt);

    for(int i=0;i<npt;++i){
        LVs[i].push_back(i,1.);
        LVs[i].set_label(labels[i]);
        if(labels[i]==0)LVs[i].set_ab(0., 0.);
        else if(labels[i]==-1)LVs[i].set_ab(-maxvalue, -rangevalue);
        else if(labels[i]== 1)LVs[i].set_ab(rangevalue, maxvalue);
    }
    Solver::solveQuadraticProgramming(K,LVs,npt,sol);

    arma::vec y(npt);
    for(int i=0;i<npt;++i)y(i) = sol.solveval[i];

    //cout<<y<<endl;
    //cout<<y.t()<<endl;
    //cout<<y<<endl;

    b = bprey * y;
    a = Minv * (y - N*b);

    //cout<<a.t()<<endl;
    //cout<<b.t()<<endl;


    cout<<endl;
    double epsilon_test = 1e-5;
    if(0)for(int i=0;i<npt;++i){
        if(labels[i]==0){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( fabs(aa) > epsilon_test){
                cout<<i<<" Wrong Solution l0: "<<aa<<endl;
            }
        }else if(labels[i]==-1){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( aa<=-maxvalue || aa >= -rangevalue+epsilon_test){
                cout<<i<<" Wrong Solution l-1: "<<aa<<endl;
            }
        }else if(labels[i]== 1){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( aa<=rangevalue-epsilon_test || aa >= maxvalue){
                cout<<i<<" Wrong Solution l+1: "<<aa<<endl;
            }
        }
    }
	
	return 1;


}

/**************************************************************/

void RBF_Core::Set_LSRBF(vector<double>&pts, vector<int>&labels, bool is_interp){
    //Init(2);
    Set_M();
    Set_N();

    if(is_interp){
        for(int i=0;i<npt;++i)if(labels[i]!=0)M(i,i)+=User_Lamnbda;
    }else{
        for(int i=0;i<npt;++i)M(i,i)+=User_Lamnbda;
    }

    Set_Minv();
    Set_Bprey();

}

int RBF_Core::Solve_RBF(){


    arma::vec y(npt);
    for(int i=0;i<npt;++i)y(i) = rangevalue*labels[i];

    b = bprey * y;
    a = Minv * (y - N*b);
	return 1;
}



/**************************************************************/

void RBF_Core::Set_InterpRBF(vector<double>&pts, vector<int>&labels){
    User_Lamnbda = 0;
    Set_LSRBF(pts,labels, false);
}
/************************************************************/
void RBF_Core::Set_VariationalRBF_P(vector<double>&pts, vector<int>&labels){

    //Init(1);
    //SetSigma(0.3);

    Set_M();
    Set_N();

    if(kernal == Gaussian){
        P.randu(npt,npt);

        double *p_pts = pts.data();
        for(int i=0;i<npt;++i){
            for(int j=i;j<npt;++j){
                P(i,j) = P(j,i) = P_Function_2p(p_pts+i*3, p_pts+j*3);
            }
        }
    }else{
        cout<<"appro"<<endl;
        P=M;
    }


}


int RBF_Core::Solve_VariationalRBF_P(){



    vector<LinearVec>LVs(npt+bsize);

    for(int i=0;i<npt;++i){
        for(int j=0;j<npt;++j)LVs[i].push_back(j,M(i,j));
        for(int j=0;j<bsize;++j)LVs[i].push_back(j+npt,N(i,j));
        LVs[i].set_label(labels[i]);
        if(labels[i]==0)LVs[i].set_ab(0., 0.);
        else if(labels[i]==-1)LVs[i].set_ab(-maxvalue, -rangevalue);
        else if(labels[i]== 1)LVs[i].set_ab(rangevalue, maxvalue);
    }
    for(int i=0;i<bsize;++i){
        for(int j=0;j<npt;++j)LVs[i+npt].push_back(j,N(j,i));
        LVs[i+npt].set_label(0);
        LVs[i+npt].set_ab(0., 0.);
    }


    Solver::solveQuadraticProgramming(P,LVs,npt+bsize,sol);


    a.randu(npt);
    b.randu(bsize);

    for(int j=0;j<npt;++j)a(j) = sol.solveval[j];
    for(int j=0;j<bsize;++j)b(j) = sol.solveval[j+npt];


    double epsilon_test = 1e-5;
    if(0)for(int i=0;i<npt;++i){
        if(labels[i]==0){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( fabs(aa) > epsilon_test){
                cout<<i<<" Wrong Solution l0: "<<aa<<endl;
            }
        }else if(labels[i]==-1){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( aa<=-maxvalue || aa >= -rangevalue+epsilon_test){
                cout<<i<<" Wrong Solution l-1: "<<aa<<endl;
            }
        }else if(labels[i]== 1){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( aa<=rangevalue-epsilon_test || aa >= maxvalue){
                cout<<i<<" Wrong Solution l+1: "<<aa<<endl;
            }
        }
    }

	return 1;
}


/**********************************************************/


void RBF_Core::Set_RayleighQuotients_P(vector<double>&pts, vector<int>&labels){

    Set_M();
    Set_N();

    Set_Minv();
    Set_Bprey();


    if(kernal == Gaussian){
        P.randu(npt,npt);

        double *p_pts = pts.data();
        for(int i=0;i<npt;++i){
            for(int j=i;j<npt;++j){
                P(i,j) = P(j,i) = P_Function_2p(p_pts+i*3, p_pts+j*3);
            }
        }
    }else{
        cout<<"appro"<<endl;
        P=M;
    }

    RQ.zeros(npt,npt);

    arma::mat Eye;
    Eye.eye(npt,npt);
    arma::vec eigval = eig_sym( P ) ;
    //    RQ.submat(0,0,npt-1,npt-1) = P+ (1-eigval(0))*Eye;
    RQ = P + (1-eigval(0))*Eye;

    RQ=-RQ;
    //arma::cx_vec eigval = eigs_gen( RQ, 1, "lr" );


    cout <<eigval<<endl;
    //cout<<RQ<<endl;



}




int RBF_Core::Solve_RayleighQuotients_P(){

    vector<LinearVec>LVs(npt+bsize);
    int tpos = npt+bsize;
    if(1){
        LVs.resize(tpos+1);
        LVs[tpos].push_back(tpos,1);
        LVs[tpos].set_ab(0.1,maxvalue);
        LVs[tpos].mt = CM_GREATER;
    }

    for(int i=0;i<npt;++i){
        double coef = labels[i]==1?-1:1;
        for(int j=0;j<npt;++j)LVs[i].push_back(j, coef*M(i,j));
        for(int j=0;j<bsize;++j)LVs[i].push_back(j+npt, coef*N(i,j));
        LVs[i].set_label(labels[i]);

        double td;
        if(labels[i]==0)td = 0;
        else if(labels[i]==-1)td = rangevalue;
        else if(labels[i]== 1)td = rangevalue;
        LVs[i].push_back(tpos, td);

        LVs[i].set_ab(0,0);
        LVs[i].mt = labels[i]==0? CM_EQUAL : CM_LESS;
    }
    for(int i=0;i<bsize;++i){
        for(int j=0;j<npt;++j)LVs[i+npt].push_back(j,N(j,i));
        LVs[i+npt].set_label(0);
        LVs[i+npt].set_ab(0., 0.);
        LVs[i+npt].mt = CM_EQUAL;
    }


    vector<QuadraticVec>Cb(1);
    for(int i=0;i<tpos;++i){
        Cb[0].push_back(i,i,1);
    }
    Cb[0].set_b(1);
    Cb[0].mt = CM_LESS;
    //Cb[0].mt = CM_GREATER;

    Solver::solveQCQP(RQ, LVs, Cb, npt+bsize+1, sol);

    a.randu(npt);
    b.randu(bsize);

    for(int j=0;j<npt;++j)a(j) = sol.solveval[j];
    for(int j=0;j<bsize;++j)b(j) = sol.solveval[j+npt];


    cout<<"RayleighQuotients: t: "<<sol.solveval[tpos]<<endl;

    cout<<endl;
    double epsilon_test = 1e-5;
    if(0)for(int i=0;i<npt;++i){
        if(labels[i]==0){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( fabs(aa) > epsilon_test){
                cout<<i<<" Wrong Solution l0: "<<aa<<endl;
            }
        }else if(labels[i]==-1){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( aa<=-maxvalue || aa >= -rangevalue+epsilon_test){
                cout<<i<<" Wrong Solution l-1: "<<aa<<endl;
            }
        }else if(labels[i]== 1){
            double aa = this->Dist_Function(pts.data()+i*3);
            if( aa<=rangevalue-epsilon_test || aa >= maxvalue){
                cout<<i<<" Wrong Solution l+1: "<<aa<<endl;
            }
        }
    }

	return 1;
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

    cout<<setprecision(4);
    cout<<endl;

    ofstream fout(fname);
    fout<<setprecision(5);
    if(!fout.fail()){
        for(int i=0;i<setup_timev.size();++i){
            fout<<npoints[i]<<'\t'<<setup_timev[i]<<"\t"<<init_timev[i]<<"\t"<<solve_timev[i]<<"\t"<<callfunc_timev[i]<<"\t"<<invM_timev[i]<<"\t"<<setK_timev[i]<<endl;
        }
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

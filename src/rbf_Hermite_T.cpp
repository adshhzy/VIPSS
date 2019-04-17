#include "rbfcore.h"
#include "mymesh/utility.h"
#include "Solver.h"
#include <armadillo>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <ctime>
#include <chrono>
#include <iomanip>
#include<algorithm>
#include "mymesh/readers.h"

typedef std::chrono::high_resolution_clock Clock;


bool RBF_Core::Set_Hermite_DesignedCurve(){

    if(curMethod!=Hermite_Tangent_UnitNormal)return false;
    cout<<"Hermite_designcurve_weight: "<<Hermite_designcurve_weight<<endl;

    double w = Hermite_designcurve_weight;
    if(w==0)return false;

    //saveK_finalH = K11;
    int ne = edges.size()/2;
    auto p_tangent = tangents.data();
    if(1)for(int i=0;i<ne;++i){
        int v1 = edges[i*2], v2 = edges[i*2+1];
        arma::sp_mat b1(3,npt*3),b2(3,npt*3);
        b1(0, v1+1*npt) = -p_tangent[v1*3+2];
        b1(0, v1+2*npt) = p_tangent[v1*3+1];

        b1(1, v1+0*npt) = p_tangent[v1*3+2];
        b1(1, v1+2*npt) = -p_tangent[v1*3+0];

        b1(2, v1+0*npt) = -p_tangent[v1*3+1];
        b1(2, v1+1*npt) = p_tangent[v1*3+0];

        b2(0, v2+1*npt) = -p_tangent[v2*3+2];
        b2(0, v2+2*npt) = p_tangent[v2*3+1];

        b2(1, v2+0*npt) = p_tangent[v2*3+2];
        b2(1, v2+2*npt) = -p_tangent[v2*3+0];

        b2(2, v2+0*npt) = -p_tangent[v2*3+1];
        b2(2, v2+1*npt) = p_tangent[v2*3+0];

        arma::sp_mat b3 = b1-b2, b4 = w* ( b3.t() * b3);
        K11 += b4;
//        K += b4;
//        finalH += b4;
//        if(isuse_sparse){
//            sp_K += b4; sp_H += b4;
//        }
    }
    else for(int i=0;i<ne;++i){
        int v1 = edges[i*2], v2 = edges[i*2+1];
        arma::sp_mat b(1,npt*3);
        b(0,v1+1*npt) = -p_tangent[v1*3+2];
        b(0,v1+2*npt) = p_tangent[v1*3+1];
        b(0,v2+1*npt) = p_tangent[v2*3+2];
        b(0,v2+2*npt) = -p_tangent[v2*3+1];
        K += w* ( b.t() * b);

        b.zeros();
        b(0,v1+0*npt) = p_tangent[v1*3+2];
        b(0,v1+2*npt) = -p_tangent[v1*3+0];
        b(0,v2+0*npt) = -p_tangent[v2*3+2];
        b(0,v2+2*npt) = p_tangent[v2*3+0];
        K += w* ( b.t() * b);

        b.zeros();
        b(0,v1+0*npt) = -p_tangent[v1*3+1];
        b(0,v1+1*npt) = p_tangent[v1*3+0];
        b(0,v2+0*npt) = p_tangent[v2*3+1];
        b(0,v2+1*npt) = -p_tangent[v2*3+0];
        K += w* ( b.t() * b);
    }
    //saveK = K;

    return true;
}


void RBF_Core::Set_Hermite_PredictNormal_TangentConstraint(vector<double>&pts, vector<double>&tgt, vector<uint> &ede){

    Set_Hermite_PredictNormal(pts);
    //Set_Hermite_DesignedCurve();


}
int RBF_Core::Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm(){

    arma::vec eigval, ny;
    arma::mat eigvec;
    arma::sp_mat c2(npt,npt*3);
    arma::sp_mat eye;
    auto p_tangent = tangents.data();
    for(int i=0;i<npt;++i){
        c2(i,i) = p_tangent[i*3];
        c2(i,i+1*npt) = p_tangent[i*3+1];
        c2(i,i+2*npt) = p_tangent[i*3+2];
    }
    eye.eye(npt*3,npt*3);
    arma::mat cct(c2*c2.t());

    arma::mat p = eye - c2.t()*arma::inv(cct)*c2;

    arma::mat AA = p*K*p;
    AA = (AA.t()+AA)*0.5;
    ny = eig_sym( eigval, eigvec, AA);

    //cout<<eigval.t()<<endl;
    //vector<pair<double,int>>projection;
    vector<bool>isNullVec(npt*3,false);
    for(int i=0;i<npt*3;++i){
        double val = arma::norm(c2*eigvec.col(i));
        //projection.emplace_back(val,i);
        if(val<1e-5)isNullVec[i]=true;
    }
    //sort(projection.begin(),projection.end());
    //for(int i=0;i<npt*3;++i)cout<<projection[i].first<<' ';cout<<endl<<endl;
    //for(int i=0;i<npt*3;++i)cout<<projection[i].second<<' ';cout<<endl<<endl;

    int smalleig = 0;
    for(int i=0;i<npt*3;++i)if(isNullVec[i]){smalleig = i;break;}
    //smalleig = 0;
    //ny = eig_sym( eigval, eigvec, K);
    //int smalleig = 0;


    //vector<double>dottangenttest(npt);
    initnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt*3;++i)y(i+npt) = eigvec(i,smalleig);
    for(int i=0;i<npt;++i){
        initnormals[i*3]   = y(npt+i);
        initnormals[i*3+1] = y(npt+i+npt);
        initnormals[i*3+2] = y(npt+i+npt*2);
        //MyUtility::normalize(normals.data()+i*3);
        //dottangenttest[i] = MyUtility::dot(p_tangent+i*3, normals.data()+i*3);
    }
    //for(auto a:dottangenttest)cout<<a<<' ';cout<<endl;
    //NormalRecification(1);

    //b = bprey * y;
    //a = Minv * (y - N*b);

    //sol.energy = arma::dot(a,M*a);
    SetInitnormal_Uninorm();

    cout<<"Solve_Hermite_PredictNormal_Tangent_UnitNorm"<<endl;
    return 1;
}




/*************************************************/
//double optfunc_Hermite(const vector<double>&x, vector<double>&grad, void *fdata);
static int countopt = 0;
double optfunc_Hermite_Tangents(const vector<double>&x, vector<double>&grad, void *fdata){

    RBF_Core *drbf = reinterpret_cast<RBF_Core*>(fdata);
    int n = drbf->npt;
    arma::vec arma_x(n*3);

    vector<double>&coff_cos = drbf->coff_cos;
    vector<double>&coff_sin = drbf->coff_sin;
    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
    vector<double>cosa(n);
    vector<double>sina(n);
    for(int i=0;i<n;++i){
        cosa[i] = cos(x[i]);
        sina[i] = sin(x[i]);
    }

    for(int i=0;i<n;++i){
        for(int j=0;j<3;++j){
            arma_x(i+n*j) = coff_cos[i*3+j] * cosa[i] + coff_sin[i*3+j] * sina[i];
        }
    }

    arma::vec a2 = drbf->finalH * arma_x;

    if (!grad.empty()) {
        grad.resize(n);
        for(int i=0;i<n;++i){
            grad[i] = 0;
            for(int j=0;j<3;++j)grad[i] += a2(i+n*j) * (-coff_cos[i*3+j] * sina[i] + coff_sin[i*3+j] * cosa[i]);
        }
    }
    //for(auto a:grad)cout<<a<<' ';cout<<endl;

    double re = arma::dot( arma_x, a2 );
    countopt++;
    //cout<<countopt++<<' '<<re<<endl;
    return re;

}


int RBF_Core::Solve_Hermite_PredictNormal_TangentConstraint_UnitNormal(){

    if(1){
        Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm();
        //initnormals = newnormals;
    }else{
        LocalEigenInit_PredictNormal_TangentConstraint();
        //newnormals = initnormals;
    }
    //initnormals = newnormals = normals;
    SetInitnormal_Uninorm();

    BuildLFrame(initnormals_uninorm);

    //cout<<"coanima"<<endl;
    // return -1;

    //LocalIterativeSolver(sol, newnormals, 1000, 1e-7);
    //LocalIterativeSolver_Tconstraint(sol, newnormals, tangents, 1000, 1e-7);
    sol.solveval.resize(npt);
    for(int i=0;i<npt;++i){
        sol.solveval[i] = 0.0;
    }
    //cout<<"smallvec: "<<smallvec<<endl;

    if(1){
        vector<double>upper(npt);
        vector<double>lower(npt);
        for(int i=0;i<npt;++i){
            upper[i] = 2 * my_PI;
            lower[i] = -2 * my_PI;
        }


        countopt = 0;

        //LocalIterativeSolver(sol,kk==0?normals:newnormals,300,1e-7);


        Solver::nloptwrapper(lower,upper,optfunc_Hermite_Tangents,this,1e-7,3000,sol);
        cout<<"number of call: "<<countopt<<endl;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

    }


    newnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;

    if(1){
        for(int i=0;i<npt;++i){
            double cosa = cos(sol.solveval[i]), sina = sin(sol.solveval[i]);
            for(int j=0;j<3;++j){
                newnormals[i*3+j]   = y(npt*(j+1)+i) = coff_cos[i*3+j]*cosa+coff_sin[i*3+j]*sina;
            }
            MyUtility::normalize(newnormals.data()+i*3);
        }
    }else{
    }




    Set_RBFCoef(y);


    //sol.energy = arma::dot(a,M*a);
    cout<<"Solve_Hermite_PredictNormal_TangentConstraint_UnitNormal"<<endl;

    return 1;

}


int RBF_Core::Opt_Hermite_PredictNormal_TangentConstraint_UnitNormal(){

    SetInitnormal_Uninorm();

    BuildLFrame(initnormals_uninorm);


    sol.solveval.resize(npt);
    for(int i=0;i<npt;++i){
        sol.solveval[i] = 0.0;
    }
    //cout<<"smallvec: "<<smallvec<<endl;

    if(1){
        vector<double>upper(npt);
        vector<double>lower(npt);
        for(int i=0;i<npt;++i){
            upper[i] = 2 * my_PI;
            lower[i] = -2 * my_PI;
        }


        countopt = 0;

        Solver::nloptwrapper(lower,upper,optfunc_Hermite_Tangents,this,1e-7,3000,sol);
        cout<<"number of call: "<<countopt<<endl;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

    }

    newnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;

    if(1){
        for(int i=0;i<npt;++i){
            double cosa = cos(sol.solveval[i]), sina = sin(sol.solveval[i]);
            for(int j=0;j<3;++j){
                newnormals[i*3+j]   = y(npt*(j+1)+i) = coff_cos[i*3+j]*cosa+coff_sin[i*3+j]*sina;
            }
            MyUtility::normalize(newnormals.data()+i*3);
        }
    }else{

    }




    Set_RBFCoef(y);


    //sol.energy = arma::dot(a,M*a);
    cout<<"Opt_Hermite_PredictNormal_TangentConstraint_UnitNormal"<<endl;

    return 1;

}


void RBF_Core::SparseMatrixTest(){

    int nnpt = 500;
    int n = nnpt*4;
    int k = 5;
    auto t1 = Clock::now();
    arma::sp_mat A = arma::sprandu<arma::sp_mat>(n*n, n*k, 1e-4);//double(k)/nnpt
    arma::vec bb;
    cout << "init time: " << (setup_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9) << endl;
    t1 = Clock::now();
    arma::sp_mat AA(A.t()*A);
    bb.ones(n*n);
    bb=A.t()*bb;
    cout << "AAbb time: " << (setup_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9) << endl;
    cout<<"SparseMatrixTest "<<AA.n_nonzero<<endl;
    t1 = Clock::now();
    arma::vec X = arma::spsolve( AA,bb );

    cout << "Setup time: " << (setup_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9) << endl;

    exit(12321);
}


void RBF_Core::LocalEigenInit_PredictNormal_TangentConstraint(){


    cout<<"LocalEigenInit_PredictNormal"<<endl;
    arma::mat eigenmatrix;
    arma::vec eigenval;

    initnormals.resize(npt*3);

    vector<double>eva(3);
    for(int i=0;i<npt;++i){
        arma::mat aa(3,3);
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(i+j*npt,i+k*npt);

        arma::mat c2(1,3);
        arma::mat eye;
        auto p_tangent = tangents.data();
        {
            c2(0,0) = p_tangent[i*3];
            c2(0,1) = p_tangent[i*3+1];
            c2(0,2) = p_tangent[i*3+2];
        }
        eye.eye(3,3);
        arma::mat cct(c2*c2.t());

        arma::mat p = eye - c2.t()*arma::inv(cct)*c2;

        arma::mat AA = p*aa*p;
        AA = (AA.t()+AA)*0.5;
        eig_sym( eigenval, eigenmatrix, AA);

        vector<bool>isNullVec(3,false);
        for(int i=0;i<3;++i){
            double val = arma::norm(c2*eigenmatrix.col(i));
            if(val<1e-5)isNullVec[i]=true;
        }

        int smalleig = 0;
        for(int i=0;i<npt*3;++i)if(isNullVec[i]){smalleig = i;break;}

        for(int j=0;j<3;++j)eva[j] = eigenval(j);
        double stdv = eigenval(2) - eigenval(1);
        //cout<<stdv<<' ';
        for(int j=0;j<3;++j)initnormals[i*3+j] = stdv*eigenmatrix(j,smalleig);
    }
    cout<<endl;

    //newnormals = initnormals;

    NormalReOrientation(initnormals,false);


}



void RBF_Core::BuildLFrame(vector<double>&nors){

    double css[3];
    auto p_nor = nors.data();
    coff_cos.resize(npt*3);
    coff_sin.resize(npt*3);
    auto p_cos = coff_cos.data();
    auto p_sin = coff_sin.data();
    for(int i=0;i<npt;++i){

        MyUtility::cross(tangents.data()+i*3,p_nor+i*3,css);
        MyUtility::normalize(css);
        MyUtility::copyVec(p_nor+i*3,p_cos+i*3);
        MyUtility::copyVec(css,p_sin+i*3);
    }

}


int RBF_Core::Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm_Iterative_MST(){


	return 1;
}


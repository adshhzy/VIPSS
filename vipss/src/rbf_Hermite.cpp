#include "rbfcore.h"
#include "utility.h"
#include "Solver.h"
#include <armadillo>
#include <fstream>
#include <limits>
#include <unordered_map>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <queue>
#include "readers.h"
//#include "mymesh/UnionFind.h"
//#include "mymesh/tinyply.h"

typedef std::chrono::high_resolution_clock Clock;
double randomdouble() {return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);}
double randomdouble(double be,double ed) {return be + randomdouble()*(ed-be);	}

void RBF_Core::NormalRecification(double maxlen, vector<double>&nors){


    double maxlen_r = -1;
    auto p_vn = nors.data();
    int  np = nors.size()/3;
    if(1){
        for(int i=0;i<np;++i){
            maxlen_r = max(maxlen_r,MyUtility::normVec(p_vn+i*3));
        }

        cout<<"maxlen_r: "<<maxlen_r<<endl;
        double ratio = maxlen / maxlen_r;
        for(auto &a:nors)a*=ratio;
    }else{
        for(int i=0;i<np;++i){
            MyUtility::normalize(p_vn+i*3);
        }

    }




}

bool RBF_Core::Write_Hermite_NormalPrediction(string fname, int mode){


//    vector<uchar>labelcolor(npt*4);
//    vector<uint>f2v;
//    uchar red[] = {255,0,0, 255};
//    uchar green[] = {0,255,0, 255};
//    uchar blue[] = {0,0,255, 255};
//    for(int i=0;i<labels.size();++i){
//        uchar *pcolor;
//        if(labels[i]==0)pcolor = green;
//        else if(labels[i]==-1)pcolor = blue;
//        else if(labels[i]==1)pcolor = red;
//        for(int j=0;j<4;++j)labelcolor[i*4+j] = pcolor[j];
//    }
    //fname += mp_RBF_METHOD[curMethod];

//    for(int i=0;i<npt;++i){
//        uchar *pcolor = green;
//        for(int j=0;j<4;++j)labelcolor[i*4+j] = pcolor[j];
//    }

    vector<double>nors;
    if(mode ==0)nors=initnormals;
    else if(mode == 1)nors=newnormals;
    else if(mode == 2)nors = initnormals_uninorm;
    NormalRecification(1.,nors);

    //for(int i=0;i<npt;++i)if(randomdouble()<0.5)MyUtility::negVec(nors.data()+i*3);
    //cout<<pts.size()<<' '<<f2v.size()<<' '<<nors.size()<<' '<<labelcolor.size()<<endl;
    //writePLYFile(fname,pts,f2v,nors,labelcolor);

//    writeObjFile_vn(fname,pts,nors);
    writePLYFile_VN(fname,pts,nors);

    return 1;
}



void RBF_Core::Set_HermiteRBF(vector<double>&pts){

    cout<<"Set_HermiteRBF"<<endl;
    //for(auto a:pts)cout<<a<<' ';cout<<endl;
    isHermite = true;

    a.set_size(npt*4);
    M.set_size(npt*4,npt*4);
    double *p_pts = pts.data();
    for(int i=0;i<npt;++i){
        for(int j=i;j<npt;++j){
            M(i,j) = M(j,i) = Kernal_Function_2p(p_pts+i*3, p_pts+j*3);
        }
    }


    //if(User_Lamnbda!=0)for(int i=0;i<npt;++i)M(i,i) += User_Lamnbda;


    double G[3];
    for(int i=0;i<npt;++i){
        for(int j=0;j<npt;++j){

            Kernal_Gradient_Function_2p(p_pts+i*3, p_pts+j*3, G);
            //            int jind = j*3+npt;
            //            for(int k=0;k<3;++k)M(i,jind+k) = -G[k];
            //            for(int k=0;k<3;++k)M(jind+k,i) = G[k];

            for(int k=0;k<3;++k)M(i,npt+j+k*npt) = G[k];
            for(int k=0;k<3;++k)M(npt+j+k*npt,i) = G[k];

        }
    }

    double H[9];
    for(int i=0;i<npt;++i){
        for(int j=i;j<npt;++j){

            Kernal_Hessian_Function_2p(p_pts+i*3, p_pts+j*3, H);
            //            int iind = i*3+npt;
            //            int jind = j*3+npt;
            //            for(int k=0;k<3;++k)
            //                for(int l=0;l<3;++l)
            //                    M(jind+l,iind+k) = M(iind+k,jind+l) = -H[k*3+l];

            for(int k=0;k<3;++k)
                for(int l=0;l<3;++l)
                    M(npt+j+l*npt,npt+i+k*npt) = M(npt+i+k*npt,npt+j+l*npt) = -H[k*3+l];
        }
    }

    //cout<<std::setprecision(5)<<std::fixed<<M<<endl;

    bsize= 4;
    N.zeros(npt*4,4);
    b.set_size(4);

    for(int i=0;i<npt;++i){
        N(i,0) = 1;
        for(int j=0;j<3;++j)N(i,j+1) = pts[i*3+j];
    }
    for(int i=0;i<npt;++i){
        //        int ind = i*3+npt;
        //        for(int j=0;j<3;++j)N(ind+j,j+1) = 1;

        for(int j=0;j<3;++j)N(npt+i+j*npt,j+1) = -1;
    }

    //cout<<N<<endl;
    //arma::vec eigval = eig_sym( M ) ;
    //cout<<eigval.t()<<endl;


    if(!isnewformula){
        cout<<"start solve M: "<<endl;
        auto t1 = Clock::now();
        if(isinv)Minv = inv(M);
        else {
            arma::mat Eye;
            Eye.eye(npt*4,npt*4);
            Minv = solve(M,Eye);
        }
        cout<<"solved M: "<<(invM_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;

        t1 = Clock::now();
        if(isinv)bprey = inv_sympd(N.t() * Minv * N) * N.t() * Minv;
        else {
            arma::mat Eye2;
            Eye2.eye(bsize,bsize);
            bprey = solve(N.t() * Minv * N, Eye2) * N.t() * Minv;
        }
        cout<<"solved bprey "<<std::chrono::nanoseconds(Clock::now() - t1).count()/1e9<<endl;
    }else{



    }
}


double Gaussian_2p(const double *p1, const double *p2, double sigma){

    return exp(-MyUtility::vecSquareDist(p1,p2)/(2*sigma*sigma));
}



void RBF_Core::Set_Actual_User_LSCoef(double user_ls){

    User_Lamnbda = User_Lamnbda_inject = user_ls > 0 ?  user_ls : 0;

}

void RBF_Core::Set_Actual_Hermite_LSCoef(double hermite_ls){

    ls_coef = Hermite_ls_weight_inject = hermite_ls > 0?hermite_ls:0;
}

void RBF_Core::Set_SparsePara(double spa){
    sparse_para = spa;
}

void RBF_Core::Set_User_Lamnda_ToMatrix(double user_ls){


    {
        Set_Actual_User_LSCoef(user_ls);
        auto t1 = Clock::now();
        cout<<"setting K, HermiteApprox_Lamnda"<<endl;
        if(User_Lamnbda>0){
            arma::sp_mat eye;
            eye.eye(npt,npt);

            dI = inv(eye + User_Lamnbda*K00);
            saveK_finalH = K = K11 - (User_Lamnbda)*(K01.t()*dI*K01);

        }else saveK_finalH = K = K11;
        cout<<"solved: "<<(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
    }

    finalH = saveK_finalH;

}

void RBF_Core::Set_HermiteApprox_Lamnda(double hermite_ls){


    {
        Set_Actual_Hermite_LSCoef(hermite_ls);
        auto t1 = Clock::now();
        cout<<"setting K, HermiteApprox_Lamnda"<<endl;
        if(ls_coef>0){
            arma::sp_mat eye;
            eye.eye(npt,npt);

            if(ls_coef > 0){
                arma:: mat tmpdI = inv(eye + (ls_coef+User_Lamnbda)*K00);
                K = K11 - (ls_coef+User_Lamnbda)*(K01.t()*tmpdI*K01);
            }else{
                K = saveK_finalH;
            }
        }
        cout<<"solved: "<<(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;    
    }


}



void RBF_Core::Set_Hermite_PredictNormal(vector<double>&pts){



    Set_HermiteRBF(pts);

    auto t1 = Clock::now();
    cout<<"setting K"<<endl;


    if(!isnewformula){
        arma::mat D = N.t()*Minv;
        K = Minv - D.t()*inv(D*N)*D;
        K = K.submat( npt, npt, npt*4-1, npt*4-1 );
        finalH = saveK_finalH = K;

    }else{
        cout<<"using new formula"<<endl;
        bigM.zeros((npt+1)*4,(npt+1)*4);
        bigM.submat(0,0,npt*4-1,npt*4-1) = M;
        bigM.submat(0,npt*4,(npt)*4-1, (npt+1)*4-1) = N;
        bigM.submat(npt*4,0,(npt+1)*4-1, (npt)*4-1) = N.t();

        //for(int i=0;i<4;++i)bigM(i+(npt)*4,i+(npt)*4) = 1;

        auto t2 = Clock::now();
        bigMinv = inv(bigM);
        cout<<"bigMinv: "<<(setK_time= std::chrono::nanoseconds(Clock::now() - t2).count()/1e9)<<endl;
		bigM.clear();
        Minv = bigMinv.submat(0,0,npt*4-1,npt*4-1);
        Ninv = bigMinv.submat(0,npt*4,(npt)*4-1, (npt+1)*4-1);

        bigMinv.clear();
        //K = Minv - Ninv *(N.t()*Minv);
        K = Minv;
        K00 = K.submat(0,0,npt-1,npt-1);
        K01 = K.submat(0,npt,npt-1,npt*4-1);
        K11 = K.submat( npt, npt, npt*4-1, npt*4-1 );

        M.clear();N.clear();
        cout<<"K11: "<<K11.n_cols<<endl;


        //Set_Hermite_DesignedCurve();

        Set_User_Lamnda_ToMatrix(User_Lamnbda_inject);

		
//		arma::vec eigval, ny;
//		arma::mat eigvec;
//		ny = eig_sym( eigval, eigvec, K);
//		cout<<ny<<endl;

        cout<<"K: "<<K.n_cols<<endl;
    }




    //K = ( K.t() + K )/2;
    cout<<"solve K total: "<<(setK_time= std::chrono::nanoseconds(Clock::now() - t1).count()/1e9)<<endl;
    return;

}



void RBF_Core::SetInitnormal_Uninorm(){

    initnormals_uninorm = initnormals;
    for(int i=0;i<npt;++i)MyUtility::normalize(initnormals_uninorm.data()+i*3);

}

int RBF_Core::Solve_Hermite_PredictNormal_UnitNorm(){

    arma::vec eigval, ny;
    arma::mat eigvec;

    if(!isuse_sparse){
        ny = eig_sym( eigval, eigvec, K);
    }else{
//		cout<<"use sparse eigen"<<endl;
//        int k = 4;
//        do{
//            ny = eigs_sym( eigval, eigvec, sp_K, k, "sa" );
//            k+=4;
//        }while(ny(0)==0);
    }


    cout<<"eigval(0): "<<eigval(0)<<endl;

    int smalleig = 0;

    initnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt*3;++i)y(i+npt) = eigvec(i,smalleig);
    for(int i=0;i<npt;++i){
        initnormals[i*3]   = y(npt+i);
        initnormals[i*3+1] = y(npt+i+npt);
        initnormals[i*3+2] = y(npt+i+npt*2);
        //MyUtility::normalize(normals.data()+i*3);
    }


    SetInitnormal_Uninorm();
    cout<<"Solve_Hermite_PredictNormal_UnitNorm finish"<<endl;
    return 1;
}



/***************************************************************************************************/
/***************************************************************************************************/
double acc_time;

static int countopt = 0;
double optfunc_Hermite(const vector<double>&x, vector<double>&grad, void *fdata){

    auto t1 = Clock::now();
    RBF_Core *drbf = reinterpret_cast<RBF_Core*>(fdata);
    int n = drbf->npt;
    arma::vec arma_x(n*3);

    //(  sin(a)cos(b), sin(a)sin(b), cos(a)  )  a =>[0, pi], b => [-pi, pi];
    vector<double>sina_cosa_sinb_cosb(n * 4);
    for(int i=0;i<n;++i){
        int ind = i*4;
        sina_cosa_sinb_cosb[ind] = sin(x[i*2]);
        sina_cosa_sinb_cosb[ind+1] = cos(x[i*2]);
        sina_cosa_sinb_cosb[ind+2] = sin(x[i*2+1]);
        sina_cosa_sinb_cosb[ind+3] = cos(x[i*2+1]);
    }

    for(int i=0;i<n;++i){
        auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;
        //        int ind = i*3;
        //        arma_x(ind) = p_scsc[0] * p_scsc[3];
        //        arma_x(ind+1) = p_scsc[0] * p_scsc[2];
        //        arma_x(ind+2) = p_scsc[1];
        arma_x(i) = p_scsc[0] * p_scsc[3];
        arma_x(i+n) = p_scsc[0] * p_scsc[2];
        arma_x(i+n*2) = p_scsc[1];
    }

    arma::vec a2;
    //if(drbf->isuse_sparse)a2 = drbf->sp_H * arma_x;
    //else
    a2 = drbf->finalH * arma_x;


    if (!grad.empty()) {

        grad.resize(n*2);

        for(int i=0;i<n;++i){
            auto p_scsc = sina_cosa_sinb_cosb.data()+i*4;

            //            int ind = i*3;
            //            grad[i*2] = a2(ind) * p_scsc[1] * p_scsc[3] + a2(ind+1) * p_scsc[1] * p_scsc[2] - a2(ind+2) * p_scsc[0];
            //            grad[i*2+1] = -a2(ind) * p_scsc[0] * p_scsc[2] + a2(ind+1) * p_scsc[0] * p_scsc[3];

            grad[i*2] = a2(i) * p_scsc[1] * p_scsc[3] + a2(i+n) * p_scsc[1] * p_scsc[2] - a2(i+n*2) * p_scsc[0];
            grad[i*2+1] = -a2(i) * p_scsc[0] * p_scsc[2] + a2(i+n) * p_scsc[0] * p_scsc[3];

        }
    }

    double re = arma::dot( arma_x, a2 );
    countopt++;

    acc_time+=(std::chrono::nanoseconds(Clock::now() - t1).count()/1e9);

    //cout<<countopt++<<' '<<re<<endl;
    return re;

}



int RBF_Core::Opt_Hermite_PredictNormal_UnitNormal(){


    sol.solveval.resize(npt * 2);

    for(int i=0;i<npt;++i){
        double *veccc = initnormals.data()+i*3;
        {
            //MyUtility::normalize(veccc);
            sol.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
            sol.solveval[i*2 + 1] = atan2( veccc[1], veccc[0]   );
        }

    }
    //cout<<"smallvec: "<<smallvec<<endl;

    if(1){
        vector<double>upper(npt*2);
        vector<double>lower(npt*2);
        for(int i=0;i<npt;++i){
            upper[i*2] = 2 * my_PI;
            upper[i*2 + 1] = 2 * my_PI;

            lower[i*2] = -2 * my_PI;
            lower[i*2 + 1] = -2 * my_PI;
        }

        countopt = 0;
        acc_time = 0;

        //LocalIterativeSolver(sol,kk==0?normals:newnormals,300,1e-7);
        Solver::nloptwrapper(lower,upper,optfunc_Hermite,this,1e-7,3000,sol);
        cout<<"number of call: "<<countopt<<" t: "<<acc_time<<" ave: "<<acc_time/countopt<<endl;
        callfunc_time = acc_time;
        solve_time = sol.time;
        //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

    }
    newnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt;++i){

        double a = sol.solveval[i*2], b = sol.solveval[i*2+1];
        newnormals[i*3]   = y(npt+i) = sin(a) * cos(b);
        newnormals[i*3+1] = y(npt+i+npt) = sin(a) * sin(b);
        newnormals[i*3+2] = y(npt+i+npt*2) = cos(a);
        MyUtility::normalize(newnormals.data()+i*3);
    }

    Set_RBFCoef(y);

    //sol.energy = arma::dot(a,M*a);
    cout<<"Opt_Hermite_PredictNormal_UnitNormal"<<endl;
    return 1;
}

void RBF_Core::Set_RBFCoef(arma::vec &y){
    cout<<"Set_RBFCoef"<<endl;
    if(curMethod==HandCraft){
        cout<<"HandCraft, not RBF"<<endl;
        return;
    }
    if(!isnewformula){
        b = bprey * y;
        a = Minv * (y - N*b);
    }else{

        if(User_Lamnbda>0)y.subvec(0,npt-1) = -User_Lamnbda*dI*K01*y.subvec(npt,npt*4-1);

        a = Minv*y;
        b = Ninv.t()*y;

    }


}



int RBF_Core::Lamnbda_Search_GlobalEigen(){

    vector<double>lamnbda_list({0, 0.001, 0.01, 0.1, 1});
    //vector<double>lamnbda_list({  0.5,0.6,0.7,0.8,0.9,1,1.1,1.5,2,3});
    //lamnbda_list.clear();
    //for(double i=1.5;i<2.5;i+=0.1)lamnbda_list.push_back(i);
    //vector<double>lamnbda_list({0});
    vector<double>initen_list(lamnbda_list.size());
    vector<double>finalen_list(lamnbda_list.size());
    vector<vector<double>>init_normallist;
    vector<vector<double>>opt_normallist;

    lamnbda_list_sa = lamnbda_list;
    for(int i=0;i<lamnbda_list.size();++i){

        Set_HermiteApprox_Lamnda(lamnbda_list[i]);

        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm();
        }

        //Solve_Hermite_PredictNormal_UnitNorm();
        OptNormal(1);

        initen_list[i] = sol.init_energy;
        finalen_list[i] = sol.energy;

        init_normallist.emplace_back(initnormals);
        opt_normallist.emplace_back(newnormals);
    }

    lamnbdaGlobal_Be.emplace_back(initen_list);
    lamnbdaGlobal_Ed.emplace_back(finalen_list);

    cout<<std::setprecision(8);
    for(int i=0;i<initen_list.size();++i){
        cout<<lamnbda_list[i]<<": "<<initen_list[i]<<" -> "<<finalen_list[i]<<endl;
    }

    int minind = min_element(finalen_list.begin(),finalen_list.end()) - finalen_list.begin();
    cout<<"min energy: "<<endl;
    cout<<lamnbda_list[minind]<<": "<<initen_list[minind]<<" -> "<<finalen_list[minind]<<endl;


    initnormals = init_normallist[minind];
    SetInitnormal_Uninorm();
    newnormals = opt_normallist[minind];
	return 1;
}




void RBF_Core::Print_LamnbdaSearchTest(string fname){


    cout<<setprecision(7);
    cout<<"Print_LamnbdaSearchTest"<<endl;
    for(int i=0;i<lamnbda_list_sa.size();++i)cout<<lamnbda_list_sa[i]<<' ';cout<<endl;
    cout<<lamnbdaGlobal_Be.size()<<endl;
    for(int i=0;i<lamnbdaGlobal_Be.size();++i){
        for(int j=0;j<lamnbdaGlobal_Be[i].size();++j){
            cout<<lamnbdaGlobal_Be[i][j]<<"\t"<<lamnbdaGlobal_Ed[i][j]<<"\t";
        }
        cout<<gtBe[i]<<"\t"<<gtEd[i]<<endl;
    }

    ofstream fout(fname);
    fout<<setprecision(7);
    if(!fout.fail()){
        for(int i=0;i<lamnbda_list_sa.size();++i)fout<<lamnbda_list_sa[i]<<' ';fout<<endl;
        fout<<lamnbdaGlobal_Be.size()<<endl;
        for(int i=0;i<lamnbdaGlobal_Be.size();++i){
            for(int j=0;j<lamnbdaGlobal_Be[i].size();++j){
                fout<<lamnbdaGlobal_Be[i][j]<<"\t"<<lamnbdaGlobal_Ed[i][j]<<"\t";
            }
            fout<<gtBe[i]<<"\t"<<gtEd[i]<<endl;
        }
    }
    fout.close();

}



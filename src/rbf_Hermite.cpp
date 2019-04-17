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
#include <algorithm>
#include <queue>
#include "mymesh/readers.h"
#include "mymesh/UnionFind.h"
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


    vector<uchar>labelcolor(npt*4);
    vector<uint>f2v;
    uchar red[] = {255,0,0, 255};
    uchar green[] = {0,255,0, 255};
    uchar blue[] = {0,0,255, 255};
//    for(int i=0;i<labels.size();++i){
//        uchar *pcolor;
//        if(labels[i]==0)pcolor = green;
//        else if(labels[i]==-1)pcolor = blue;
//        else if(labels[i]==1)pcolor = red;
//        for(int j=0;j<4;++j)labelcolor[i*4+j] = pcolor[j];
//    }
    //fname += mp_RBF_METHOD[curMethod];

    for(int i=0;i<npt;++i){
        uchar *pcolor = green;
        for(int j=0;j<4;++j)labelcolor[i*4+j] = pcolor[j];
    }

    vector<double>nors;
    if(mode ==0)nors=initnormals;
    else if(mode == 1)nors=newnormals;
    else if(mode == 2)nors = initnormals_uninorm;
    NormalRecification(1.,nors);

    //for(int i=0;i<npt;++i)if(randomdouble()<0.5)MyUtility::negVec(nors.data()+i*3);

    //cout<<pts.size()<<' '<<f2v.size()<<' '<<nors.size()<<' '<<labelcolor.size()<<endl;
    writePLYFile(fname,pts,f2v,nors,labelcolor);

    return 1;
}

bool RBF_Core::Write_Hermite_MST(string fname){

    return writeObjFile_line(fname,pts,MST_edges);

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



int RBF_Core::Solve_HermiteRBF(vector<double>&vn){

    cout<<"Solve_HermiteRBF "<<npt<<' '<<vn.size()<<endl;
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;
    //for(int i=0;i<npt*3;++i)y(i+npt) = normals[i];
    for(int i=0;i<npt;++i)for(int j=0;j<3;++j)y(npt+i+j*npt) = -vn[i*3+j];

    //cout<<y<<endl;

    b = bprey * y;
    a = Minv * (y - N*b);
    //    cout<<a.t()<<endl;
    //    cout<<b.t()<<endl;
    sol.energy = arma::dot(a,M*a);
    cout<<"Solve_HermiteRBF"<<endl;
    return 1;
}


double Gaussian_2p(const double *p1, const double *p2, double sigma){

    return exp(-MyUtility::vecSquareDist(p1,p2)/(2*sigma*sigma));
}


void MinandAveDist(vector<double>&pts, double &mindist, double &avedist){

    int n = pts.size()/3;
    if(n==0)return;
    mindist = 1e30;
    avedist = 0;
    auto p_v = pts.data();
    for(int i=0;i<n;++i){
        for(int j=i+1;j<n;++j){

            double dist = MyUtility::_VerticesDistance(p_v+i*3,p_v+j*3);
            mindist = min(mindist, dist);
            avedist += dist;




        }
    }
    avedist/=n;



}


void RBF_Core::Set_Actual_User_LSCoef(double user_ls){

    User_Lamnbda = User_Lamnbda_inject = user_ls > 0 ?  user_ls : 0;

}

void RBF_Core::Set_Actual_Hermite_LSCoef(double hermite_ls){

    ls_coef = Hermite_ls_weight_inject = hermite_ls > 0?hermite_ls:0;
}

void RBF_Core::Set_SparsePara(double spa){

    sparse_para = spa;
    if(isuse_sparse){
        SparsifyKH(finalH, sp_H);
        sp_K = sp_H;
    }
}

void RBF_Core::SparsifyKH(arma::mat &KH, arma::sp_mat &sp_KH){

    if(isuse_sparse){
        arma::mat spKH(KH);
        double cutele;
        int spmethod = 3;
        if(spmethod == 0){
            cout<<"spmethod 0"<<endl;
            cutele = max(fabs(spKH.max()),fabs(spKH.min())) * sparse_para;
            for(int i=0;i<npt*3;++i)for(int j=0;j<npt*3;++j)if(fabs(spKH(i,j))<cutele){
                spKH(i,j) = 0;
            }
        }else if(spmethod == 1){
            cout<<"spmethod 1"<<endl;
            Coherence_M.set_size(npt,npt);
            cutele = max(fabs(spKH.max()),fabs(spKH.min())) * sparse_para;
            //cutele =  1e-2;
            for(int p1=0; p1<npt;++p1)for(int p2=p1+1; p2<npt;++p2){
                arma::mat aa(3,3);
                for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = spKH(p1+j*npt,p2+k*npt);
                double ll  = Coherence_M(p1,p2) = Coherence_M(p2,p1) =  arma::norm(aa);
                if(ll<cutele)for(int j=0;j<3;++j)for(int k=0;k<3;++k){spKH(p1+j*npt,p2+k*npt) = spKH(p2+j*npt,p1+k*npt) = 0;}
            }
        }else if(spmethod == 2){
            cout<<"spmethod 2"<<endl;
            Coherence_M.set_size(npt,npt);

            for(int p1=0; p1<npt;++p1)Coherence_M(p1,p1) = 0;
            for(int p1=0; p1<npt;++p1)for(int p2=p1+1; p2<npt;++p2){
                arma::mat aa(3,3);
                for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = spKH(p1+j*npt,p2+k*npt);
                Coherence_M(p1,p2) = Coherence_M(p2,p1) =  arma::norm(aa);
            }
            cutele = Coherence_M.max() * sparse_para;
            for(int p1=0; p1<npt;++p1)for(int p2=p1+1; p2<npt;++p2){
                if(Coherence_M(p1,p2)<cutele)for(int j=0;j<3;++j)for(int k=0;k<3;++k){spKH(p1+j*npt,p2+k*npt) = spKH(p2+j*npt,p1+k*npt) = 0;}
            }
        }else if(spmethod == 3){
            cout<<"spmethod 3"<<endl;
            Coherence_M.set_size(npt,npt);
            for(int p1=0; p1<npt;++p1)Coherence_M(p1,p1) = 0;
            for(int p1=0; p1<npt;++p1)for(int p2=p1+1; p2<npt;++p2){
                arma::mat aa(3,3);
                for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = spKH(p1+j*npt,p2+k*npt);
                Coherence_M(p1,p2) = Coherence_M(p2,p1) =  arma::norm(aa);
            }
            vector<double>rowmax(npt),colmax(npt);
            for(int p1=0; p1<npt;++p1){
                arma::rowvec vv_r = Coherence_M.row(p1);
                rowmax[p1] = vv_r.max();
                arma::vec vv_c = Coherence_M.col(p1);
                colmax[p1] = vv_c.max();
            }
            for(int p1=0; p1<npt;++p1)for(int p2=p1+1; p2<npt;++p2){
                cutele = min(rowmax[p1],colmax[p2]) * sparse_para;
                if(Coherence_M(p1,p2)<cutele)for(int j=0;j<3;++j)for(int k=0;k<3;++k){spKH(p1+j*npt,p2+k*npt) = spKH(p2+j*npt,p1+k*npt) = 0;}
            }

        }
       // spKH.for_each( [&cutele](arma::mat::elem_type& val) { if(fabs(val)<cutele)val=0; } );
        sp_KH = arma::sp_mat(spKH);
        int nsparse = sp_KH.n_nonzero;
        cout<<nsparse<<' '<<double(nsparse)/sp_KH.n_elem<<' '<<sp_KH.max()<<' '<<cutele<<endl;
    }
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

    SparsifyKH(finalH, sp_H);
    if(isuse_sparse)sp_K = sp_H;

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

    SparsifyKH(K,sp_K);

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
        SparsifyKH(K,sp_K);
        SparsifyKH(finalH, sp_H);
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


        Set_Hermite_DesignedCurve();

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


void RBF_Core::KMatrixAnalysis(){

    int pickpnt = 60;

    //cout<<"sum: "<<arma::sum(K.col(0))<<endl;


    vector<double>weights(npt);
    for(int i=0;i<npt;++i){
        if(i!=pickpnt){
            weights[i] = 0;
            for(int j=0;j<3;++j)for(int k=0;k<3;++k)weights[i]+=pow(K(pickpnt+j*npt,i+k*npt),2);
        }else{
            weights[i] = -1;
        }
    }
    auto sortweights = weights;
    sort(sortweights.begin(),sortweights.end());
    for(auto a:sortweights)cout<<a<<' ';cout<<endl;
    int pickmaxNei = max_element(weights.begin(),weights.end()) - weights.begin();

    arma::vec eigval, ny;
    arma::mat eigvec;


    vector<arma::mat>eigenmatrix(npt);
    vector<arma::vec>eigenval(npt);
    for(int i=0;i<npt;++i){
        arma::mat aa(3,3);
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(i+j*npt,i+k*npt);
        eig_sym( eigenval[i], eigenmatrix[i], aa);
        //cout<<eigenval[i].t();
    }


    arma::mat aa(3,3);
    for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(pickpnt+j*npt,pickmaxNei+k*npt);



    cout<<eigenmatrix[pickpnt]<<endl;
    cout<<eigenval[pickpnt].t()<<endl;
    cout<<eigenmatrix[pickmaxNei]<<endl;
    cout<<eigenval[pickmaxNei].t()<<endl;
    cout<<normalise( aa*eigenmatrix[pickmaxNei] )<<endl;
    cout<<aa<<endl;
    //cout<<aa*aa.t()<<endl;


    exit(0);

}

class Triple{
public:
    double val;
    int i,j;
    Triple();
    Triple(double val, int i, int j):val(val),i(i),j(j){}
};
bool operator <(const Triple& l, const Triple& r){return l.val < r.val;}
bool operator >(const Triple& l, const Triple& r){return l.val > r.val;}


void BuildMST(vector<uint>&MST_edges, vector<vector<int>>&MST_v2v, vector<double>&MST_edgesWeight, vector<Triple>&eweights, int n_v){


    MST_edges.clear();
    MST_v2v.clear();
    MST_edgesWeight.clear();
    MST_v2v.resize(n_v);

    sort(eweights.begin(),eweights.end(),greater<Triple>());

    int n_neg = 0;
    vector<int>ele(n_v);
    for(int i=0;i<n_v;++i)ele[i]=i;
    UnionFind unifind;
    unifind.SetElements(ele);
    double lastval;

    //for(auto &a:eweights)cout<<a.val<<' ';cout<<endl;
    for(auto &a:eweights){
        //cout<<a.val<<' ';
        if(!unifind.Is_Connected(a.i,a.j)){
            lastval=a.val;
            unifind.Union(a.i,a.j);
            MST_edges.push_back(a.i);MST_edges.push_back(a.j);
            MST_v2v[a.i].push_back(a.j);MST_v2v[a.j].push_back(a.i);
            MST_edgesWeight.push_back(a.val);
            ++n_neg;
        }
    }
    //cout<<"lastval: "<<lastval<<endl;
    //cout<<"n_neg: "<<n_neg<<endl;
    assert(unifind.GetNumOfComponents()==1);


}
void MST_Orientate(vector<bool>&isflip,vector<vector<double>>&edgelength, vector<vector<int>>&MST_v2v, int npt, int startpnt){

    int nflip = 0;
    vector<bool>isvisited(npt,false);
    isflip.clear();
    isflip.resize(npt,false);
    queue<pair<int,int>>q;
    for(auto a:MST_v2v[startpnt])q.emplace(a,startpnt);
    isvisited[startpnt]=true;
    while(!q.empty()){
        pair<int,int> aa = q.front();
        int curpnt = aa.first;
        q.pop();
        if(isvisited[curpnt])continue;
        isvisited[curpnt] = true;
        if((edgelength[curpnt][aa.second]<0 && isflip[aa.second]) || (edgelength[curpnt][aa.second]>0 && !isflip[aa.second])){
            isflip[curpnt] = true;++nflip;
        }
        for(auto a:MST_v2v[curpnt])q.emplace(a,curpnt);
    }

    //cout<<"nflip: "<<nflip<<endl;


}


void RBF_Core::NormalReOrientation(vector<double>&nors, bool isnormalized){


    cout<<"NormalReOrientation"<<endl;
    vector<double>nors_strength(npt);
    //arma::mat edgelength(npt,npt);
    vector<vector<double>>edgelength(npt,vector<double>(npt));

    vector<arma::vec>arma_normals(npt,arma::vec(3));
    auto p_nor = nors.data();


    for(int p1=0; p1<npt;++p1){
        for(int j=0;j<3;++j)arma_normals[p1](j) = p_nor[p1*3+j];
        nors_strength[p1] = arma::norm(arma_normals[p1] );
        arma_normals[p1] = arma::normalise( arma_normals[p1] );
    }

    vector<Triple>eabslength;
    int n_neg = 0;
    for(int p1=0; p1<npt;++p1)for(int p2=p1+1; p2<npt;++p2){
        arma::mat aa(3,3);
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(p1+j*npt,p2+k*npt);
        edgelength[p1][p2] = edgelength[p2][p1]  = arma::dot(arma_normals[p1], aa*arma_normals[p2]);
        if(edgelength[p1][p2]<0)n_neg++;
        if(  fabs(edgelength[p1][p2]) > 1e-4 )eabslength.emplace_back(fabs(edgelength[p1][p2]),p1,p2);
    }

    vector<vector<int>>MST_v2v;
    BuildMST(MST_edges, MST_v2v,MST_edgesWeight, eabslength, npt);



    vector<bool>isflip(npt,false);
    int startpnt = max_element(nors_strength.begin(),nors_strength.end())-nors_strength.begin();
    MST_Orientate(isflip,edgelength,MST_v2v,npt,startpnt);

    for(int p1=0; p1<npt;++p1)if(isflip[p1]){
        arma_normals[p1] = -arma_normals[p1];
    }

    for(int p1=0; p1<npt;++p1){
        for(int j=0;j<3;++j)p_nor[p1*3+j] = arma_normals[p1](j);
    }
    if(!isnormalized){
        for(int p1=0; p1<npt;++p1){
            for(int j=0;j<3;++j)p_nor[p1*3+j] = nors_strength[p1] * p_nor[p1*3+j];
        }
    }

}

void RBF_Core::NormalReOrientation_WithGT(vector<double> &nors, bool isnormalized){

    for(int i=0;i<npt;++i){
        auto p_gt = normals.data()+i*3;
        auto p_nor = nors.data()+i*3;

        if(MyUtility::dot(p_gt,p_nor)<0)for(int j=0;j<3;++j)p_nor[j] = -p_nor[j];
    }

    if(isnormalized)for(int i=0;i<npt;++i){
        MyUtility::normalize(nors.data()+i*3);
    }

}

void RBF_Core::BuildCoherentGraph(){

    cout<<"BuildCoherentGraph"<<endl;
    vector<Triple>eabslength;
    Coherence_M.set_size(npt,npt);


    for(int p1=0; p1<npt;++p1)for(int p2=p1; p2<npt;++p2){
        arma::mat aa(3,3);
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(p1+j*npt,p2+k*npt);
        double ll  = Coherence_M(p1,p2) = Coherence_M(p2,p1) =  arma::norm(aa);
        if(p2!=p1)
            if(  fabs(ll) > 1e-4 )
                eabslength.emplace_back(fabs(ll),p1,p2);
    }

    sort(eabslength.begin(),eabslength.end(),greater<Triple>());

    Coherence_Graph_Weight.clear();
    Coherence_Graph.clear();
    for(auto &a:eabslength){
        Coherence_Graph_Weight.push_back(a.val);
        Coherence_Graph.push_back(a.i);
        Coherence_Graph.push_back(a.j);
    }
}


double stdvar(vector<double>&v){

    int nv = v.size();
    double ave = 0;
    for(auto a:v)ave+=a;
    ave/=nv;

    double re = 0;
    for(auto a:v)re+= pow(ave-a,2);
    return sqrt(re/nv);


}


void RBF_Core::LocalEigenInit_PredictNormal(){


    cout<<"LocalEigenInit_PredictNormal"<<endl;
    arma::mat eigenmatrix;
    arma::vec eigenval;

    initnormals.resize(npt*3);

    vector<double>eva(3);
    for(int i=0;i<npt;++i){
        arma::mat aa(3,3);
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(i+j*npt,i+k*npt);
        eig_sym( eigenval, eigenmatrix, aa);
        for(int j=0;j<3;++j)eva[j] = eigenval(j);
        //double stdv = stdvar(eva);
        double stdv = eva[1] - eva[0];
        //cout<<stdv<<' ';
        for(int j=0;j<3;++j)initnormals[i*3+j] = stdv*eigenmatrix(j,0);
    }

    NormalReOrientation(initnormals,false);
    SetInitnormal_Uninorm();

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
		cout<<"use sparse eigen"<<endl;
        int k = 4;
        do{
            ny = eigs_sym( eigval, eigvec, sp_K, k, "sa" );
            k+=4;
        }while(ny(0)==0);
    }




    cout<<"eigval(0): "<<eigval(0)<<endl;

    int smalleig = 0;
    //    for(int i=0;i<eigval.size();++i)if(eigval(i)>0){
    //        smalleig = i;
    //        break;
    //    }
    //    if(smalleig!=0)if( fabs(eigval(smalleig-1)) < fabs(eigval(smalleig)) )smalleig-=1;

    //smalleig = 9;

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

    //b = bprey * y;
    //a = Minv * (y - N*b);

    //sol.energy = arma::dot(a,M*a);

    cout<<"Solve_Hermite_PredictNormal_UnitNorm finish"<<endl;
    return 1;
}


void WritePointPLY(string fname, vector<double>&out_Vs, vector<double>&normals, vector<int>&labels);
void RBF_Core::WriteSeletivePLY(string fname, vector<double>&allnormals, vector<int>&pickInd){
    vector<double>npts;
    vector<double>nnor;
    vector<int>labels(pickInd.size(),0);
    for(auto ind:pickInd){
        for(int j=0;j<3;++j)npts.push_back(pts[ind*3+j]);
        for(int j=0;j<3;++j)nnor.push_back(allnormals[ind*3+j]);
    }
    WritePointPLY(fname,npts,nnor,labels);
}

void RBF_Core::ClearSavedIterBatchInit(){

    InitBatch_Vs.clear();
    InitBatch_VNs.clear();
    InitBatch_VUNs.clear();
    InitBatch_VInds.clear();
    InitBatch_Coherence_Graph.clear();
    InitBatch_Coherence_Graph_Weight.clear();
}


void RBF_Core::SaveIterBatchInit(vector<double>&allnormals, vector<double>&allnormals_uninorm, vector<int>&pickInd){
    int nnInd = InitBatch_Vs.size();
    InitBatch_Vs.resize(nnInd+1);
    InitBatch_VNs.resize(nnInd+1);
    InitBatch_VUNs.resize(nnInd+1);
    InitBatch_Coherence_Graph.resize(nnInd+1);
    InitBatch_Coherence_Graph_Weight.resize(nnInd+1);

    InitBatch_VInds.emplace_back(pickInd);
    vector<double>&npts = InitBatch_Vs[nnInd];
    vector<double>&nor = InitBatch_VNs[nnInd];
    vector<double>&nnor = InitBatch_VUNs[nnInd];
    vector<uint>&gg = InitBatch_Coherence_Graph[nnInd];
    vector<double>&ggw = InitBatch_Coherence_Graph_Weight[nnInd];

    vector<bool>ispick(npt,false);
    for(auto ind:pickInd){
        for(int j=0;j<3;++j)npts.push_back(pts[ind*3+j]);
        for(int j=0;j<3;++j)nor.push_back(allnormals[ind*3+j]);
        for(int j=0;j<3;++j)nnor.push_back(allnormals_uninorm[ind*3+j]);
        ispick[ind]=true;
    }

    vector<int>sortInd(pickInd);
    sort(sortInd.begin(),sortInd.end());
    vector<int>invertInd(npt,-1);
    for(int p1=0; p1<pickInd.size();++p1)invertInd[pickInd[p1]] = p1;
    vector<Triple>eeT;
    for(int p1=0; p1<sortInd.size();++p1)for(int p2=p1+1; p2<sortInd.size();++p2){
        int i=sortInd[p1],j=sortInd[p2];
        eeT.emplace_back(Coherence_M(i,j),invertInd[i],invertInd[j]);
    }

    sort(eeT.begin(),eeT.end(),greater<Triple>());

    gg.clear();ggw.clear();
    gg.reserve(eeT.size()*2);ggw.reserve(eeT.size());
    for(auto &a:eeT){
        ggw.push_back(a.val);
        gg.push_back(a.i);
        gg.push_back(a.j);
    }

}

int RBF_Core::Solve_Hermite_PredictNormal_UnitNorm_Iterative(){

    vector<int>labels(npt,0);
    double thres_coef = 0.5;

    arma::vec eigval, ny;
    arma::mat eigvec;

    int iter = 0;
    int maxiter = 20;
    vector<int>activeInd(npt);
    vector<int>fixedInd;
    for(int i=0;i<npt;++i)activeInd[i] = i;

    arma::mat subactiveMat;
    arma::vec fixVec,addVec;

    initnormals.clear();
    initnormals.resize(npt*3,0);
    vector<double>norm_normals(npt*3,0);
    ClearSavedIterBatchInit();
    while(true){

        ++iter;
        int n_Active = activeInd.size();
        int n_Fixed = fixedInd.size();
        bool prebreak = (iter>=maxiter) /*|| (n_Active < npt * 0.05)*/;

        arma::uvec row_ind(n_Active*3);
        for(int i=0;i<n_Active;++i){
            for(int j=0;j<3;++j)row_ind(i*3+j) = activeInd[i]+j*npt;
        }
        //cout<<row_ind.t()<<endl;

        subactiveMat = K.submat(row_ind,row_ind);
        //cout<<K.n_rows<<" "<<K.n_cols<<endl;

        ny = eig_sym( eigval, eigvec, subactiveMat);
        arma::vec firstEigvec(eigvec.col(0));
        //cout<<firstEigvec.n_elem<<endl;

        double thres = sqrt(1./npt) * thres_coef ;
        vector<int>newactive,addfixInd;
        //cout<<"thres: "<<thres<<endl;
        for(int i=0;i<n_Active;++i){
            auto p_nor = initnormals.data()+activeInd[i]*3;
            auto p_nor_n = norm_normals.data()+activeInd[i]*3;
            for(int j=0;j<3;++j)p_nor_n[j] = p_nor[j] = firstEigvec(i*3+j);
            //cout<<i<<endl;
            if(!prebreak && MyUtility::normVec(p_nor)<thres){
                newactive.push_back(activeInd[i]);
                for(int j=0;j<3;++j)p_nor_n[j] = 0;
            }else{
                //cout<<MyUtility::normVec(p_nor)<<' ';
                addfixInd.push_back(activeInd[i]);
                fixedInd.push_back(activeInd[i]);
                MyUtility::normalize(p_nor_n);
            }
        }



        cout<<iter<<" addfix: "<<addfixInd.size()<<' '<<fixedInd.size()<<endl;
        if(n_Fixed>0){
            int n_Addfix = fixedInd.size()-n_Fixed;
            arma::uvec addfix_row_ind(n_Addfix*3);
            addVec.set_size(n_Addfix*3);
            for(int i=0;i<n_Addfix;++i){
                for(int j=0;j<3;++j)addVec(i*3+j) = norm_normals[addfixInd[i]*3+j];
                for(int j=0;j<3;++j)addfix_row_ind(i*3+j) = addfixInd[i]+j*npt;
            }

            arma::uvec fix_row_ind(n_Fixed*3);
            fixVec.set_size(n_Fixed*3);
            for(int i=0;i<n_Fixed;++i){
                for(int j=0;j<3;++j)fix_row_ind(i*3+j) = fixedInd[i]+j*npt;
                for(int j=0;j<3;++j)fixVec(i*3+j) = norm_normals[fixedInd[i]*3+j];
            }


            arma::mat subInteractiveMat = K.submat(addfix_row_ind,fix_row_ind);
            double inter_en = arma::dot(addVec,subInteractiveMat*fixVec);

            cout<<"inter_en: "<<inter_en<<endl;
            if(inter_en>0){
                cout<<"flip"<<endl;
                for(int i=0;i<n_Addfix;++i){
                    int ind = addfixInd[i]*3;
                    for(int j=0;j<3;++j)initnormals[ind+j] = -initnormals[ind+j];
                    for(int j=0;j<3;++j)norm_normals[ind+j] = -norm_normals[ind+j];
                }
            }
        }
        //        WritePointPLY(string("/Users/Research/Geometry/RBF/data/iter/nn_")+to_string(iter),pts,norm_normals,labels);
        //        WritePointPLY(string("/Users/Research/Geometry/RBF/data/iter/un_")+to_string(iter),pts,newnormals,labels);
        WriteSeletivePLY(string("/Users/Research/Geometry/RBF/data/iter/nn_")+to_string(iter),norm_normals,activeInd);
        WriteSeletivePLY(string("/Users/Research/Geometry/RBF/data/iter/un_")+to_string(iter),initnormals,activeInd);

        SaveIterBatchInit(initnormals,norm_normals,activeInd);
        activeInd = newactive;
        if(activeInd.size()==0)break;
        if(prebreak)break;

    }

    initnormals_uninorm = initnormals = norm_normals;

    cout<<"Solve_Hermite_PredictNormal_UnitNorm_Iterative"<<endl;
    return 1;



}

int ExhaustedCombination(vector<bool>&isflip, vector<vector<double>>&en_mat, int n_v){
    isflip.clear();
    isflip.resize(n_v,false);
    if(n_v<2)return 1;
    cout<<"ExhaustedCombination"<<endl;

    arma::mat arma_enmat(n_v,n_v);
    for(int i=0;i<n_v;++i)for(int j=0;j<n_v;++j)arma_enmat(i,j) = en_mat[i][j];
    for(int i=0;i<n_v;++i)arma_enmat(i,i) = 0;

    double min_en = DBL_MAX;
    int min_ind = -1;
    arma::vec fvec;
    fvec.ones(n_v);
    int npow = pow(2,n_v-1);
    for(int i=0;i<npow;++i){
        //cout<<i<<' ';
        int maker = 1;
        for(int j=0;j<n_v-1 && maker<=i+1;++j){
            cout<<j<<' ';
            fvec(j+1) = i & maker ? -1:1;
            maker<<=1;
        }
        //cout<<endl;
        //if(i==63)cout<<fvec.t()<<endl;
        double cur_en = arma::dot(fvec,arma_enmat*fvec);
        if(cur_en<min_en){
            min_ind = i;min_en = cur_en;
        }
    }

    int maker = 1;
    for(int j=0;maker<=min_ind+1;++j){
        isflip[j+1] = min_ind & maker ? true : false;
        maker<<=1;
    }
	return 1;
}


int RBF_Core::Solve_Hermite_PredictNormal_UnitNorm_Iterative_MST(){

    vector<int>labels(npt,0);
    double thres_coef = 0.9;

    arma::vec eigval, ny;
    arma::mat eigvec;

    int iter = 0;
    int maxiter = 100;
    vector<int>activeInd(npt);
    vector<int>fixedInd;
    for(int i=0;i<npt;++i)activeInd[i] = i;

    arma::mat subactiveMat;
    arma::vec fixVec,addVec;

    initnormals.clear();
    initnormals.resize(npt*3,0);
    vector<double>norm_normals(npt*3,0);
    vector<vector<int>>hyper_patch;

    ClearSavedIterBatchInit();
    while(true){

        ++iter;
        int n_Active = activeInd.size();
        int n_Fixed = fixedInd.size();
        bool prebreak = (iter>=maxiter) /*|| (n_Active < npt * 0.05)*/;

        arma::uvec row_ind(n_Active*3);
        for(int i=0;i<n_Active;++i){
            for(int j=0;j<3;++j)row_ind(i*3+j) = activeInd[i]+j*npt;
        }
        //cout<<row_ind.t()<<endl;

        subactiveMat = K.submat(row_ind,row_ind);
        //cout<<K.n_rows<<" "<<K.n_cols<<endl;

        ny = eig_sym( eigval, eigvec, subactiveMat);
        arma::vec firstEigvec(eigvec.col(0));
        //cout<<firstEigvec.n_elem<<endl;

        double thres = sqrt(1./npt) * thres_coef ;
        //double thres = sqrt(1./n_Active) * thres_coef ;
        hyper_patch.resize(hyper_patch.size()+1);
        vector<int>newactive;
        vector<int>&addfixInd = hyper_patch[hyper_patch.size()-1];
        //cout<<"thres: "<<thres<<endl;
        for(int i=0;i<n_Active;++i){
            auto p_nor = initnormals.data()+activeInd[i]*3;
            auto p_nor_n = norm_normals.data()+activeInd[i]*3;
            for(int j=0;j<3;++j)p_nor_n[j] = p_nor[j] = firstEigvec(i*3+j);
            //cout<<i<<endl;
            if(!prebreak && MyUtility::normVec(p_nor)<thres){
                newactive.push_back(activeInd[i]);
                for(int j=0;j<3;++j)p_nor_n[j]  = 0;
            }else{
                //cout<<MyUtility::normVec(p_nor)<<' ';
                addfixInd.push_back(activeInd[i]);
                fixedInd.push_back(activeInd[i]);
                MyUtility::normalize(p_nor_n);
            }
        }

        cout<<iter<<" addfix: "<<addfixInd.size()<<' '<<fixedInd.size()<<endl;
        //        WritePointPLY(string("/Users/Research/Geometry/RBF/data/iter/nn_")+to_string(iter),pts,norm_normals,labels);
        //        WritePointPLY(string("/Users/Research/Geometry/RBF/data/iter/un_")+to_string(iter),pts,newnormals,labels);
        //WriteSeletivePLY(string("/Users/Research/Geometry/RBF/data/iter/nn_")+to_string(iter),norm_normals,activeInd);
        //WriteSeletivePLY(string("/Users/Research/Geometry/RBF/data/iter/un_")+to_string(iter),initnormals,activeInd);
        SaveIterBatchInit(initnormals,norm_normals,activeInd);
        activeInd = newactive;
        if(activeInd.size()==0)break;
        if(prebreak)break;
    }


    vector<uint>MST_edges;
    vector<vector<int>>MST_v2v;
    vector<Triple>eweights;
    int n_hyper = hyper_patch.size();

    vector<arma::vec>hyperb_vec(hyper_patch.size());
    vector<arma::uvec>hyperb_rowind(hyper_patch.size());
    for(int i=0;i<hyper_patch.size();++i){
        int b_nv = hyper_patch[i].size();
        auto &hyperVec = hyperb_vec[i];
        auto &row_ind = hyperb_rowind[i];
        auto &hyper_ind = hyper_patch[i];
        hyperVec.set_size(b_nv*3);
        row_ind.set_size(b_nv*3);

        for(int k=0;k<b_nv;++k){
            for(int j=0;j<3;++j)hyperVec(k*3+j) = norm_normals[hyper_ind[k]*3+j];
            for(int j=0;j<3;++j)row_ind(k*3+j) = hyper_ind[k]+j*npt;
        }
    }
    cout<<"adsada"<<endl;

    vector<vector<double>>interactiveEn_mat(n_hyper,vector<double>(n_hyper));
    for(int i=0;i<hyper_patch.size();++i){
        for(int j=i+1;j<hyper_patch.size();++j){
            arma::mat subInteractiveMat = K.submat(hyperb_rowind[i],hyperb_rowind[j]);
            //cout<<"aaa"<<endl;
            interactiveEn_mat[j][i] = interactiveEn_mat[i][j] = arma::dot(hyperb_vec[i],subInteractiveMat*hyperb_vec[j]);
            //cout<<i<<" "<<j<<" "<<interactiveEn_mat[i][j]<<endl;
            eweights.emplace_back(fabs(interactiveEn_mat[i][j]),i,j);
        }
    }


    vector<bool>isflip(n_hyper,false);

    if(1){
        cout<<"BuildMST"<<endl;
        BuildMST(MST_edges, MST_v2v,MST_edgesWeight, eweights, n_hyper);
        MST_Orientate(isflip,interactiveEn_mat,MST_v2v,n_hyper,0);

    }else{
        ExhaustedCombination(isflip,interactiveEn_mat,n_hyper);
    }


    for(int i=0;i<n_hyper;++i)if(isflip[i]){
        int b_nv = hyper_patch[i].size();
        auto &hyper_ind = hyper_patch[i];
        for(int k=0;k<b_nv;++k){
            for(int j=0;j<3;++j)norm_normals[hyper_ind[k]*3+j] = -norm_normals[hyper_ind[k]*3+j];
        }
    }


    initnormals_uninorm = initnormals = norm_normals;
    this->MST_edges.clear();
    this->MST_edgesWeight.clear();
    cout<<"Solve_Hermite_PredictNormal_UnitNorm_Iterative"<<endl;
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
    if(drbf->isuse_sparse)a2 = drbf->sp_H * arma_x;
    else a2 = drbf->finalH * arma_x;


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


int RBF_Core::Solve_Hermite_PredictNormal_UnitNormal(){

    cout<<"ads"<<endl;
    sol.solveval.resize(npt * 2);

    //Solve_Hermite_PredictNormal_UnitNorm();
    //NormalReOrientation();


    for(int kk=2;kk<3;++kk){
        auto t1 = Clock::now();
        if(kk==1){
            //Solve_Hermite_PredictNormal_UnitNorm();
            //Solve_Hermite_PredictNormal_UnitNorm_Iterative();
            Solve_Hermite_PredictNormal_UnitNorm_Iterative_MST();
            //NormalReOrientation();
            //initnormals = newnormals;
        }else if(kk==2){
            LocalEigenInit_PredictNormal();

            //for(int i=0;i<npt;++i)MyUtility::normalize(newnormals.data()+i*3);
        }else if(kk==0){
            initnormals = normals;
            cout<<"initnormals: "<<normals.size()<<endl;
        }
        //
        init_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9;

        for(int i=0;i<npt;++i){
            double veccc[3];
            for(int j=0;j<3;++j)veccc[j] = initnormals[i*3+j];
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
    }

    //exit(-1);

    newnormals.resize(npt*3);
    arma::vec y(npt*4);
    for(int i=0;i<npt;++i)y(i) = 0;
    for(int i=0;i<npt;++i){

        double a = sol.solveval[i*2], b = sol.solveval[i*2+1];
        //        int ind = i*3 + npt;
        //        y(ind) = sin(a) * cos(b);
        //        y(ind+1) = sin(a) * sin(b);
        //        y(ind+2) = cos(a);
        newnormals[i*3]   = y(npt+i) = sin(a) * cos(b);
        newnormals[i*3+1] = y(npt+i+npt) = sin(a) * sin(b);
        newnormals[i*3+2] = y(npt+i+npt*2) = cos(a);
        MyUtility::normalize(newnormals.data()+i*3);
    }




    Set_RBFCoef(y);

    //sol.energy = arma::dot(a,M*a);
    cout<<"Solve_Hermite_PredictNormal_UnitNormal"<<endl;
	
	return 1;
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



int RBF_Core::Solve_Hermite_PredictNormal_UnitNormal_InitializationTest(){

    sol.solveval.resize(npt * 2);

    int global_eigenIndex  = 1;
    int local_eigenIndex  = 2;
    int gt_Index = 0;



    for(int kk=0;kk<2;++kk){

        if(kk==global_eigenIndex){
            Solve_Hermite_PredictNormal_UnitNorm();
            NormalReOrientation_WithGT(initnormals,false);
        }else if(kk==local_eigenIndex){
            LocalEigenInit_PredictNormal();
        }
        double veccc[3];
        for(int i=0;i<npt;++i){
            if(kk!=gt_Index)for(int j=0;j<3;++j)veccc[j] = newnormals[i*3+j];
            else for(int j=0;j<3;++j)veccc[j] = normals[i*3+j];

            {
                //MyUtility::normalize(veccc);
                sol.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
                sol.solveval[i*2 + 1] = atan2( veccc[1], veccc[0]   );
            }
        }

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

            Solver::nloptwrapper(lower,upper,optfunc_Hermite,this,1e-7,3000,sol);
            cout<<"number of call: "<<countopt<<endl;
            //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

        }

        if(kk==global_eigenIndex){
            eigenBe.push_back(sol.init_energy);
            eigenEd.push_back(sol.energy);
        }else if(kk==local_eigenIndex){
            local_eigenBe.push_back(sol.init_energy);
            local_eigenEd.push_back(sol.energy);
        }else if(kk==gt_Index){
            gtBe.push_back(sol.init_energy);
            gtEd.push_back(sol.energy);
        }

        continue;

    }
    //sol.energy = arma::dot(a,M*a);
    cout<<"Solve_Hermite_PredictNormal_UnitNormal"<<endl;
    return 1;
}


int RBF_Core::Solve_Hermite_PredictNormal_UnitNormal_SolverTest(){



    int global_solver  = 1;
    int local_solver  = 0;



    LocalEigenInit_PredictNormal();
    for(int kk=0;kk<2;++kk){

        if(kk==global_solver){
            sol.solveval.resize(npt * 2);
            double veccc[3];
            for(int i=0;i<npt;++i){
                for(int j=0;j<3;++j)veccc[j] = newnormals[i*3+j];
                {
                    //MyUtility::normalize(veccc);
                    sol.solveval[i*2] = atan2(sqrt(veccc[0]*veccc[0]+veccc[1]*veccc[1]),veccc[2] );
                    sol.solveval[i*2 + 1] = atan2( veccc[1], veccc[0]   );
                }
            }

            {
                vector<double>upper(npt*2);
                vector<double>lower(npt*2);
                for(int i=0;i<npt;++i){
                    upper[i*2] = 2 * my_PI;
                    upper[i*2 + 1] = 2 * my_PI;

                    lower[i*2] = -2 * my_PI;
                    lower[i*2 + 1] = -2 * my_PI;
                }

                countopt = 0;

                Solver::nloptwrapper(lower,upper,optfunc_Hermite,this,1e-7,3000,sol);
                cout<<"number of call: "<<countopt<<endl;
                //for(int i=0;i<npt;++i)cout<< sol.solveval[i]<<' ';cout<<endl;

            }
        }else if(kk == local_solver){
            LocalIterativeSolver(sol,newnormals,300,1e-7);
        }

        if(kk==local_solver){
            local_eigenBe.push_back(sol.init_energy);
            local_eigenEd.push_back(sol.energy);
        }else if(kk==global_solver){
            gtBe.push_back(sol.init_energy);
            gtEd.push_back(sol.energy);
        }

        continue;

    }
    //sol.energy = arma::dot(a,M*a);
    cout<<"Solve_Hermite_PredictNormal_UnitNormal"<<endl;
    return 1;
}


void RBF_Core::ClearInitializationTest(){

    lamnbdaGlobal_Be.clear();
    lamnbdaGlobal_Ed.clear();
    local_eigenBe.clear();
    local_eigenEd.clear();
    eigenEd.clear();
    eigenBe.clear();
    eigenEd.clear();
    gtBe.clear();
    gtEd.clear();
}

void RBF_Core::Print_InitializationTest(string fname){

    bool isglobalEigen = true;
    cout<<setprecision(5);
    cout<<"Eigen \t\t\t Gt"<<endl;
    if(isglobalEigen){ for(int i=0;i<gtBe.size();++i){
            cout<<eigenBe[i]<<" -> "<<eigenEd[i]<<"\t"<<gtBe[i]<<" -> "<<gtEd[i]<<endl;
        }
    }else{
        for(int i=0;i<gtBe.size();++i){
            cout<<local_eigenBe[i]<<" -> "<<local_eigenEd[i]<<"\t"<<gtBe[i]<<" -> "<<gtEd[i]<<endl;
        }
    }

    //    int nclose=0;
    //    for(int i=0;i<gtBe.size();++i){
    //        if(fabs(gtEd[i]-local_eigenEd[i]) < 10)nclose++;
    //    }
    //cout<<"nclose: "<<nclose<<endl;

    ofstream fout(fname);
    fout<<setprecision(5);
    if(!fout.fail()){
        if(isglobalEigen){  for(int i=0;i<eigenBe.size();++i){
                fout<<eigenBe[i]<<"\t"<<eigenEd[i]<<"\t"<<gtBe[i]<<"\t"<<gtEd[i]<<endl;
            }
        }else{
            for(int i=0;i<gtBe.size();++i){
                fout<<local_eigenBe[i]<<"\t"<<local_eigenEd[i]<<"\t"<<gtBe[i]<<"\t"<<gtEd[i]<<endl;
            }
        }
    }
    fout.close();

}



int RBF_Core::Solve_Hermite_PredictNormal_UnitNorm_Cluster_MST(double graph_cut_percentage, double localstv_maxcut_percentage, int visualmethod){


    double thres_coef = 0.9;

    arma::vec eigval, ny;
    arma::mat eigvec;

    int iter = 0;
    int maxiter = 100;

    vector<int>fixedInd;

    arma::mat subactiveMat;

    vector<vector<int>>hyper_patch;

    ClearSavedIterBatchInit();


    LocalEigenInit_PredictNormal();
    vector<double>norm_normals(initnormals);
    vector<double>nors_strength(npt);


    auto p_nor = initnormals.data();

    for(int p1=0; p1<npt;++p1){
        nors_strength[p1] = MyUtility::normVec(p_nor+p1*3);
    }

    double maxstrength = *max_element(nors_strength.begin(),nors_strength.end()) * localstv_maxcut_percentage;

    vector<int>headInd;

    vector<int>ccInd;
    vector<bool>headpick(npt,false);
    for(int p1=0; p1<npt;++p1){
        if(nors_strength[p1]>maxstrength){
            headInd.push_back(p1);
            headpick[p1] = true;
        }
        else ccInd.push_back(p1);
    }

    UnionFind uni_find_head;
    UnionFind uni_find_cc;
    uni_find_cc.SetElements(ccInd);
    uni_find_head.SetElements(headInd);

    double westrength = Coherence_Graph_Weight[0]*graph_cut_percentage;
    double n_reserve = 2 * npt * graph_cut_percentage;
    if(0){

        for(int i=0;i<Coherence_Graph_Weight.size();++i){
            if(Coherence_Graph_Weight[i]>westrength){
                int a = Coherence_Graph[i*2],b = Coherence_Graph[i*2+1];
                if(!headpick[a] && !headpick[b]){
                    uni_find_cc.Union(a,b);
                }
                if(headpick[a] && headpick[b]){
                    uni_find_head.Union(a,b);
                }
            }else break;
        }
    }else{

        for(int i=0;i<n_reserve;++i){
            int a = Coherence_Graph[i*2],b = Coherence_Graph[i*2+1];
            if(!headpick[a] && !headpick[b]){
                uni_find_cc.Union(a,b);
            }
            if(headpick[a] && headpick[b]){
                uni_find_head.Union(a,b);
            }
        }
    }


    vector<vector<int>>ccHyper,headHyper;
    uni_find_head.ExtractComponents(headHyper);
    uni_find_cc.ExtractComponents(ccHyper);
    cout<<"Hyper: "<<headHyper.size()<<' '<<ccHyper.size()<<endl;

    for(auto &a:headHyper)hyper_patch.emplace_back(a);
    for(auto &a:ccHyper)hyper_patch.emplace_back(a);


    for(auto &activeInd:headHyper){
        int n_Active = activeInd.size();
        for(int i=0;i<n_Active;++i)MyUtility::normalize(norm_normals.data()+activeInd[i]*3);
        //SaveIterBatchInit(initnormals,norm_normals,activeInd);
    }
    SaveIterBatchInit(initnormals,norm_normals,headInd);
    //    for(auto &activeInd:headHyper){
    //        SaveIterBatchInit(initnormals,norm_normals,activeInd);
    //    }
    cout<<"n_Active: ";
    for(auto &activeInd:headHyper){
        int n_Active = activeInd.size();

        vector<uint>MST_edges;
        vector<vector<int>>MST_v2v;
        vector<double>MST_edgesWeight;
        vector<int>ind(npt,-10000);
        for(int i=0;i<n_Active;++i)ind[activeInd[i]] = i;
        vector<vector<double>>edgelength(npt,vector<double>(npt));
        vector<arma::vec>arma_normals(n_Active,arma::vec(3));
        auto p_nor = norm_normals.data();
        for(int p1=0; p1<n_Active;++p1){
            for(int j=0;j<3;++j)arma_normals[p1](j) = p_nor[activeInd[p1]*3+j];
            arma_normals[p1] = arma::normalise( arma_normals[p1] );
        }

        vector<Triple>eabslength;
        int n_neg = 0;
        for(int m=0; m<n_Active;++m)for(int n=m+1; n<n_Active;++n){
            int p1 = activeInd[m], p2 = activeInd[n];
            arma::mat aa(3,3);
            for(int j=0;j<3;++j)for(int k=0;k<3;++k)aa(j,k) = K(p1+j*npt,p2+k*npt);
            edgelength[m][n] = edgelength[m][n]  = arma::dot(arma_normals[m], aa*arma_normals[n]);
            /*if(  fabs(edgelength[m][n]) > 1e-4 )*/eabslength.emplace_back(fabs(edgelength[m][n]),m,n);
        }
        vector<bool>isflip(n_Active);
        BuildMST(MST_edges, MST_v2v,MST_edgesWeight, eabslength, n_Active);
        MST_Orientate(isflip,edgelength,MST_v2v,n_Active,0);
        for(int i=0;i<n_Active;++i)if(isflip[i]){
            for(int j=0;j<3;++j)initnormals[activeInd[i]*3+j] = -initnormals[activeInd[i]*3+j];
            for(int j=0;j<3;++j)norm_normals[activeInd[i]*3+j] = -norm_normals[activeInd[i]*3+j];
        }



        if(0){
            arma::uvec row_ind(n_Active*3);
            for(int i=0;i<n_Active;++i){
                for(int j=0;j<3;++j)row_ind(i*3+j) = activeInd[i]+j*npt;
            }
            //cout<<row_ind.t()<<endl;
            subactiveMat = K.submat(row_ind,row_ind);
            //cout<<K.n_rows<<" "<<K.n_cols<<endl;
            ny = eig_sym( eigval, eigvec, subactiveMat);
            arma::vec firstEigvec(eigvec.col(0));
            for(int i=0;i<n_Active;++i){
                auto p_nor = initnormals.data()+activeInd[i]*3;
                auto p_nor_n = norm_normals.data()+activeInd[i]*3;
                for(int j=0;j<3;++j)p_nor_n[j] = p_nor[j] = firstEigvec(i*3+j);
                MyUtility::normalize(p_nor_n);
            }
        }
        cout<<n_Active<<' ';
        SaveIterBatchInit(initnormals,norm_normals,activeInd);
    }
    cout<<" || ";
    vector<double>saveInit(initnormals),save_norm(norm_normals);



    for(auto &ccInd:ccHyper){
        UnionFind uni_find;
        vector<int>ele;
        for(auto &b:headHyper)for(auto a:b)ele.push_back(a);
        for(auto b:ccInd)ele.push_back(b);
        uni_find.SetElements(ele);
        for(auto &b:headHyper)for(int i=1;i<b.size();++i)uni_find.Union(b[i],b[i-1]);
        for(int i=1;i<ccInd.size();++i)uni_find.Union(ccInd[i],ccInd[i-1]);
        vector<bool>ccpick(npt,false);
        for(auto b:ccInd)ccpick[b]=true;
        for(int i=0;i<n_reserve;++i){
            int a = Coherence_Graph[i*2],b = Coherence_Graph[i*2+1];
            if((headpick[a] && ccpick[b])|| (headpick[b] && ccpick[a])){
                uni_find.Union(a,b);
            }
        }
        vector<int>activeInd(ccInd);
        for(auto &b:headHyper)if(uni_find.Is_Connected(b[0],ccInd[0]))for(auto a:b)activeInd.push_back(a);

        int n_Active = activeInd.size();
        arma::uvec row_ind(n_Active*3);
        for(int i=0;i<n_Active;++i){
            for(int j=0;j<3;++j)row_ind(i*3+j) = activeInd[i]+j*npt;
        }
        //cout<<row_ind.t()<<endl;
        subactiveMat = K.submat(row_ind,row_ind);
        //cout<<K.n_rows<<" "<<K.n_cols<<endl;
        ny = eig_sym( eigval, eigvec, subactiveMat);
        arma::vec firstEigvec(eigvec.col(0));
        vector<int>*drawInd;
        if(visualmethod==0)drawInd = &activeInd;//activeInd ccInd
        else if(visualmethod==1)drawInd = &ccInd;
        for(int i=0;i<drawInd->size();++i){
            auto p_nor = initnormals.data()+activeInd[i]*3;
            auto p_nor_n = norm_normals.data()+activeInd[i]*3;
            for(int j=0;j<3;++j)p_nor_n[j] = p_nor[j] = firstEigvec(i*3+j);
            MyUtility::normalize(p_nor_n);
        }
        cout<<n_Active<<' ';
        SaveIterBatchInit(initnormals,norm_normals,*drawInd);
    }
    cout<<endl;
    for(auto &activeInd:headHyper){
        int n_Active = activeInd.size();
        for(int i=0;i<n_Active;++i){
            int offset = activeInd[i]*3;
            MyUtility::copyVec(saveInit.data()+offset, initnormals.data()+offset);
            MyUtility::copyVec(save_norm.data()+offset, norm_normals.data()+offset);
        }
    }



    vector<Triple>eweights;
    int n_hyper = hyper_patch.size();

    vector<arma::vec>hyperb_vec(hyper_patch.size());
    vector<arma::uvec>hyperb_rowind(hyper_patch.size());
    for(int i=0;i<hyper_patch.size();++i){
        int b_nv = hyper_patch[i].size();
        auto &hyperVec = hyperb_vec[i];
        auto &row_ind = hyperb_rowind[i];
        auto &hyper_ind = hyper_patch[i];
        hyperVec.set_size(b_nv*3);
        row_ind.set_size(b_nv*3);

        for(int k=0;k<b_nv;++k){
            for(int j=0;j<3;++j)hyperVec(k*3+j) = norm_normals[hyper_ind[k]*3+j];
            for(int j=0;j<3;++j)row_ind(k*3+j) = hyper_ind[k]+j*npt;
        }
    }
    cout<<"adsada"<<endl;

    vector<vector<double>>interactiveEn_mat(n_hyper,vector<double>(n_hyper));
    for(int i=0;i<hyper_patch.size();++i){
        for(int j=i+1;j<hyper_patch.size();++j){
            arma::mat subInteractiveMat = K.submat(hyperb_rowind[i],hyperb_rowind[j]);
            //cout<<"aaa"<<endl;
            interactiveEn_mat[j][i] = interactiveEn_mat[i][j] = arma::dot(hyperb_vec[i],subInteractiveMat*hyperb_vec[j]);
            //cout<<i<<" "<<j<<" "<<interactiveEn_mat[i][j]<<endl;
            eweights.emplace_back(fabs(interactiveEn_mat[i][j]),i,j);
        }
    }


    vector<bool>isflip(n_hyper,false);

    if(1){
        cout<<"BuildMST"<<endl;
        vector<uint>MST_edges;
        vector<vector<int>>MST_v2v;
        vector<double>MST_edgesWeight;
        BuildMST(MST_edges, MST_v2v,MST_edgesWeight, eweights, n_hyper);
        MST_Orientate(isflip,interactiveEn_mat,MST_v2v,n_hyper,0);

    }else{
        ExhaustedCombination(isflip,interactiveEn_mat,n_hyper);
    }


    for(int i=0;i<n_hyper;++i)if(isflip[i]){
        int b_nv = hyper_patch[i].size();
        auto &hyper_ind = hyper_patch[i];
        for(int k=0;k<b_nv;++k){
            for(int j=0;j<3;++j)norm_normals[hyper_ind[k]*3+j] = -norm_normals[hyper_ind[k]*3+j];
        }
    }


    initnormals_uninorm = initnormals = norm_normals;
    this->MST_edges.clear();
    this->MST_edgesWeight.clear();
    cout<<"Solve_Hermite_PredictNormal_UnitNorm_Cluster_MST"<<endl;
    return 1;


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
        }else if(curMethod==Hermite_Tangent_UnitNormal){
            Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm();
        }else if(curMethod==HandCraft){
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


int RBF_Core::Solve_Hermite_PredictNormal_UnitNormal_LamnbdaSearchTest(){

    cout<<"Solve_Hermite_PredictNormal_UnitNormal_LamnbdaSearchTest"<<endl;
    curMethod = Hermite_UnitNormal;
    Set_Hermite_PredictNormal(pts);
    Lamnbda_Search_GlobalEigen();

    RBF_Paras para;
    para.InitMethod = GT_NORMAL;
    InitNormal(para);
    OptNormal(1);

    gtBe.push_back(sol.init_energy);
    gtEd.push_back(sol.energy);

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



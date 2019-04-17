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
#include "ImplicitedSurfacing.h"
typedef std::chrono::high_resolution_clock Clock;
int RBF_Core::Invoke(vector<double>&pts, vector<int>&labels,  vector<double>&normals, vector<double>&tangents, vector<uint>&edges, RBF_Paras para){


    //SparseMatrixTest();
    isuse_sparse = true;

    this->pts = pts;
    this->labels = labels;
    this->normals = normals;
    this->tangents = tangents;
    this->edges = edges;
    npt = this->pts.size()/3;
    curMethod = para.Method;

    polyDeg = para.polyDeg;
    User_Lamnbda = para.user_lamnbda;
    rangevalue = para.rangevalue;
    maxvalue = 10000;
    Hermite_weight_smoothness = para.Hermite_weight_smoothness;

    Hermite_designcurve_weight = para.Hermite_designcurve_weight;

    Set_Actual_Hermite_LSCoef(para.Hermite_ls_weight);

    cout<<"number of points: "<<pts.size()/3<<endl;
    cout<<"normals: "<<this->normals.size()<<endl;
    sol.Statue = 1;
    Init(para.Kernal);

    SetSigma(para.sigma);

    SetThis();

    double ttime;
    auto t1 = Clock::now();

    switch(para.Method){
    case Variational:
        Set_VariationalRBF(pts,labels);
        Solve_VariationalRBF();
        break;
    case Variational_P:
        Set_VariationalRBF_P(pts, labels);
        Solve_VariationalRBF_P();
        break;
    case LS:
        Set_LSRBF(pts, labels, false);
        Solve_RBF();
        break;
    case LSinterp:
        Set_LSRBF(pts, labels, true);
        Solve_RBF();
        break;
    case Interp:
        Set_InterpRBF(pts,labels);
        Solve_RBF();
        break;
    case RayleighQuotients:
        Set_RayleighQuotients(pts,labels);
        Solve_RayleighQuotients();
        break;
    case RayleighQuotients_P:
        Set_RayleighQuotients_P(pts,labels);
        Solve_RayleighQuotients_P();
        break;
    case RayleighQuotients_I:
        Set_RayleighQuotients(pts,labels);
        IterativelySolveRQ();
        break;
    case Hermite:
        Set_HermiteRBF(pts);
        Solve_HermiteRBF(normals);
        break;
    case Hermite_UnitNorm:
        Set_Hermite_PredictNormal(pts);
        Solve_Hermite_PredictNormal_UnitNorm();
        break;
    case Hermite_UnitNormal:
        Set_Hermite_PredictNormal(pts);
        cout << "Setup time: " << (setup_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9) << endl;
        //LocalEigenInit_PredictNormal();
        Solve_Hermite_PredictNormal_UnitNormal();
        break;
    case Hermite_InitializationTest:

        Solve_Hermite_PredictNormal_UnitNormal_LamnbdaSearchTest();
//        Set_Hermite_PredictNormal(pts);
//        cout << "Setup time: " << (setup_time = std::chrono::nanoseconds(Clock::now() - t1).count()/1e9) << endl;
//        Solve_Hermite_PredictNormal_UnitNormal_InitializationTest();
//        //Solve_Hermite_PredictNormal_UnitNormal_SolverTest();
        break;
    case Hermite_Tangent_UnitNorm:
        Set_Hermite_PredictNormal_TangentConstraint(pts,tangents,edges);
        Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm();
        break;
    case Hermite_Tangent_UnitNormal:
        Set_Hermite_PredictNormal_TangentConstraint(pts,tangents,edges);
        Solve_Hermite_PredictNormal_TangentConstraint_UnitNormal();
        break;
    }
    auto t2 = Clock::now();
    cout << "Total opt time: " << (ttime = std::chrono::nanoseconds(t2 - t1).count()/1e9) << endl<< endl;

    //cout<<a<<endl;
    //cout<<b<<endl;

    Record(para.Method,para.Kernal,sol,ttime);

    return sol.Statue;

}


void RBF_Core::BuildK(RBF_Paras para){

    isuse_sparse = para.isusesparse;
    sparse_para = para.sparse_para;
    Hermite_weight_smoothness = para.Hermite_weight_smoothness;
    Hermite_designcurve_weight = para.Hermite_designcurve_weight;
    handcraft_sigma = para.handcraft_sigma;
    wDir = para.wDir;
    wOrt = para.wOrt;
    wFlip = para.wFlip;
    curMethod = para.Method;

    Set_Actual_Hermite_LSCoef( para.Hermite_ls_weight );
    Set_Actual_User_LSCoef(  para.user_lamnbda  );
    isNewApprox = true;
    isnewformula = true;

    auto t1 = Clock::now();

    switch(curMethod){
    case Variational:
        Set_VariationalRBF(pts,labels);
        break;
    case Variational_P:
        Set_VariationalRBF_P(pts, labels);
        break;
    case LS:
        Set_LSRBF(pts, labels, false);
        break;
    case LSinterp:
        Set_LSRBF(pts,labels, true);
        break;
    case Interp:
        Set_InterpRBF(pts,labels);
        break;
    case RayleighQuotients:
        Set_RayleighQuotients(pts,labels);
        break;
    case RayleighQuotients_P:
        Set_RayleighQuotients_P(pts,labels);
        break;
    case RayleighQuotients_I:
        Set_RayleighQuotients(pts,labels);
        break;
    case Hermite:
        Set_HermiteRBF(pts);
        break;
    case Hermite_UnitNorm:
        Set_Hermite_PredictNormal(pts);
        break;
    case Hermite_UnitNormal:
        Set_Hermite_PredictNormal(pts);
        break;
    case Hermite_InitializationTest:
        Set_Hermite_PredictNormal(pts);
        break;
    case Hermite_Tangent_UnitNorm:
        Set_Hermite_PredictNormal_TangentConstraint(pts,tangents,edges);
        break;
    case Hermite_Tangent_UnitNormal:
        Set_Hermite_PredictNormal_TangentConstraint(pts,tangents,edges);
        break;
    case HandCraft:
        HandCraftK();
        break;
    }
    auto t2 = Clock::now();
    cout << "Build Time: " << (setup_time = std::chrono::nanoseconds(t2 - t1).count()/1e9) << endl<< endl;


    if(0)BuildCoherentGraph();
}

void RBF_Core::InitNormal(RBF_Paras para){


    auto t1 = Clock::now();
    curInitMethod = para.InitMethod;
    cout<<"Init Method: "<<mp_RBF_INITMETHOD[curInitMethod]<<endl;
    switch(curInitMethod){

    case RBF_Init_EMPTY:
        cout<<"No Initialization"<<endl;
        break;
    case GlobalEigen:
        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm();
        }else if(curMethod==Hermite_Tangent_UnitNormal){
            Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm();
        }else if(curMethod==HandCraft){
            Solve_Hermite_PredictNormal_UnitNorm();
        }
        break;
    case GlobalEigenWithMST:
        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm();
        }else if(curMethod==Hermite_Tangent_UnitNormal){
            Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm();
        }else if(curMethod==HandCraft){
            Solve_Hermite_PredictNormal_UnitNorm();
        }
        NormalReOrientation(initnormals,false);
        break;
    case GlobalEigenWithGT:
        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm();
        }else if(curMethod==Hermite_Tangent_UnitNormal){
            Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm();
        }else if(curMethod==HandCraft){
            Solve_Hermite_PredictNormal_UnitNorm();
        }
        NormalReOrientation_WithGT(initnormals,false);
        break;
    case LocalEigen:
        if(curMethod==Hermite_UnitNormal){
            LocalEigenInit_PredictNormal();
        }else if(curMethod==Hermite_Tangent_UnitNormal){
            LocalEigenInit_PredictNormal_TangentConstraint();
        }else if(curMethod==HandCraft){
            LocalEigenInit_PredictNormal();
        }
        break;
    case IterativeEigen:
        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm_Iterative_MST();
        }else if(curMethod==Hermite_Tangent_UnitNormal){
            Solve_Hermite_PredictNormal_TangentConstraint_UnitNorm_Iterative_MST();
        }else if(curMethod==HandCraft){
            Solve_Hermite_PredictNormal_UnitNorm_Iterative_MST();
        }
        break;
    case ClusterEigen:
        if(curMethod==Hermite_UnitNormal){
            Solve_Hermite_PredictNormal_UnitNorm_Cluster_MST(para.ClusterCut_percentage,para.ClusterCut_LocalMax_percentage,para.ClusterVisualMethod);
        }else if(curMethod==Hermite_Tangent_UnitNormal){


        }else if(curMethod==HandCraft){
            Solve_Hermite_PredictNormal_UnitNorm_Cluster_MST(para.ClusterCut_percentage,para.ClusterCut_LocalMax_percentage,para.ClusterVisualMethod);
        }
        break;

    case GT_NORMAL:
        initnormals = normals;
        SetInitnormal_Uninorm();
        newnormals = initnormals_uninorm;
        break;

    case Lamnbda_Search:
        Lamnbda_Search_GlobalEigen();
        break;

    case GlobalMRF:
        GlobalMRFOrientation();
        break;

    case Voronoi_Covariance:
        Compute_Voronoi_Covariance();
        break;
    case CNN:
        ComputeCNN();
        break;
    case PCA:
        ComputePCA();
        break;
    }



    auto t2 = Clock::now();
    cout << "Init Time: " << (init_time = std::chrono::nanoseconds(t2 - t1).count()/1e9) << endl<< endl;

    mp_RBF_InitNormal[curMethod==HandCraft?0:1][curInitMethod] = initnormals;

}

void RBF_Core::OptNormal(int method){

    cout<<"OptNormal"<<endl;
    auto t1 = Clock::now();


    switch(curMethod){

    case HandCraft:
    case Hermite_UnitNormal:
        Opt_Hermite_PredictNormal_UnitNormal();
        break;
    case Hermite_Tangent_UnitNormal:
        Opt_Hermite_PredictNormal_TangentConstraint_UnitNormal();
        break;
    }
    auto t2 = Clock::now();
    cout << "Opt Time: " << (solve_time = std::chrono::nanoseconds(t2 - t1).count()/1e9) << endl<< endl;
    if(method==0)mp_RBF_OptNormal[curMethod==HandCraft?0:1][curInitMethod] = newnormals;
}


void RBF_Core::Surfacing(int method){

    n_evacalls = 0;
    Surfacer sf;
    double re_time;
    if(method==0)re_time = sf.Surfacing_Implicit(pts,labels,true,RBF_Core::Dist_Function);
    else if(method==1)re_time = sf.Surfacing_CGALImplicit(RBF_Core::Dist_Function);
    else if(method==2)re_time = sf.Surfacing_Spectral(pts,newnormals);
    else if(method==3)re_time = sf.Surfacing_Poisson(pts,newnormals);
    else if(method==4)re_time = sf.Surfacing_PowerCrust(pts);
    else if(method==5)re_time = sf.Surfacing_ScreenedPoisson(pts,newnormals);
    else if(method==6)re_time = sf.Surfacing_NASR(pts);
    else if(method==7){
        isnewformula = false;
        Set_HermiteRBF(pts);
        Solve_HermiteRBF(newnormals);
        re_time = sf.Surfacing_Implicit(pts,labels,true,RBF_Core::Dist_Function);
        isnewformula = true;
    }

    sf.WriteSurface(finalMesh_v,finalMesh_fv);

    cout<<"n_evacalls: "<<n_evacalls<<"   ave: "<<re_time/n_evacalls<<endl;


}

int RBF_Core::InjectData(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para){

    isuse_sparse = para.isusesparse;
    sparse_para = para.sparse_para;
    //isuse_sparse = false;
    this->pts = pts;
    this->labels = labels;
    this->normals = normals;
    this->tangents = tangents;
    this->edges = edges;
    npt = this->pts.size()/3;
    curMethod = para.Method;
    curInitMethod = para.InitMethod;

    polyDeg = para.polyDeg;
    User_Lamnbda = para.user_lamnbda;
    rangevalue = para.rangevalue;
    maxvalue = 10000;

    cout<<"number of points: "<<pts.size()/3<<endl;
    cout<<"normals: "<<this->normals.size()<<endl;
    sol.Statue = 1;
    Init(para.Kernal);

    SetSigma(para.sigma);

    SetThis();

	return 1;
}

int RBF_Core::ThreeStep(vector<double>&pts, vector<int>&labels, vector<double>&normals, vector<double>&tangents,  vector<uint>&edges, RBF_Paras para){

    InjectData(pts, labels, normals, tangents, edges,  para);
    BuildK(para);
    InitNormal(para);
    OptNormal(0);
	
	return 1;
}


int RBF_Core::AllStep(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para){

    InjectData(pts, labels, normals, tangents, edges,  para);
    BuildK(para);
    InitNormal(para);
    OptNormal(0);
    Surfacing(0);
	return 1;

}

void RBF_Core::BatchInitEnergyTest(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para){

    InjectData(pts, labels, normals, tangents, edges,  para);
    BuildK(para);
    para.ClusterVisualMethod = 0;//RBF_Init_EMPTY
    for(int i=0;i<RBF_Init_EMPTY;++i){
        para.InitMethod = RBF_InitMethod(i);
        InitNormal(para);
        OptNormal(0);
        Record();
    }
    Print_Record_Init();
}


vector<double>* RBF_Core::ExportPts(){

    return &pts;



}

vector<double>* RBF_Core::ExportPtsNormal(int normal_type){

    if(normal_type==0)return &normals;
    else if(normal_type==1)return &initnormals;
    else if(normal_type==2)return &initnormals_uninorm;
    else if(normal_type==3)return &newnormals;

	return NULL;
}

vector<double>* RBF_Core::ExportGraphWeight(){
    return &Coherence_Graph_Weight;
}

vector<uint>* RBF_Core::ExportGraph(){
    return &Coherence_Graph;
}


vector<double>* RBF_Core::ExportMSTWeight(){

    return &MST_edgesWeight;

}

vector<uint>* RBF_Core::ExportMST(){
    return &MST_edges;

}

vector<double>* RBF_Core::ExportInitNormal(int kmethod, RBF_InitMethod init_type){

    if(mp_RBF_InitNormal[kmethod].find(init_type)!=mp_RBF_InitNormal[kmethod].end())return &(mp_RBF_InitNormal[kmethod][init_type]);
    else return NULL;
}

vector<double>* RBF_Core::ExportOptNormal(int kmethod, RBF_InitMethod init_type){

    if(mp_RBF_OptNormal[kmethod].find(init_type)!=mp_RBF_OptNormal[kmethod].end())return &(mp_RBF_OptNormal[kmethod][init_type]);
    else return NULL;
}



void RBF_Core::Print_Record_Init(){

    cout<<"InitMethod"<<string(30-string("InitMethod").size(),' ')<<"InitEn\t\t FinalEn"<<endl;
    cout<<std::setprecision(8);
    {
        for(int i=0;i<record_initmethod.size();++i){
            cout<<record_initmethod[i]<<string(30-record_initmethod[i].size(),' ')<<record_initenergy[i]<<"\t\t"<<record_energy[i]<<endl;
        }
    }

}

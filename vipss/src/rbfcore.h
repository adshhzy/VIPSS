#ifndef RBFCORE_H
#define RBFCORE_H




#include <iostream>
#include <vector>
#include "Solver.h"
#include "ImplicitedSurfacing.h"
//#include "eigen3/Eigen/Dense"
#include <armadillo>
#include <unordered_map>
using namespace std;

enum RBF_INPUT{
    ON,
    ONandNORMAL,
    ALL,
    INandOUT,
};


enum RBF_METHOD{
    Variational,
    Variational_P,
    LS,
    LSinterp,
    Interp,
    RayleighQuotients,
    RayleighQuotients_P,
    RayleighQuotients_I,
    Hermite,
    Hermite_UnitNorm,
    Hermite_UnitNormal,
    Hermite_Tangent_UnitNorm,
    Hermite_Tangent_UnitNormal,
    HandCraft,
    Hermite_InitializationTest
};

enum RBF_InitMethod{
    GT_NORMAL,
    GlobalEigen,
    GlobalEigenWithMST,
    GlobalEigenWithGT,
    LocalEigen,
    IterativeEigen,
    ClusterEigen,
    Lamnbda_Search,
    GlobalMRF,
    Voronoi_Covariance,
    CNN,
    PCA,
    RBF_Init_EMPTY
};

enum RBF_Kernal{
    XCube,
    ThinSpline,
    XLinear,
    Gaussian,
};

class RBF_Paras{
public:
    RBF_METHOD Method;
    RBF_Kernal Kernal;
    RBF_InitMethod InitMethod;
    bool isusesparse;
    int polyDeg;
    double sigma;
    double user_lamnbda;
    double rangevalue;
    double sparse_para = 1e-3;
    double Hermite_weight_smoothness;
    double Hermite_ls_weight;
    double Hermite_designcurve_weight;
    double ClusterCut_percentage;
    double ClusterCut_LocalMax_percentage;
    int ClusterVisualMethod;
    double wDir,wOrt,wFlip,handcraft_sigma;
    RBF_Paras(RBF_METHOD Method, RBF_Kernal Kernal, int polyDeg, double sigma, double user_lamnbda,double rangevalue):\
        Method(Method),Kernal(Kernal),InitMethod(RBF_Init_EMPTY),polyDeg(polyDeg),sigma(sigma),user_lamnbda(user_lamnbda),rangevalue(rangevalue),Hermite_weight_smoothness(0){}
    RBF_Paras(){}
};



class RBF_Core{

public:



    int npt;
    int polyDeg = 2;
    int bsize;

    bool isinv = true;
    bool isnewformula = true;
    double User_Lamnbda;

    RBF_Kernal kernal;
    RBF_METHOD curMethod;
    RBF_InitMethod curInitMethod;

    double rangevalue = 0.2;
    double maxvalue = 10000;

    vector<double>pts;
    vector<double>normals;
    vector<double>tangents;
    vector<uint>edges;
    vector<int>labels;

private:

    vector<double>initnormals;
    vector<double>initnormals_uninorm;

public:
    vector<double>newnormals;


public:
    vector<double>coff_cos;
    vector<double>coff_sin;


public:
    vector<double>finalMesh_v;
    vector<uint>finalMesh_fv;

public:

    arma::mat M;
    arma::mat N;

    arma::vec a;
    arma::vec b;

    arma::mat Minv;
    arma::mat P;
    arma::mat K;
    arma::mat bprey;
    arma::mat saveK;
    arma::mat saveK_finalH;
    arma::mat finalH;

    arma::mat RQ;

    arma::mat bigM;
    arma::mat bigMinv;
    arma::mat Ninv;
    arma::mat K00;
    arma::mat K01;
    arma::mat K11;
    arma::mat dI;


    bool isuse_sparse = false;
    double sparse_para = 1e-3;

public:
    unordered_map<int, string>mp_RBF_INITMETHOD;
    unordered_map<int, string>mp_RBF_METHOD;
    unordered_map<int, string>mp_RBF_Kernal;

public:

    double (*Kernal_Function)(const double x);

    double (*Kernal_Function_2p)(const double *p1, const double *p2);

    double (*P_Function_2p)(const double *p1, const double *p2);

public:

    bool isNewApprox;
    bool isHermite;
    double Hermite_weight_smoothness;
    double Hermite_ls_weight_inject, User_Lamnbda_inject;
    double Hermite_designcurve_weight;

    double ls_coef;
    void (*Kernal_Gradient_Function_2p)(const double *p1, const double *p2, double *G);
    void (*Kernal_Hessian_Function_2p)(const double *p1, const double *p2, double *H);

private:
    vector<double>local_eigenBe, local_eigenEd, eigenBe, eigenEd, gtBe, gtEd;

private:
    vector<vector<double> >lamnbdaGlobal_Be,lamnbdaGlobal_Ed;
    vector<double>lamnbda_list_sa;
private:
    unordered_map<int,unordered_map<int, vector<double> > >mp_RBF_InitNormal;
    unordered_map<int,unordered_map<int, vector<double> > >mp_RBF_OptNormal;

public:

    RBF_Core();
    RBF_Core(RBF_Kernal kernal);
    void Init(RBF_Kernal kernal);


public:

    double Dist_Function(const double x, const double y, const double z);
    inline double Dist_Function(const double *p);

public:
    static double Dist_Function(const R3Pt &in_pt);
    //static FT Dist_Function(const Point_3 in_pt);
    int n_evacalls;
public:
    void SetThis();
public:

    void SetSigma(double x);


public:


    void NormalRecification(double maxlen, vector<double> &nors);



public:
    void Set_HermiteRBF(vector<double>&pts);
    int Solve_HermiteRBF(vector<double>&vn);

public:
    void Set_Hermite_PredictNormal(vector<double>&pts);

public:

    int Solve_Hermite_PredictNormal_UnitNorm();


    int Lamnbda_Search_GlobalEigen();


public:
    int Solve_Hermite_PredictNormal_UnitNormal();

public:

    void Print_LamnbdaSearchTest(string fname);



public:


    void Set_RBFCoef(arma::vec &y);

    void Set_Actual_Hermite_LSCoef(double hermite_ls);
    void Set_HermiteApprox_Lamnda(double hermite_ls);
    void Set_Actual_User_LSCoef(double user_ls);
    void Set_User_Lamnda_ToMatrix(double user_ls);

    void Set_SparsePara(double spa);


public:

    bool Write_Hermite_NormalPrediction(string fname,int mode);
    bool Write_Hermite_MST(string fname);
    void WriteSeletivePLY(string fname, vector<double>&allnormals, vector<int>&pickInd);

    void SetInitnormal_Uninorm();
    void ClearSavedIterBatchInit();
    void SaveIterBatchInit(vector<double>&allnormals, vector<double>&allnormals_uninorm, vector<int>&pickInd);

    void Write_Surface(string fname);

public:

    int Opt_Hermite_PredictNormal_UnitNormal();

public:

    int ThreeStep(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para);
    int AllStep(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para);

    int InjectData(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para);

    int InjectData(vector<double> &pts, RBF_Paras para);

    void BuildK(RBF_Paras para);

    void InitNormal(RBF_Paras para);

    void OptNormal(int method);

    void Surfacing(int method, int n_voxels_1d);

    void BuildCoherentGraph();

    void BatchInitEnergyTest(vector<double> &pts, vector<int> &labels, vector<double> &normals, vector<double> &tangents, vector<uint> &edges, RBF_Paras para);


public:
    vector<double>* ExportPts();
    vector<double>* ExportPtsNormal(int normal_type);


    vector<double>* ExportInitNormal(int kmethod, RBF_InitMethod init_type);
    vector<double>* ExportOptNormal(int kmethod, RBF_InitMethod init_type);



public:

    Solution_Struct sol;

    vector<int>record_partition;
    vector<string>record_partition_name;

    vector<int>npoints;

    vector<string>record_initmethod;
    vector<string>record_method;
    vector<string>record_kernal;
    vector<double>record_initenergy;
    vector<double>record_energy;
    vector<double>record_time;

    double setup_time, init_time, solve_time, callfunc_time,invM_time, setK_time, surf_time;
    vector<double>setup_timev, init_timev, solve_timev, callfunc_timev,invM_timev,setK_timev;

    void Record();
    void Record(RBF_METHOD method, RBF_Kernal kernal, Solution_Struct &rsol, double time);
    void AddPartition(string pname);
    void Print_Record();
    void Print_TimerRecord(string fname);
    void Clear_TimerRecord();
    void Print_Record_Init();

    void Print_TimerRecord_Single(string fname);

};














#endif // RBFCORE_H

#include <iostream>
#include <unistd.h>
#include "src/rbfcore.h"
#include "src/readers.h"
using namespace std;



void SplitPath(const std::string& fullfilename,std::string &filepath);
void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname);
RBF_Paras Set_RBF_PARA();
int main(int argc, char** argv)
{
    cout << argc << endl;



    string infilename;
    string outpath, pcname, ext, inpath;

    int n_voxel_line = 100;

    double user_lambda = 0;

    bool is_surfacing = false;
    bool is_outputtime = false;

    int c;
    optind=1;
    while ((c = getopt(argc, argv, "i:o:l:s:t")) != -1) {
        switch (c) {
        case 'i':
            infilename = optarg;
            break;
        case 'o':
            outpath = string(optarg);
            break;
        case 'l':
            user_lambda = atof(optarg);
            break;
        case 's':
            is_surfacing = true;
            n_voxel_line = atoi(optarg);
            break;
        case 't':
            is_outputtime = true;
            break;
        case '?':
            cout << "Bad argument setting!" << endl;
            break;
        }
    }

    if(outpath.empty())SplitFileName(infilename,outpath,pcname,ext);
    else SplitFileName(infilename,inpath,pcname,ext);
    cout<<"input file: "<<infilename<<endl;
    cout<<"output path: "<<outpath<<endl;

    cout<<"user lambda: "<<user_lambda<<endl;
    cout<<"is surfacing: "<<is_surfacing<<endl;

    cout<<"number of voxel per D: "<<n_voxel_line<<endl;


    vector<double>Vs;
    RBF_Core rbf_core;
    RBF_Paras para = Set_RBF_PARA();
    para.user_lamnbda = user_lambda;

    readXYZ(infilename,Vs);
    rbf_core.InjectData(Vs,para);
    rbf_core.BuildK(para);
    rbf_core.InitNormal(para);
    rbf_core.OptNormal(0);

    rbf_core.Write_Hermite_NormalPrediction(outpath+pcname+"_normal", 1);

    if(is_surfacing){
        rbf_core.Surfacing(0,n_voxel_line);
        rbf_core.Write_Surface(outpath+pcname+"_surface");
    }

    if(is_outputtime){
        rbf_core.Print_TimerRecord_Single(outpath+pcname+"_time.txt");
    }








    return 0;
}







RBF_Paras Set_RBF_PARA(){

    RBF_Paras para;
    RBF_InitMethod initmethod = Lamnbda_Search;

    RBF_Kernal Kernal = XCube;
    int polyDeg = 1;
    double sigma = 0.9;
    double rangevalue = 0.001;

    para.Kernal = Kernal;para.polyDeg = polyDeg;para.sigma = sigma;para.rangevalue = rangevalue;
    para.Hermite_weight_smoothness = 0.0;
    para.Hermite_ls_weight = 0;
    para.Hermite_designcurve_weight = 00.0;
    para.Method = RBF_METHOD::Hermite_UnitNormal;


    para.InitMethod = initmethod;

    para.user_lamnbda = 0;

    para.isusesparse = false;


    return para;
}


inline void SplitFileName (const std::string& fullfilename,std::string &filepath,std::string &filename,std::string &extname) {
    int pos;
    pos = fullfilename.find_last_of('.');
    filepath = fullfilename.substr(0,pos);
    extname = fullfilename.substr(pos);
    pos = filepath.find_last_of("\\/");
    filename = filepath.substr(pos+1);
    pos = fullfilename.find_last_of("\\/");
    filepath = fullfilename.substr(0,pos+1);
    //cout<<modname<<' '<<extname<<' '<<filepath<<endl;
}

void SplitPath(const std::string& fullfilename,std::string &filepath){

    std::string filename;
    std::string extname;
    SplitFileName(fullfilename, filepath, filename,extname);
}

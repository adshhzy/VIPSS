#include"readers.h"
#include<iostream>
#include <iomanip>
#include<fstream>
#include<sstream>
#include<assert.h>
using namespace std;

bool readOffFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices){
    ifstream reader(filename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open the OFF file " << filename << endl;
        return false;
    }else {
        cout << "Reading: "<<filename<<endl;
    }

    string ss;

    int ivalue,ibuf[20];
    reader>>ss;
    //cout<<ss<<endl;
    if(ss!="OFF"){
        cout << "Not OFF file: " << filename << endl;
        return false;
    }


    int n_vertices,n_faces,n_edges;
    reader>>n_vertices;
    reader>>n_faces;
    reader>>n_edges;

    cout<<n_vertices<<' '<< n_faces<<' '<<n_edges<<endl;
    vertices.resize(n_vertices*3);
    //faces2vertices.resize(n_faces*3);

    for(int i =0;i<vertices.size();i++){
        reader>>vertices[i];
    }

    faces2vertices.clear();
    for(int i =0;i<n_faces;i++){
        reader>>ivalue;
        vector<int>tempvlist(ivalue);
        for(int j=0;j<ivalue;++j)reader>>tempvlist[j];
        for(int i=2;i<tempvlist.size();++i){
            faces2vertices.push_back(tempvlist[0]);
            faces2vertices.push_back(tempvlist[i-1]);
            faces2vertices.push_back(tempvlist[i]);
        }
    }




    reader.close();
    return true;



}


bool writeOffFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices){
    filename = filename + ".off";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output OFF file " << filename << endl;
        return false;
    }


    int n_vertices = vertices.size()/3;
    int n_faces = faces2vertices.size()/3;

    outer<<"OFF"<<endl;
    outer<<n_vertices<<' '<<n_faces<<' '<<0<<endl;
    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }
    for(int i=0;i<n_faces;++i){
        auto p_fv = faces2vertices.data()+i*3;
        outer << "3 " << p_fv[0]+1<< " "<< p_fv[1]+1 << " "<< p_fv[2]+1 << endl;
    }

    outer.close();
    cout<<"saving finish: "<<filename<<endl;
    return true;

}

bool writePLYFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices,
                  const vector<double>&vertices_normal,const vector<unsigned char>&vertices_color){
    filename = filename + ".ply";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output PLY file " << filename << endl;
        return false;
    }


    int n_vertices = vertices.size()/3;
    int n_faces = faces2vertices.size()/3;
    outer << "ply" <<endl;
    outer << "format ascii 1.0"<<endl;
    outer << "element vertex " << n_vertices <<endl;
    outer << "property float x" <<endl;
    outer << "property float y" <<endl;
    outer << "property float z" <<endl;
    outer << "property float nx" <<endl;
    outer << "property float ny" <<endl;
    outer << "property float nz" <<endl;
    outer << "property uchar red" <<endl;
    outer << "property uchar green" <<endl;
    outer << "property uchar blue" <<endl;
    outer << "property uchar alpha" <<endl;
    outer << "element face " << n_faces <<endl;
    outer << "property list uchar int vertex_indices" <<endl;
    outer << "end_header" <<endl;

    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        auto p_vn = vertices_normal.data()+i*3;
        auto p_vc = vertices_color.data()+i*4;
        for(int j=0;j<3;++j)outer << p_v[j] << " ";
        for(int j=0;j<3;++j)outer << p_vn[j] << " ";
        for(int j=0;j<4;++j)outer << int(p_vc[j]) << " ";
        outer << endl;
    }

    for(int i=0;i<n_faces;++i){
        auto p_fv = faces2vertices.data()+i*3;
        outer << "3 ";
        for(int j=0;j<3;++j)outer << p_fv[j] << " ";
        outer << endl;
    }
    outer.close();
    cout<<"saving finish: "<<filename<<endl;
    return true;
}

bool readPLYFile(string filename,  vector<double>&vertices, vector<double> &vertices_normal){
    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    vertices.clear();
    vertices_normal.clear();
    auto readVerticesAndNormal = [&vertices,&vertices_normal](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
        for(int i=0;i<3;++i){strs>>dvalue;vertices_normal.push_back(dvalue);}
    };

    string oneline;

    cout<<"reading: "<<filename<<endl;
    bool isstart = false;
    while( getline( fin, oneline ) ){
        stringstream strs( oneline );
        string prefix;

        if(isstart){
            readVerticesAndNormal( strs ); continue;
        }else{
            strs >> prefix;
            if( prefix == "end_header"  ) {isstart = true; continue; }
        }
    }


    fin.close();
    return true;



}


bool readObjFile(string filename, vector<double>&vertices, vector<unsigned int>&faces2vertices, vector<double>&vertices_normal){

    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    vertices.clear();
    faces2vertices.clear();
    vertices_normal.clear();
    auto readVertices = [&vertices](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
    };
    auto readFaces = [&faces2vertices](stringstream &strs){
        string oneset,indstring;
        int ivalue;
        int nv = 0;
        vector<int>tempvlist;
        while (strs>>oneset){
            stringstream oneset_ss( oneset );
            getline( oneset_ss, indstring, '/' );
            stringstream indstring_ss(indstring);
            indstring_ss >> ivalue;
            nv++;
            tempvlist.push_back(ivalue-1);
        }
        assert(tempvlist.size()>=3);
        for(int i=2;i<tempvlist.size();++i){
            faces2vertices.push_back(tempvlist[0]);
            faces2vertices.push_back(tempvlist[i-1]);
            faces2vertices.push_back(tempvlist[i]);
        }
    };
    auto readVerticesNormal = [&vertices_normal](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices_normal.push_back(dvalue);}
    };

    string oneline;

    cout<<"reading: "<<filename<<endl;

    while( getline( fin, oneline ) ){
        stringstream strs( oneline );
        string prefix;

        strs >> prefix;

        if( prefix == "v"  ) { readVertices( strs ); continue; } // vertex
        if( prefix == "vt" ) {  continue; } // texture coordinate
        if( prefix == "vn" ) {  readVerticesNormal(strs);continue; } // vertex normal
        if( prefix == "vf" ) { /*readVerticesField( ss );*/ continue; } // tangent vector
        if( prefix == "f"  ) { readFaces( strs ); continue; } // face
        if( prefix[0] == '#' ) continue; // comment
        if( prefix == "o" ) continue; // object name
        if( prefix == "g" ) continue; // group name
        if( prefix == "s" ) continue; // smoothing group
        if( prefix == "mtllib" ) continue; // material library
        if( prefix == "usemtl" ) continue; // material
        if( prefix == "k" ) continue; // field degree
        if( prefix == "fs" ) continue; // field singularity
        if( prefix == "" ) continue; // empty string
        if( prefix == "c" ) continue;

        cout << "Error: not a valid curf file!" << endl;
        cout << "(Offending line: " << oneline << ")" << endl;
        return false;
    }


    fin.close();
    return true;



}


bool readObjFile_Line(string filename,vector<double>&vertices,vector<unsigned int>&edges2vertices){


    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        exit(-1213);
        return false;
    }

    vertices.clear();
    edges2vertices.clear();
    auto readVertices = [&vertices](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
    };
    auto readEdges = [&edges2vertices](stringstream &strs){
        string oneset,indstring;
        unsigned int ivalue,ivalue2;
        strs>>ivalue;
        while (strs>>ivalue2){
            edges2vertices.push_back(ivalue-1);
            edges2vertices.push_back(ivalue2-1);
            ivalue = ivalue2;
        }

        //        for(int i=0;i<3;++i){
        //            strs>>ivalue;faces2vertices.push_back(ivalue-1);
        //        }
    };


    string oneline;

    cout<<"reading: "<<filename<<endl;

    while( getline( fin, oneline ) ){
        stringstream strs( oneline );
        string prefix;

        strs >> prefix;

        if( prefix == "v"  ) { readVertices( strs ); continue; } // vertex
        if( prefix == "l"  ) { readEdges( strs ); continue; } // edges
        if( prefix == "vt" ) {  continue; } // texture coordinate
        if( prefix == "vn" ) {  continue; } // vertex normal
        if( prefix == "vf" ) { /*readVerticesField( ss );*/ continue; } // tangent vector
        if( prefix == "f"  ) {  continue; } // face
        if( prefix[0] == '#' ) continue; // comment
        if( prefix == "o" ) continue; // object name
        if( prefix == "g" ) continue; // group name
        if( prefix == "s" ) continue; // smoothing group
        if( prefix == "mtllib" ) continue; // material library
        if( prefix == "usemtl" ) continue; // material
        if( prefix == "k" ) continue; // field degree
        if( prefix == "fs" ) continue; // field singularity
        if( prefix == "" ) continue; // empty string
        if( prefix == "c" ) continue;

        cout << "Error: not a valid curf file!" << endl;
        cout << "(Offending line: " << oneline << ")" << endl;
        return false;
    }


    fin.close();
    return true;




}

bool writeObjFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices){
    filename = filename + ".obj";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output Obj file " << filename << endl;
        return false;
    }

    outer << setprecision(13);
    int n_vertices = vertices.size()/3;
    int n_faces = faces2vertices.size()/3;
    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }

    for(int i=0;i<n_faces;++i){
        auto p_fv = faces2vertices.data()+i*3;
        outer << "f " << p_fv[0]+1<< " "<< p_fv[1]+1 << " "<< p_fv[2]+1 << endl;
    }


    outer.close();
    cout<<"saving finish: "<<filename<<endl;
    return true;



}

bool writeObjFile_line(string filename, const vector<double>&vertices, const vector<unsigned int> &edge2vertices){

    filename = filename + ".obj";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output Obj file " << filename << endl;
        return false;
    }



    int n_vertices = vertices.size()/3;

    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }


    {
        int n_edges = edge2vertices.size()/2;
        for(int i=0;i<n_edges;++i){
            auto p_ev = edge2vertices.data()+i*2;
            outer << "l " << p_ev[0]+1<< " "<< p_ev[1]+1 << endl;
        }
    }


    outer.close();
    cout<<"saving finish: "<<filename<<endl;
    return true;
}


bool readSurfFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices,vector<double>&vertices_field){
    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    vertices.clear();
    faces2vertices.clear();
    vertices_field.clear();
    auto readVertices = [&vertices](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
    };
    auto readFaces = [&faces2vertices](stringstream &strs){
        int ivalue;
        for(int i=0;i<3;++i){
            strs>>ivalue;faces2vertices.push_back(ivalue-1);
        }
    };
    auto readVerticesField = [&vertices_field](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices_field.push_back(dvalue);}
    };

    string oneline;

    cout<<"reading: "<<filename<<endl;

    while( getline( fin, oneline ) ){
        stringstream strs( oneline );
        string prefix;

        strs >> prefix;

        if( prefix == "v"  ) { readVertices( strs ); continue; } // vertex
        if( prefix == "vt" ) {  continue; } // texture coordinate
        if( prefix == "vn" ) { /* readVerticesNormal(strs);*/continue; } // vertex normal
        if( prefix == "vf" ) { readVerticesField( strs ); continue; } // tangent vector
        if( prefix == "f"  ) { readFaces( strs ); continue; } // face
        if( prefix[0] == '#' ) continue; // comment
        if( prefix == "o" ) continue; // object name
        if( prefix == "g" ) continue; // group name
        if( prefix == "s" ) continue; // smoothing group
        if( prefix == "mtllib" ) continue; // material library
        if( prefix == "usemtl" ) continue; // material
        if( prefix == "k" ) continue; // field degree
        if( prefix == "fs" ) continue; // field singularity
        if( prefix == "" ) continue; // empty string
        if( prefix == "c" ) continue;

        cout << "Error: not a valid curf file!" << endl;
        cout << "(Offending line: " << oneline << ")" << endl;
        return false;
    }


    fin.close();
    return true;


}

bool writeSurfFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices,const vector<double>&vertices_field){



    filename = filename + ".surf";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output Surf file " << filename << endl;
        return false;
    }

    int n_vertices = vertices.size()/3;
    int n_faces = faces2vertices.size()/3;
    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << "v " << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }

    for(int i=0;i<n_faces;++i){
        auto p_fv = faces2vertices.data()+i*3;
        outer << "f " << p_fv[0]+1<< " "<< p_fv[1]+1 << " "<< p_fv[2]+1 << endl;
    }

    for(int i=0;i<n_vertices;++i){
        auto p_vvec = vertices_field.data()+i*3;
        outer << "vf " << p_vvec[0] << " "<< p_vvec[1] << " "<< p_vvec[2] << endl;
    }

    outer.close();
    cout<<"saving finish: "<<filename<<endl;

    return true;
}


bool readCurfFile(string filename,vector<double>&vertices,vector<unsigned int>&edges2vertices,vector<double>&vertices_field,vector<double>&vertices_tangent){

    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    vertices.clear();
    edges2vertices.clear();
    vertices_field.clear();
    auto readVertices = [&vertices](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
    };
    auto readEdges = [&edges2vertices](stringstream &strs){
        int ivalue;
        for(int i=0;i<2;++i){
            strs>>ivalue;edges2vertices.push_back(ivalue-1);
        }
    };
    auto readVerticesField = [&vertices_field](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices_field.push_back(dvalue);}
    };
    auto readVerticesTangent = [&vertices_tangent](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices_tangent.push_back(dvalue);}
    };


    string oneline;

    cout<<"reading: "<<filename<<endl;

    while( getline( fin, oneline ) ){
        stringstream strs( oneline );
        string prefix;

        strs >> prefix;

        if( prefix == "v"  ) { readVertices( strs ); continue; } // vertex
        if( prefix == "e" ) {  readEdges(strs);continue; } // texture coordinate
        if( prefix == "vn" ) {  readVerticesTangent(strs);continue; } // vertex normal(surface/volume)/tangent(curve)/LookAt Vector(Curve)
        if( prefix == "vf" ) { readVerticesField( strs ); continue; } // vertices field
        if( prefix[0] == '#' ) continue; // comment line

        cout << "Error: not a valid curf file!" << endl;
        cout << "(Offending line: " << oneline << ")" << endl;
        return false;
    }


    fin.close();
    return true;


}


bool writeCurfFile(string filename, const vector<double>&vertices, const vector<unsigned int>&edges2vertices , vector<double> &vertices_field, vector<double>&vertices_tangent){

    filename = filename + ".curf";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }
    int n_vertices = vertices.size()/3;
    int n_edges = edges2vertices.size()/2;

    auto p_vd = vertices.data();
    auto p_evd = edges2vertices.data();
    auto p_vvecd = vertices_field.data();
    auto p_vtd = vertices_tangent.data();
    for(int i =0; i<n_vertices;++i){
        auto p_v = p_vd+i*3;
        fout<<"v "<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<endl;

    }
    for(int i=0;i<n_edges;++i){
        auto p_ev = p_evd+i*2;
        fout<<"e "<<p_ev[0]+1<<' '<<p_ev[1]+1<<endl;
    }

    if(vertices_field.size()!=vertices.size()){
        cout<<"Output vertices field is invalid or empty!"<<endl;
    }else{
        for(int i =0; i<n_vertices;++i){
            auto p_vvec = p_vvecd+i*3;
            fout<<"vf "<<p_vvec[0]<<' '<<p_vvec[1]<<' '<<p_vvec[2]<<endl;
        }
    }

    if(vertices_tangent.size()!=vertices.size()){
        cout<<"Output vertices field is invalid or empty!"<<endl;
    }else{
        for(int i =0; i<n_vertices;++i){
            auto p_vt = p_vtd+i*3;
            fout<<"vn "<<p_vt[0]<<' '<<p_vt[1]<<' '<<p_vt[2]<<endl;
        }
    }
    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;


}


bool readVolfFile(string filename,vector<double>&vertices,vector<unsigned int>&tets2vertices,vector<double> &vertices_normal,vector<double> &vertices_field){

    ifstream fin(filename.data());
    if(fin.fail()){
        cout<<"Fail to open input file: "<<filename<<endl;
        return false;
    }

    vertices.clear();
    tets2vertices.clear();
    vertices_field.clear();
    auto readVertices = [&vertices](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices.push_back(dvalue);}
    };
    auto readTets = [&tets2vertices](stringstream &strs){
        int ivalue;
        for(int i=0;i<4;++i){
            strs>>ivalue;tets2vertices.push_back(ivalue-1);
        }
    };
    auto readVerticesField = [&vertices_field](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices_field.push_back(dvalue);}
    };
    auto readVerticesNormal = [&vertices_normal](stringstream &strs){
        double dvalue;
        for(int i=0;i<3;++i){strs>>dvalue;vertices_normal.push_back(dvalue);}
    };

    string oneline;

    cout<<"reading: "<<filename<<endl;

    while( getline( fin, oneline ) ){
        stringstream strs( oneline );
        string prefix;

        strs >> prefix;

        if( prefix == "v"  ) { readVertices( strs ); continue; } // vertex
        if( prefix == "vn" ) {  readVerticesNormal(strs);continue; } // vertex normal
        if( prefix == "vf" ) { readVerticesField( strs ); continue; } // tangent vector
        if( prefix == "tet"  ) { readTets( strs ); continue; } // face
        if( prefix[0] == '#' ) continue; // comment

        cout << "Error: not a valid curf file!" << endl;
        cout << "(Offending line: " << oneline << ")" << endl;
        return false;
    }


    fin.close();
    return true;

}


bool writeVolfFile(string filename, const vector<double>&vertices, const vector<unsigned int>&tets2vertices, vector<double> &vertices_normal, vector<double> &vertices_field){
    filename = filename + ".volf";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }
    int n_vertices = vertices.size()/3;
    int n_tets = tets2vertices.size()/4;

    auto p_vd = vertices.data();
    auto p_tvd = tets2vertices.data();
    auto p_vvecd = vertices_field.data();
    auto p_vnd = vertices_normal.data();
    for(int i =0; i<n_vertices;++i){
        auto p_v = p_vd+i*3;
        fout<<"v "<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<endl;

    }
    for(int i=0;i<n_tets;++i){
        auto p_tv = p_tvd+i*4;
        fout<<"tet "<<p_tv[0]+1<<' '<<p_tv[1]+1<<' '<<p_tv[2]+1<<' '<<p_tv[3]+1<<endl;
    }

    if(vertices_normal.size()/3!=n_vertices){
        cout<<"Output vertices normal is invalid or empty!"<<endl;
    }else{
        for(int i =0; i<n_vertices;++i){
            auto p_vn = p_vnd+i*3;
            fout<<"vn "<<p_vn[0]<<' '<<p_vn[1]<<' '<<p_vn[2]<<endl;
        }
    }
    if(vertices_field.size()!=vertices.size()){
        cout<<"Output vertices field is invalid or empty!"<<endl;
    }else{
        for(int i =0; i<n_vertices;++i){
            auto p_vvec = p_vvecd+i*3;
            fout<<"vf "<<p_vvec[0]<<' '<<p_vvec[1]<<' '<<p_vvec[2]<<endl;
        }
    }


    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;
}
bool readContourEdgeTxtFile(string filename, vector<int>&edges2vertices){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    edges2vertices.resize(nnum*2);
    for(int i=0;i<nnum*2;++i)reader>>edges2vertices[i];
    reader.close();
    return true;
}
bool writeContourEdgeTxtFile(string filename, const vector<unsigned int>&edges2vertices){
    filename = filename + "_ContuorEdges.txt";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }

    int numofE = edges2vertices.size()/2;
    fout<<numofE<<endl;
    for(int i=0;i<numofE;++i){
        fout<<edges2vertices[i*2]<<' '<<edges2vertices[i*2+1]<<endl;
    }
    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;
}
bool writeVecFile(string filename, const vector<int> &vec){
    filename = filename + "_vec.txt";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }

    int numof = vec.size();
    fout<<numof<<endl;
    if(numof==0)return true;
    for(auto a:vec)fout<<a<<endl;
    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;
}
bool readVecFile(string filename, vector<int> &vec){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    vec.resize(nnum);
    for(int i=0;i<nnum;++i)reader>>vec[i];
    reader.close();
    return true;
}
bool readVVecFile(string filename, vector<vector<int>> &vvec){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    vector<int>vvnum(nnum,0);
    for(int i=0;i<nnum;++i)reader>>vvnum[i];

    vvec.resize(nnum);
    for(int i=0;i<nnum;++i){
        vvec[i].resize(vvnum[i]);
        for(int j=0;j<vvnum[i];++j)reader>>vvec[i][j];
    }
    reader.close();
    return true;
}

bool writeVVecFile(string filename, const vector< vector<int> > &vvec){
    filename = filename + "_vvec.txt";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }

    int numof = vvec.size();

    fout<<numof<<endl;

    if(numof==0)return true;
    for(int i=0;i<numof-1;i++){
        fout<<vvec[i].size()<<' ';
    }
    fout<<vvec[numof-1].size()<<endl;

    for(auto &a:vvec){
        for(int i=0;i<a.size()-1;i++){
            fout<<a[i]<<' ';
        }
        fout<<a[a.size()-1]<<endl;
    }
    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;
}
bool readVVecFile(string filename, vector<vector<double>> &vvec){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    vector<int>vvnum(nnum,0);
    for(int i=0;i<nnum;++i)reader>>vvnum[i];

    vvec.resize(nnum);
    for(int i=0;i<nnum;++i){
        vvec[i].resize(vvnum[i]);
        for(int j=0;j<vvnum[i];++j)reader>>vvec[i][j];
    }
    reader.close();
    return true;
}

bool writeVVecFile(string filename, const vector< vector<double> > &vvec){
    filename = filename + "_vvec.txt";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }

    int numof = vvec.size();

    fout<<numof<<endl;

    if(numof==0)return true;
    for(int i=0;i<numof-1;i++){
        fout<<vvec[i].size()<<' ';
    }
    fout<<vvec[numof-1].size()<<endl;

    for(auto &a:vvec){
        for(int i=0;i<a.size()-1;i++){
            fout<<a[i]<<' ';
        }
        fout<<a[a.size()-1]<<endl;
    }
    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;
}
bool writeVecFile(string filename, const vector<float>&vec){
    filename = filename + "_vec.txt";
    ofstream fout(filename.data());
    if(fout.fail()){
        cout<<"Fail to create output file: "<<filename<<endl;
        return false;
    }

    int numof = vec.size();
    fout<<numof<<endl;
    for(auto a:vec)fout<<a<<endl;
    fout.close();

    cout<<"Write: "<<filename<<endl;
    return true;
}
bool readVecFile(string filename,vector<float>&vec){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    vec.resize(nnum);
    for(int i=0;i<nnum;++i)reader>>vec[i];

    reader.close();
    return true;

}
bool readVecFile(string filename,vector<double>&vec){
    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }
    int nnum = 0;
    reader>>nnum;
    vec.resize(nnum);
    for(int i=0;i<nnum;++i)reader>>vec[i];

    reader.close();
    return true;

}
bool writeSufFile(string filename, const vector<double>&vertices, const vector<unsigned int>&faces2vertices, const vector<int>&facesMat, const vector<int> &CtrEdges){

    filename = filename + ".suf";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output suf file " << filename << endl;
        return false;
    }
    int n_vertices = vertices.size()/3;
    int n_faces = faces2vertices.size()/3;


    outer<<n_vertices<<' '<<n_faces<<endl;
    for(int i=0;i<n_vertices;++i){
        auto p_v = vertices.data()+i*3;
        outer << p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }
    for(int i=0;i<n_faces;++i){
        auto p_fv = faces2vertices.data()+i*3;
        auto p_fm = facesMat.data()+i*2;
        outer<< p_fv[0]<< " "<< p_fv[1] << " "<< p_fv[2] << " "<< p_fm[0]<< " "<<p_fm[1] <<endl;
    }

    outer<<CtrEdges.size()/2<<endl;
    for(int i=0;i<CtrEdges.size()/2;++i){
        outer<<CtrEdges[i*2]<<' '<<CtrEdges[i*2+1]<<endl;
    }
    outer.close();
    return 1;
}
bool writeCtrGraphFile(string filename,const vector<float>&vertices,const vector<vector<int>>&edge2vertices,const vector<vector<int>>&edgeMat, const vector<vector<float>>&planepara){
    filename = filename + ".CtrGraph";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output CtrGraph file " << filename << endl;
        return false;
    }

    int n_vertices = vertices.size()/3;
    int n_plane = edge2vertices.size();
    outer<<"n "<<n_vertices<<endl;
    auto p_vd = vertices.data();
    for(int i=0;i<n_vertices;++i){
        auto p_v = p_vd+ i*3;
        outer << "v "<<p_v[0] << " "<< p_v[1] << " "<< p_v[2] << endl;
    }

    outer<<"n "<<n_plane<<endl;
    for(int i=0;i<n_plane;++i){
        outer<<"p ";
        for(int j=0;j<4;++j)outer<<planepara[i][j]<<' ';

        auto &p_ev = edge2vertices[i];
        auto &p_em = edgeMat[i];
        int ne = p_ev.size()/2;
        outer<<ne<<endl;

        for(int j=0;j<ne;++j){
            auto ind = j*2;
            outer<<"e "<<p_ev[ind]<<' '<<p_ev[ind+1]<<' '<<p_em[ind]<<' '<<p_em[ind+1]<<endl;
        }
    }

    outer.close();


    cout<<"finish: "<<filename<<endl;
    return true;





}


bool writeCurNetFile(string filename, const vector<double> &vertices, const vector<vector<int> > &edge2vertices, const vector<vector<int> > &edgeMat, const vector<vector<double> > &planepara, const vector<double>&verticesNor){

    filename = filename + ".CurNet";
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output CurNet file " << filename << endl;
        return false;
    }

    int n_vertices = vertices.size()/3;
    int n_plane = edge2vertices.size();
    outer<<"n "<<n_vertices<<endl;
    auto p_vd = vertices.data();
    auto p_vnd = verticesNor.data();
    for(int i=0;i<n_vertices;++i){
        auto p_v = p_vd+ i*3;
        auto p_vn = p_vnd+i*3;
        outer << "v "<<p_v[0] << " "<< p_v[1] << " "<< p_v[2]  << " " <<p_vn[0] << " "<< p_vn[1] << " "<< p_vn[2] << endl;
    }

    outer<<"n "<<n_plane<<endl;
    for(int i=0;i<n_plane;++i){
        outer<<"p ";
        for(int j=0;j<4;++j)outer<<planepara[i][j]<<' ';

        auto &p_ev = edge2vertices[i];
        auto &p_em = edgeMat[i];
        int ne = p_ev.size()/2;
        outer<<ne<<endl;

        for(int j=0;j<ne;++j){
            auto ind = j*2;
            outer<<"e "<<p_ev[ind]<<' '<<p_ev[ind+1]<<' '<<p_em[ind]<<' '<<p_em[ind+1]<<endl;
        }
    }

    outer.close();

    return true;




}


bool readXYZ(string filename, vector<double>&v){

    ifstream reader(filename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }else {
        cout << "Reading: "<<filename<<endl;
    }
    v.clear();
    double val;
    while(!reader.eof()){
        reader>>val;
        v.push_back(val);
    }
    reader.close();
    return true;
}

bool readXYZnormal(string filename, vector<double>&v, vector<double>&vn){

    ifstream reader(filename.data(), ofstream::in);
    if (!reader.good()) {
        cout << "Can not open the file " << filename << endl;
        return false;
    }else {
        cout << "Reading: "<<filename<<endl;
    }
    v.clear();
    vn.clear();
    double val;
    while(!reader.eof()){
        for(int i=0;i<3;++i){reader>>val;v.push_back(val);}
        for(int i=0;i<3;++i){reader>>val;vn.push_back(val);}
    }
    reader.close();

    return true;
}

bool writeXYZ(string filename, vector<double>&v){

    int npt = v.size()/3;
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output file " << filename << endl;
        return false;
    }

    for(int i=0;i<npt;++i){
        auto p_v = v.data()+i*3;
        outer<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<endl;
    }
    outer.close();
    return true;

}

bool writeXYZnormal(string filename, vector<double>&v, vector<double>&vn){

    int npt = v.size()/3;
    ofstream outer(filename.data(), ofstream::out);
    if (!outer.good()) {
        cout << "Can not create output file " << filename << endl;
        return false;
    }

    for(int i=0;i<npt;++i){
        auto p_v = v.data()+i*3;
        auto p_vn = vn.data()+i*3;
        outer<<p_v[0]<<' '<<p_v[1]<<' '<<p_v[2]<<' ';
        outer<<p_vn[0]<<' '<<p_vn[1]<<' '<<p_vn[2]<<endl;
    }
    outer.close();
    return true;

}

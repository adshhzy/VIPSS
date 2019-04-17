#ifndef READERS_H
#define READERS_H


#include<vector>
#include<string>
using namespace std;

bool readOffFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices);


bool writeOffFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices);

bool writePLYFile(string filename, const vector<double>&vertices, const vector<unsigned int>&faces2vertices,
                  const vector<double> &vertices_normal, const vector<unsigned char>&vertices_color);

bool readPLYFile(string filename,  vector<double>&vertices, vector<double> &vertices_normal);

bool readObjFile(string filename, vector<double>&vertices, vector<unsigned int>&faces2vertices, vector<double> &vertices_normal);
bool readObjFile_Line(string filename,vector<double>&vertices,vector<unsigned int>&edges2vertices);

bool writeObjFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices);
bool writeObjFile_line(string filename,const vector<double>&vertices,const vector<unsigned int>&edge2vertices);

bool readSurfFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices,vector<double>&vertices_field);

bool writeSurfFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices,const vector<double>&vertices_field);


bool readCurfFile(string filename, vector<double>&vertices, vector<unsigned int>&edges2vertices, vector<double>&vertices_field, vector<double> &vertices_tangent);


bool writeCurfFile(string filename, const vector<double>&vertices, const vector<unsigned int>&edges2vertices, vector<double>&vertices_field , vector<double> &vertices_tangent);


bool readVolfFile(string filename, vector<double>&vertices, vector<unsigned int>&tets2vertices, vector<double> &vertices_normal, vector<double> &vertices_field);


bool writeVolfFile(string filename, const vector<double>&vertices, const vector<unsigned int>&tets2vertices, vector<double> &vertices_normal, vector<double> &vertices_field);

bool readContourEdgeTxtFile(string filename, vector<int>&edges2vertices);
bool writeContourEdgeTxtFile(string filename, const vector<unsigned int>&edges2vertices);



bool writeVecFile(string filename, const vector<int>&vec);
bool readVecFile(string filename,vector<int>&vec);



bool readVVecFile(string filename, vector<vector<int>> &vvec);
bool writeVVecFile(string filename, const vector<vector<int>> &vvec);

bool writeVecFile(string filename, const vector<float>&vec);
bool readVecFile(string filename,vector<float>&vec);


bool readVecFile(string filename,vector<double>&vec);

bool readVVecFile(string filename, vector<vector<double>> &vvec);
bool writeVVecFile(string filename, const vector<vector<double>> &vvec);

bool writeSufFile(string filename, const vector<double>&vertices, const vector<unsigned int>&faces2vertices, const vector<int>&facesMat, const vector<int> &CtrEdges);
bool writeCtrGraphFile(string filename, const vector<float> &vertices, const vector<vector<int>>&edge2vertices, const vector<vector<int>>&edgeMat, const vector<vector<float>>&planepara);

bool writeCurNetFile(string filename, const vector<double> &vertices, const vector<vector<int>>&edge2vertices, const vector<vector<int>>&edgeMat, const vector<vector<double>>&planepara, const vector<double>&verticesNor);

bool readXYZ(string filename, vector<double>&v);
bool readXYZnormal(string filename, vector<double>&v, vector<double>&vn);
bool writeXYZ(string filename, vector<double>&v);
bool writeXYZnormal(string filename, vector<double>&v, vector<double>&vn);

#endif // READERS_H

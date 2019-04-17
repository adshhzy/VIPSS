#ifndef _CTRSUFMANAGER_H_
#define _CTRSUFMANAGER_H_

#include "../../src/utils/Ctr2Suf/config.h"
#include "../../src/utils/Ctr2Suf/ContourHandler/ContourHandler.h"
#include "../../src/utils/Ctr2Suf/SpacePartitioner/SpacePartitioner.h"
#include "../../src/utils/Ctr2Suf/MAGenerator/MAGenerator.h"
#include "../../src/utils/Ctr2Suf/Util/MapArraySR.h"
#include "../../src/utils/Ctr2Suf/SubspaceContour/sscontour.h"
class SSConotur;

#include "../../src/utils/Ctr2Suf/Projection/Projector.h"
#include "../../src/utils/Ctr2Suf/Projection/ProjSplitter.h"

#include "../../src/utils/Ctr2Suf/RegionHandler/regionHandler.h"

#include "../../src/utils/Ctr2Suf/FaceGenerator/FaceGenerator.h"

#include "../../src/utils/Ctr2Suf/LuMesh/LuMesh.h"

/**
* Function 
* @param 
* @param
* @param
* @param
* @param
*/
//extern const float cols[ 18 ][ 3 ];
//const float cols[18][3] =	{ {79, 184, 196},{255, 0, 0}, {0, 255, 0}, {0, 255, 255}, 	{0, 0, 255}, {255, 0, 255}, 
//{128, 	255, 0}, {255, 128, 0}, {0, 128, 255}, {128, 0, 255}, {255, 128, 128}, {128, 128, 255}, 
//{243, 117, 76}, {127, 145, 110}, {128, 23, 82}, {249, 249, 19},//gold
//{102,2, 102}, {240, 240, 240}};

//const float DIM = 1.5;
//const float PROCSIZE = 1000;
class Ctr2SufManager
{
	

	float height;
	float width;
	float nearplane;
	float farplane;

public:
	LuMesh* mesh;

	void setGLParam(float width, float height, float nearplane, float farplane)
	{
		this->width = width;
		this->height = height;
		this->nearplane = nearplane;
		this->farplane = farplane;
		if( mesh!=NULL)
			mesh->setGLParam( width, height, nearplane, farplane );
	}

	Ctr2SufManager();
	~Ctr2SufManager();

	//for show
	float** showctrvers;		//contour vertices to show
	int showplanenum;
	//int** showctrvernum;
	int* showctredgenum;
	int** showctredges;

	bool inresize;
	void resize();			
	void render();
	void clearContour();

	//process all of them
	//void getMeshFromContour();

	//!READ IN FILES AND PROCESS!
	floatvector param;				//the parameters of the planes
    vector<floatvector> ctrvers;	//vertices of the contours
    vector<intvector> ctredges;		//edges of the contours
	float bbox[ 6 ];			//bounding box of the contours
	bool bboxset;				//if bounding box is set or not
	void readContourFromMemory();
	void readContour(const int filenum, char** fnames);
	void addContour(const int filenum, char** fnames);
	bool hasInterPtInfo;
	void addContour2(const int filenum, char** fnames);
	void readContour2(const int filenum, char** fnames);
	
	void readContour3(const int filenum, char** fnames);

	void clearOldContour();
	bool showContour;
	void renderContour();

	//new for 1 set of radial contours, and 1 set of random contours
	int radialContourNum;
	
	//!PREPROCESS THE DATA!
	bool uniform;			//if uniform the contours or not
	int expectlen;			//if uniform, what is the expect length
	float center[ 3 ];
	float unitlen;			//the biggest length should be 2.
    float pbbox[ 6 ];		//bouding box of the processed contour
	int planenum;			//number of the planes
	float* pparam;			//parameters of the plane after process
	float** pctrvers;		//contour vertices after process
	int* pctrvernum;		//contour vertices number on each plane
	int** pctredges;		//contour edges after process
	int* pctredgenum;		//number of contour edges on each plane
	bool ispreproced;

	//ver2planelist
	int** preProcData(bool ISCOMMONLINE = false);
	void writePreProcData();
	int** ver2planelistpos;
	vector<intvector> ver2planelist;

	void clearPData();		//clear all the processed data

	//for common edge case
	bool isComnCase;										//true- is cutting plane sharing one common line case
	float comndir[ 3 ];
	float comnpt[3 ];
	int comnedgei;
	int comnveri[ 2 ];
	void toggleIsComnCase( bool isComnCase ){ this->isComnCase = isComnCase; }
	
	//!PARTITION SPACE!
	int ssvernum;
	float* ssver;
	int ssedgenum;
	int* ssedge;
	int ssfacenum;
	int* ssfaceedgenum;
	int** ssface;
	int* ssface_planeindex;
	int ssspacenum;
	int* ssspacefacenum;
	int** ssspace;
	int** ssspace_planeside;
	float enlargeratio[ 3 ];	
	void partition();
	void clearPartition();
	bool showPartition;
	void renderPartition();
	void writePartitionOut(const char* fname);

	int* ctrfvernum;
	float** ctrfverpos;
	int** ctrfvertype;
	int** ctrfverval;
	int* ctrfedgenum;
	int** ctrfedge;
	int** ctrfedgetype;
	int** ctrfedgeval;
	int** ctrfedgeancestor;
	void putContourIntoFace(int** ctrverMark);		//only used in common line case
	void clearFaceContour();
	bool gatherSubspaceCtr(vector<SSPCTRVERVEC>& sspctrver_vec,vector<SSPCTREDGEVEC>& sspctredge_vec, int spacei );
	/*void writeGatheredCtr_DB(vector<SSPCTRVERVEC>& sspctrver_vec,	vector<SSPCTREDGEVEC>& sspctredge_vec, int spacei);*/
	//void divideOneFaceContour(
	//	int spacei,
	//	int facei,
	//	//contour
	//	SSPCTRVERVEC& sspctrver_vec,	//temporary result
	//	SSPCTRVERVEC& sspctredge_vec,
	//	//ma
	//	int subvernum,float* subver,int subedgenum,int* subedge,int subfacenum,int* subfaceedgenum,int** subface,
	//	float* subparam,int* subver2wver,int* subedge2wedg,int majptnum,float* majpt,int maseamnum,
	//	int* maseam,MapArraySR& doubleface2sheet,int* seamonsheetnum,int** seamonsheet,int* ver2jpt
	//	);
	void divideFaceContourByProjMA(
		int spacei,
		//contour
	//	SUBSPACECTRVER**& sspctrver,int& sspctrvernum,	//final result
	//	SUBSPACECTREDGE**& sspctredge,int& sspctredgenum,
		vector<SSPCTRVERVEC>& sspctrver_vec,	//temporary result
		vector<SSPCTREDGEVEC>& sspctredge_vec,
		//ma
		int subvernum,float* subver,int subedgenum,int* subedge,int subfacenum,int* subfaceedgenum,int** subface,
		float* subparam,int* subver2wver,int* subedge2wedge,int majptnum,float* majpt,int maseamnum,
		int* maseam,MapArraySR& doubleface2sheet,int* seamonsheetnum,int** seamonsheet,float* sheettab,int* ver2jpt);

	//!CTR TO SURFACE PROCESS!
	void generateMA(int spaci,
		//subspace info
		int& subvernum, float*& subver, int& subedgenum,int*& subedge,int& subfacenum,
		int*& subfaceedgenum, int**& subface, float*& subparam, int*& subver2wver,int*& subedge2wedge,
		//ma info
		int& majptnum, float*& majpt, int& maseamnum, int*& maseam, 
		MapArraySR& doubleface2sheet, int*& seamonsheetnum, int**& seamonsheet,  float*& sheettab,
		int*& ver2jpt);			//generate MA

	//when stitching along vertex, edge, and face, need these registration info
	//one for each subspace vertex
	//in each vector, saves sequence of (spacei, subspace mesh vertex index)
	intvector* sverreg;
	//one for each edge
	//double array, each array for one vertex on the edge
	vector<intvector>* sedgereg;
	//mark the subedge is in the final mesh or not
	intvector* sedgesubedgemark;
	//one for one subsapce face
	//one intset saves the ancestor contour edge numbers on this face
	intset* sfacectrei;
	//two for one subspace face, save the space index share this face
	int* sfacespaci;
	//one for each face
	//double array, first dimension means different contour edge, two for one contour edge
	//second dimension, saves all the faces incident with the contour edge
	vector<intvector>* sfaceregface;
	//one vector for one face
	//in each vector, saves the list of pair of (v1, v2) from (spaci1, spaci2) in sfacespaci
	intvector* sfaceregver;
	//in each vector, saves the list of one edge (ind1 ind2) in sfacespaci
	//ind1 ind2 are the two positions of two vertices in sfaceregver
	intvector* sfaceregedgever;

	//generate mesh in one subspace
	void ctr2sufSubspaceProc(int spaci, floatvector& meshVer, intvector& meshEdge, intvector& meshFace );	//process each subspace
//	void ctr2sufSubspaceProc(int spaci, floatvector& meshVer, intvector& meshFace );	//process each subspace

	//result mesh
	floatvector mver;		//three is a group, the position of the vertex
	//intvector medge;
	intvector mface;		//five is a group, v1, v2, v3, and m1 and m2
	intvector ctrmedge;	//two is one edge, mark the contour edge

	//stitch the meshes from subspaces
	void sortFace(intvector& meshFace,intvector& regface,int*& sregface);
	void splitCtrEdgeOnFace(
		floatvector*& subMeshVer,intvector*& subMeshFace);
	void stitchMesh(floatvector*& subMeshVer,intvector*& subMeshFace);

	void ctr2sufProc();			//generate mesh from contour
	
	//render
	void renderSufMesh();

	//debug
	void writeSheetTab(		float* masheettab,		int masheetnum,		int spaci		);


	///for debug reason
	vector<intset> ancestorlist;
	  vector<intset> ancestorlist_split;

	  bool dbOneSpaceMode;
	  int dbCurSpace;
	  int dbSpaceNum;
	  void setMode( bool onespacmode ){ dbOneSpaceMode = onespacmode; }
	  void setCurSpace(int curspace ){ dbCurSpace = curspace; }
	  int getSpaceNum(){return dbSpaceNum;}
};




#endif
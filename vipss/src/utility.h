#ifndef PREDEFINE_H
#define PREDEFINE_H

//#define USINGQVECTOR
//#ifdef USINGQVECTOR
//#define Vector QVector
//#else
//#include<vector>
//#define Vector vector
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned long ulong;
typedef unsigned long long ulonglong;

#define my_PI 3.141592653589793238463

#include <climits>
#include<limits>
#include<math.h>
#include <stdio.h>
#include <string>
#include <vector>
#include "dirent.h"
namespace MyUtility {

#define UINTFLAG  std::numeric_limits<unsigned int>::max()

//double randomdouble() {return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);}
//double randomdouble(double be,double ed) {return be + randomdouble()*(ed-be);	}
/****************************************************************/

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

inline void GetExtName (const std::string& fullfilename,std::string &extname) {
    int pos;
    pos = fullfilename.find_last_of('.');
    extname = fullfilename.substr(pos);
}

/*************************************************************/
template <class T>
inline T dot(const T *e1,const T *e2,const int dim = 3){
    T d = 0.0;
    for(int i =0;i<dim;++i)d+=e1[i]*e2[i];
    return d;
}
template <class T>
inline T normalize(T *face_normal,const int dim = 3){
    T len = sqrt(dot(face_normal,face_normal,dim));
    for(int i =0;i<dim;++i)face_normal[i] /= len;
    return len;
}
template <class T>
inline void inversevec(T *p_s, T *p_d, const int dim =3){
    for(int i =0;i<dim;++i)p_d[i] = -p_s[i];
}
template <class T>
inline void product(const T a,const T *vin, T *vout,int dim = 3){
    for(int i =0;i<dim;++i)vout[i]=a*vin[i];
}
template <class T>
inline void cross(const T *e1,const T *e2,T *pn){
    pn[0] = e1[1] * e2[2] - e2[1] * e1[2];
    pn[1] = e1[2] * e2[0] - e2[2] * e1[0];
    pn[2] = e1[0] * e2[1] - e2[0] * e1[1];
}
//    static T dot(T *e1,T *e2){
//        return e1[0]*e2[0]+e1[1]*e2[1]+e1[2]*e2[2];
//    }
template <class T>
inline void negVec(T *p_v1, int vdim=3){
    for(int i = 0;i<vdim;++i)p_v1[i] = -p_v1[i];
}

template <class T>
inline void add(const T *p_v1,const T *p_v2, T *p_vc, int vdim=3){
    for(int i = 0;i<vdim;++i)p_vc[i] = (p_v1[i] + p_v2[i]);
}
template <class T>
inline void minusVec(const T *v1,const T *v2,T *e,const int dim = 3){
    for(int i =0;i<dim;++i)e[i]= v1[i]-v2[i];
}
template <class T>
inline T normVec(const T *e){return sqrt(dot(e,e));}
template <class T>
inline T len(const T *e){
    return dot(e,e);
}
template <class T>
inline void copyVec(T *vo, T *vd,int dim = 3){
    for(int i =0;i<dim;i++)vd[i] = vo[i];
}
template <class T>
inline T cosine(const T *e1,const T *e2){
    T a = dot(e1,e2)/normVec(e1)/normVec(e2);
    if(a>1)a=1;
    if(a<-1)a=-1;
    return a;
}
template <class T>
inline T angleNor(const T *e1,const T *e2){
    return acos(cosine(e1,e2));
}
template <class T>
inline T _VerticesDistance(const T *p_v1,const T *p_v2,const int vdim=3){
    T dist = 0.0;
    for(int i = 0;i<vdim;++i)dist += (p_v1[i] - p_v2[i] )*(p_v1[i] - p_v2[i] );
    return sqrt(dist);
}
template <class T>
inline T vecSquareDist(const T *p_v1,const T *p_v2, const int vdim = 3){
    T dist = 0.0;
    for(int i = 0;i<vdim;++i)dist += (p_v1[i] - p_v2[i] )*(p_v1[i] - p_v2[i] );
    return dist;
}
template <class T>
inline void _VerticesMidpoint(const T *p_v1,const T *p_v2, T *p_vc, const int vdim = 3){
    for(int i = 0;i<vdim;++i)p_vc[i] = (p_v1[i] + p_v2[i]) /2;
}
template <class T>
inline void _TriangleMidpoint(const T *p_v1,const T *p_v2,const  T *p_v3,T *p_vc,const int vdim = 3){
    for(int i = 0;i<vdim;++i)p_vc[i] = (p_v1[i] + p_v2[i] + p_v3[i]) /3;
}
template <class T>
T inline _TriangleArea(const T *v1,const T *v2,const T *v3){
    T e1[3], e2[3], p_n[3];
    for (uint j = 0; j < 3; ++j){
        e1[j] = v2[j] - v1[j];
        e2[j] = v3[j] - v2[j];
    }
    p_n[0] = e1[1] * e2[2] - e2[1] * e1[2];
    p_n[1] = e1[2] * e2[0] - e2[2] * e1[0];
    p_n[2] = e1[0] * e2[1] - e2[0] * e1[1];
    return(0.5*sqrt(p_n[0] * p_n[0] + p_n[1] * p_n[1] + p_n[2] * p_n[2]));
}
template <class T>
void inline _TriangleCircumcenter(const T *pv1,const T *pv2,const T *pv3,T *cc,T &r){

    T n[3],e1[3],e2[3];
    for(int i = 0;i<3;++i){e1[i]=pv1[i]-pv2[i];e2[i]=pv1[i]-pv3[i];}
    cross(e1,e2,n);normalize(n);
    T c1[3],c2[3];
    _VerticesMidpoint(pv1,pv2,c1,3);
    _VerticesMidpoint(pv1,pv3,c2,3);
    T n1[3],n2[3];
    cross(e1,n,n1);cross(e2,n,n2);
    T t = ((c2[0]-c1[0])*n2[1] - (c2[1]-c1[1])*n2[0])/(n1[0]*n2[1] - n1[1]*n2[0]);
    for(int i = 0;i<3;++i)cc[i] = c1[i] + n1[i]*t;
    r = _VerticesDistance(pv1,cc,3);
    //cout<<VerticesDistance(pv1,cc,3)<<' '<<VerticesDistance(pv2,cc,3)<<' '<<VerticesDistance(pv3,cc,3)<<endl;
}
template <class T>
T inline _TriangleLeastAngle(const T *pv1, const T *pv2, const T *pv3,const int dim){
    T leastangle;
    T vdis1,vdis2,vdis3;
    T dot1,dot2,dot3;
    T e1[3],e2[3],e3[3];
    vdis1 = _VerticesDistance(pv1,pv2,dim);
    vdis2 = _VerticesDistance(pv1,pv3,dim);
    vdis3 = _VerticesDistance(pv2,pv3,dim);
    minusVec(pv1,pv2,e1,dim);minusVec(pv1,pv3,e2,dim);minusVec(pv2,pv3,e3,dim);
    dot1 = dot(e1,e2,dim);dot2 = -dot(e1,e3,dim);dot3 = dot(e2,e3,dim);
    leastangle = fmax(fmax(dot1/vdis1/vdis2,dot2/vdis1/vdis3),dot3/vdis2/vdis3);
    //        angle1 = acosf(dot1/vdis1/vdis2);
    //        angle2 = acosf(dot2/vdis1/vdis3);
    //        angle3 = acosf(dot3/vdis2/vdis3);
    //        return min(min(angle1,angle2),angle3);
    return acosf(leastangle);
}

template <class T>
inline T projectVectorUNor(const T *proj, const T *normal,T *vec){
    T dis = dot(proj,normal);
    T tmp[3];
    product(dis,normal,tmp);
    minusVec(proj,tmp,vec);
    return (dis);
}
template <class T>
inline T projectVectorNor(const T *proj, const T *normal,T *vec){
    T dis = projectVectorUNor(proj,normal,vec);
    normalize(vec);
    return (dis);
}
template <class T>
inline T projectVectorLeng(const T *proj, const T *normal,T *vec){
    T dis = projectVectorNor(proj,normal,vec);
    T lengg = normVec(proj);
    product(lengg,vec,vec);
    return (dis);
}
template <class T>
inline void weightedAddVec(const T weight, const T *orivec,T *desvec){
    for(int i =0;i<3;++i)desvec[i] += weight*orivec[i];
}
template <class T>
inline void weightedAddVec(const T weight1,const T weight2, const T *orivec1,const T *orivec2,T *desvec){
    for(int i =0;i<3;++i)desvec[i] = weight1*orivec1[i] + weight2*orivec2[i];
}
template <class T>
inline T computePoint2LineDistance(const T *p, const T *lst, const T *ldir){
    T vec[3];
    minusVec(p,lst,vec);
    T d2st =dot(vec,ldir);
    product(d2st,ldir,vec);
    add(lst,vec,vec);
    return _VerticesDistance(p,vec);

}
template <class T>
inline T computePoint2LineDistance(const T *p, const T *lst, const T *ldir, T& d2st){
    T vec[3];
    minusVec(p,lst,vec);
    d2st =dot(vec,ldir);
    product(d2st,ldir,vec);
    add(lst,vec,vec);
    return _VerticesDistance(p,vec);

}

template <class T>
inline T relativeAngle(const T *vec,const T *ref, const T *normal){
    T tmpcross[3];
    T dd = angleNor(vec,ref);
    cross(ref,vec,tmpcross);
    if(dot(tmpcross,normal)<0)dd=-dd;
    return dd;

}
template <class T>
inline T point2TriSquareDistance(const T *P, const T *v1, const T *v2, const T *v3,T *cp){
    T vec[3],E0[3],E1[3];
    const T *B = v1;
    minusVec(v2,v1,E0);
    minusVec(v3,v1,E1);
    T a = dot(E0,E0);
    T b = dot(E0,E1);
    T c = dot(E1,E1);
    minusVec(B,P,vec);
    T d = dot(E0,vec);
    T e = dot(E1,vec);
    T f = dot(vec,vec);

    T det = a*c-b*b, s = b*e-c*d, t = b*d-a*e;

    auto f0 = [&s,&t,&det](){T invdet = 1/det; s*=invdet;t*=invdet;};
    auto f1 = [&s,&t,&det,&a,&b,&c,&d,&e,&f](){
        T numer = c+e-b-d;
        if ( numer <= 0 ){
            s = 0;
        }else{
            T denom = a-2*b+c;
            s = ( numer >= denom ? 1 : numer/denom );
        }
        t = 1-s;
    };
    auto f2 = [&s,&t,&det,&a,&b,&c,&d,&e,&f](){
        T tmp0 = b+d;
        T tmp1 = c+e;
        if ( tmp1 > tmp0 ){
            T numer = tmp1 - tmp0;
            T denom = a-2*b+c;
            s = ( numer >= denom ? 1 : numer/denom );
            t = 1-s;
        }else {
            s = 0;
            t = ( tmp1 <= 0 ? 1 : ( e >= 0 ? 0 : -e/c ) );
        }
    };
    auto f3 = [&s,&t,&det,&e,&c](){
        s = 0;
        t = ( e >= 0 ? 0 : ( -e >= c ? 1 : -e/c ) );
    };
    auto f4 = [&s,&t,&det,&a,&b,&c,&d,&e,&f](){
        if (d < 0){
            t = 0;
            s = ( -d >= a ? 1 : -d/a );

        }else{
            s = 0;
            t = ( e >= 0 ? 0 : ( -e >= c ? 1 : -e/c ) );
        }


    };
    auto f5 = [&s,&t,&det,&a,&d](){
        t = 0;
        s = ( d >= 0 ? 0 : ( -d >= a ? 1 : -d/a ) );
    };
    auto f6 = [&s,&t,&det,&a,&b,&c,&d,&e,&f](){
        T tmp0 = b+e;
        T tmp1 = a+d;
        if ( tmp1 > tmp0 ){
            T  numer = tmp1 - tmp0;
            T denom = a-2*b+c;
            t = ( numer >= denom ? 1 : numer/denom );
            s = 1-t;
        }else {
            t = 0;
            s = ( tmp1 <= 0 ? 1 : ( d >= 0 ? 0 : -d/a ) );
        }
    };




    if ( s+t <= det )
    {
        if ( s < 0 ) { if ( t < 0 ) { f4(); } else { f3(); } }
        else if ( t < 0 ) { f5(); }
        else { f0(); }
    }
    else
    {
        if ( s < 0 ) { f2();}
        else if ( t < 0 ) { f6(); }
        else { f1();}
    }

    product(s,E0,E0);product(t,E1,E1);
    add(B,E0,cp);add(cp,E1,cp);
    minusVec(P,cp,vec);
    return dot(vec,vec);
}
template <class T>
inline T point2TriDistance(const T *P, const T *v1, const T *v2, const T *v3,T *cp){
    return sqrt(point2TriSquareDistance(P,v1,v2,v3,cp));
}
template <class T>
inline T threeDet(const T *p1,const T *p2,const T *p3){
   return p1[0]*(p2[1]*p3[2]-p2[2]*p3[1]) - p1[1]*(p2[0]*p3[2]-p2[2]*p3[0]) + p1[2]*(p2[0]*p3[1]-p2[1]*p3[0]);
}
template <class T>
inline bool threeplaneIntersection(const T *p1,const T *p2,const T *p3,T *point){

    T det = threeDet(p1,p2,p3);
    if(fabs(det)<1e-5)return false;
    T css[3];
    for(int i=0;i<3;++i)point[i] =0;
    cross(p2,p3,css);
    weightedAddVec(-p1[3],css,point);
    cross(p3,p1,css);
    weightedAddVec(-p2[3],css,point);
    cross(p1,p2,css);
    weightedAddVec(-p3[3],css,point);

    product(1/det,point,point);

    return true;


}
template <class T>
inline T pointNor2para4(const T *point, const T *nor, T *para){
    for(int i=0;i<3;++i)para[i] = nor[i];
    normalize(para);
    para[3] = -dot(point,para);

	return para[3];
}

template <class T>
inline bool planeSegIntersectionTest(const T *e1p1, const T *e1p2,const T *e2p1, const T *e2p2 ){

    T vec1[3],vec2[3],vec3[3],d1[3],d2[3];

    minusVec(e1p2,e1p1,vec1);
    minusVec(e2p1,e1p1,vec2);
    minusVec(e2p2,e1p1,vec3);
    cross(vec1,vec2,d1);
    cross(vec1,vec3,d2);
    T a = dot(d1,d2);

    //T a = threeDet(vec1,vec2,vec3);

    minusVec(e2p2,e2p1,vec1);
    minusVec(e1p1,e2p1,vec2);
    minusVec(e1p2,e2p1,vec3);

    cross(vec1,vec2,d1);
    cross(vec1,vec3,d2);
    T b = dot(d1,d2);

    if(a<0 && b<0)return true;
    else return false;
}

template <class T>
bool isPointOnSeg(const T *v, const T *v1, const T *v2, const T THRES){

    double vec1[3],vec2[3];

    minusVec(v,v1,vec1);
    minusVec(v2,v,vec2);

    if( normVec(vec1)<THRES ||  normVec(vec2)<THRES) return true;

    return fabs(cosine(vec1,vec2)-1)<THRES;
}


}//namespace


using namespace std;
inline void GetFolders(string path, vector<string>&files){
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (path.data())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            if(ent->d_type == DT_DIR)if(ent->d_name[0]!='.'){
                files.emplace_back(ent->d_name);
            }
        }
        closedir (dir);
    } else {
        /* could not open directory */
        perror ("could not open directory");
    }

}

inline void GetFiles(string path, vector<string>&files){
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (path.data())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            if(ent->d_type == DT_REG)files.emplace_back(ent->d_name);
        }
        closedir (dir);
    } else {
        /* could not open directory */
        perror ("could not open directory");
    }
}

inline void GetFiles(string path, vector<string>&files, string match_ext){
    GetFiles(path,files);
    vector<string>rev_files;
    for(auto &a: files){
        string ext;
        MyUtility::GetExtName(a,ext);
        if(ext==match_ext)rev_files.push_back(a);
    }
    swap(files,rev_files);
}


#endif




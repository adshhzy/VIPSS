/* Copyright 1991, Brown Computer Graphics Group.  
   All Rights Reserved. */


/* -----------------------------------------------------------------
 *
 * Several objects in Euclidean space
 *
 * -------------------------------------------------------------- */

#ifndef _RN_DEFS_DEFS
#define _RN_DEFS_DEFS

#include <WINSystemDefines.H>

/* -------------------  Constants  ------------------------------- */

/** \addtogroup HelperGeometry */
//@{

/// Default precision for double precision operations (default is $1e^{-15}$)
extern double RNEpsilon_d;
/// Default precision for floating point operations (default is $1e^{-6}$)
extern float RNEpsilon_f;

#ifndef M_PI
/// Pi
const double M_PI =3.1415926535897932384626433832795028842;
#endif


/* -------------------  Inline functions  ------------------------------- */

/** \brief Is zero within given epsilon, double precision.
 @param in_d number to compare
 @param in_dEps anything smaller is zero
 @returns TRUE/FALSE
 */
inline WINbool RNIsZero( const double in_d, const double in_dEps = RNEpsilon_d ) 
{ return (in_d < in_dEps && in_d > -in_dEps)? TRUE : FALSE; }

/** \brief Is zero within given epsilon, floating point.
 @param in_f number to compare
 @param in_fEps anything smaller is zero
 @returns TRUE/FALSE
 */
inline WINbool RNIsZero(const float in_f, const float in_fEps = RNEpsilon_f ) 
{ return (in_f < in_fEps && in_f > -in_fEps)? TRUE : FALSE; }


/** \brief Is the point/vector not at infinity?
    @param in_a is one of point, vector, covector
*/
template <class Affine>
inline
WINbool RNIsValid(const Affine &in_a) 
{ 
    for (int i = 0; i < in_a.Dim(); i++)
        if (fabs(in_a[i]) < 1e30) return TRUE;
    return FALSE;
}

/** \brief Zero out the point/vector
    @param out_a is one of point, vector, covector
*/
template <class Affine>
inline
Affine RNZero(Affine &out_a) 
{ 
    for (int i = 0; i < out_a.Dim(); i++)
        out_a[i] = 0;
    return out_a;
}

/** \brief Are the two numbers within epsilon of each other? Double precision
  @param in_d1 First number
  @param in_d2 Second number
  @param in_dEpsilon default is RNEpsilon_d
  */
inline WINbool RNApproxEqual( const double in_d1, const double in_d2,
                              const double in_dEpsilon = RNEpsilon_d)
{
    if (fabs(in_d1 - in_d2) > in_dEpsilon) return FALSE;
    return TRUE;
}

/** \brief Are the two numbers within epsilon of each other? Floating point
  @param in_f1 First number
  @param in_f2 Second number
  @param in_fEpsilon default is RNEpsilon_f
  */
inline WINbool RNApproxEqual( const float in_f1, const float in_f2,
                              const float in_fEpsilon = RNEpsilon_f)
{
    if (fabs(in_f1 - in_f2) > in_fEpsilon) return FALSE;
    return TRUE;
}

/** \brief Make sure dot product values are within angle ranges */
inline double DotToAngle( const double in_dDot ) 
{ return acos( WINminmax( in_dDot, -1.0 + RNEpsilon_d, 1.0 - RNEpsilon_d ) ); }


/** \brief Read binary  
 Read a single number in binary from ifstream. Note: won't work across operating systems that have different byte order (little indian versus big) 
*/
inline void ReadBinary( double &in_d, ifstream &in ) 
{ in.read( (char *) &in_d, sizeof(double) ); }

/** \brief Write binary  
Write a single number in binary to ofstream. Note: won't work across operating systems that have different byte order (little indian versus big) 
*/
inline void WriteBinary( const double in_d, ofstream &out ) 
{ out.write( (const char *) &in_d, sizeof(double) ); }


//@}

#include <utils/Rn_Affine1.H>
#include <utils/Rn_Affine2.H>
#include <utils/Rn_Affine3.H>
#include <utils/Rn_Affine4.H>
#include <utils/Rn_Projective.H>

/// Zero out the 2D homgoneous point/vector
template <class Coord>
inline
R2HPointTC<Coord> RNZero(R2HPointTC<Coord> &io_a) 
{ 
    io_a.u = 0;
    io_a.v = 0;
    io_a.w = 1;
    return io_a;
}

/// Zero out the 3D homgoneous point/vector
template <class Coord>
inline
R3HPointTC<Coord> RNZero(R3HPointTC<Coord> &io_a) 
{ 
    io_a.x = 0;
    io_a.y = 0;
    io_a.z = 0;
    io_a.w = 1;
    return io_a;
}



/* -------------------  External functions  ------------------------------- */

/**\addtogroup LinearGeometry */
//@{
/** \brief Closest point to a plane.

   Plane is given by normal and distance along normal. 

<br> Files: include/utils/Rn_defs.H src/utils/geometry/RNDefs.cpp <br>
   @param in_vecPlaneNormal The normal of the plane.
   @param in_dDistAlong Where along the normal the plane starts.
   @param in_pt Point to project
   @returns The closest point
*/
extern R3Pt 
RNClosestPointToPlane( const R3Vec &in_vecPlaneNormal,
                       double in_dDistAlong,
                       const R3Pt &in_pt );


/** \brief Intersect vector with plane. 

   Plane is given by normal and distance along normal. 

<br> Files: include/utils/Rn_defs.H src/utils/geometry/RNDefs.cpp <br>
   @param in_vecPlaneNormal The normal of the plane.
   @param in_dDistAlong Where along the normal the plane starts.
   @param in_ptRay Start point of ray to intersect with plane
   @param in_vecRay Direction of ray to intersect with plane
   @param out_ptIntersect The intersection point, if any
   @returns True if intersection
*/
extern WINbool 
RNIntersectRayWithPlane( const R3Vec &in_vecPlaneNormal,
                         double       in_dDistAlong,
                         const R3Pt  &in_ptRay,
                         const R3Vec &in_vecRay,
                         R3Pt        &out_ptIntersect );


/** \brief Calculate all curvature information, given all of the derivative information of a surface S(s,t) 

<br> Files: include/utils/Rn_defs.H src/utils/geometry/RNDefs.cpp <br>

  @param in_pt Point in space, S(s,t)
  @param in_vecDs First derivative of S(s,t) wrt s, d S(s,t) / ds
  @param in_vecDt First derivative of S(s,t) wrt t
  @param in_vecDss Second derivative of S(s,t) wrt s
  @param in_vecDtt Second derivative of S(s,t) wrt t
  @param in_vecDst Mixed derivative of S(s,t) wrt s and t
  @param out_dArea Infitesmal measure of area
  @param out_dCrv1 Minimum principle curvature
  @param out_dCrv2 Maximum principle curvature 
  @param out_vec1 Vector in (s,t) of minimum principle curvature
  @param out_vec2 Vector in (s,t) of maximum principle curvature
  @param out_vecMin Tangent vector of minimum principle curvature (out_vec1 mapped to the tangent plane)
  @param out_vecMax Tangent vector of maximum principle curvature (out_vec1 mapped to the tangent plane)
  @param out_vecNorm Unit surface normal at S(s,t)
*/
extern void
Curvature( const R3Pt &in_pt, const R3Vec &in_vecDs, const R3Vec &in_vecDt, 
           const R3Vec &in_vecDss, const R3Vec &in_vecDtt, 
           const R3Vec &in_vecDst, double &out_dArea, 
           double &out_dCrv1, double &out_dCrv2, 
           R2Vec &out_vec1, R2Vec &out_vec2, 
           R3Vec &out_vecMin, R3Vec &out_vecMax, R3Vec &out_vecNorm );

//@}
#endif


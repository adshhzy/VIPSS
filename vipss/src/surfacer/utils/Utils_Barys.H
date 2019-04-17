/* $Id: Utils_Barys.H,v 1.3 2007/10/03 18:53:11 cmg Exp $ */

/* Copyright 1991, Brown Computer Graphics Group.  All Rights Reserved. */



/* -------------------------------------------------------------------------
 *
 * Bary centric coords on a polygon
 *
 * ------------------------------------------------------------------------- */

#ifndef UTILS_BARY_DEFS
#define UTILS_BARY_DEFS

#include <utils/Rn_Polygon.H>

/* -----------------------  Constants  ------------------------------- */

/* -----------------------  Classes  ------------------------------- */

/** \class UTILSBarys Utils_Barys.H utils/Utils_Barys.H
\ingroup UtilitiesGeom
\brief Barycentric coords of a triangle, sphere, polygon

  Most of these methods are static functions. Can also create an instance
  by passing the class a polygon, then call operators.

b[0], b[1], b[2] such that sum b[i] = 1 and b[i] * p[i] = point.

<br>   Files: 
- include/utils/Utils_Barys.H
- src/utils/utils/UTILSBarys.cpp

*/
class UTILSBarys {
private:
    UTILSBarys( UTILSBarys & );
    UTILSBarys &operator=( UTILSBarys & );
    
    R2Polygon m_poly; 
    
public:
    /// Take barycentric coordinates and return a point in the polygon
    R2Pt operator()(const R3Pt &in_ptBary) const;
    
    /// Take a point in the polygon and return barycentric coordinates
    R3Pt Inverse(const R2Pt &in_ptPoly) const;
    
    /// BaryCoords, triangle
    static R3Pt BaryCoords( const R2Polygon &in_poly, const R2Pt &in_pt );    
    /// Bary coords, polygon
    static Array<double> BaryCoordsFull( const R2Polygon &in_poly, const R2Pt &in_pt );
    /// Reconstruct point, triangle
    static R2Pt ReconstructBaryCoords( const R2Polygon &in_poly, const R3Pt &in_ptBary );    ///
    /// Reconstruct point, polygon
    static R2Pt ReconstructBaryCoords( const R2Polygon &in_poly, const Array<double> &in_adBary );
    /// Spherical bary coords (triangle)
    static R3Pt SphericalBaryCoords( const Array<R3Vec> &in_avecTri, const R3Vec & );
    /// Reconstruct point on sphere (triangle)
    static R3Vec ReconstructSphericalBaryCoords( const Array<R3Vec> &in_avecTri, const R3Pt &in_ptBary );

    ///
    UTILSBarys &operator=( const R2Polygon &in_poly ) 
        { ASSERT( in_poly.Num_pts() == 3 ); m_poly = in_poly; return *this; }

    /// A three sided polygon
    UTILSBarys( const R2Polygon &in_poly ) : m_poly(in_poly) 
        { ASSERT( m_poly.Num_pts() == 3 ); }
    ///
    virtual ~UTILSBarys() {}
    
    ///
    void Print() const { m_poly.Print(); }
    ///
    static WINbool Test();
};

inline
R2Pt UTILSBarys::operator()(const R3Pt &in_ptB) const
{
   return R2Pt( m_poly[0][0] * in_ptB[0] + m_poly[1][0] * in_ptB[1] + m_poly[2][0] * in_ptB[2],
                m_poly[0][1] * in_ptB[0] + m_poly[1][1] * in_ptB[1] + m_poly[2][1] * in_ptB[2] );
                
}


#endif


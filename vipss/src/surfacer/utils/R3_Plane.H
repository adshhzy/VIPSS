#ifndef _R3_PLANE_DEFS_H
#define _R3_PLANE_DEFS_H


#include <utils/Rn_Defs.H>

class R3Line;
class R3Line_seg;

/** \class R3Plane R3_Plane.H utils/R3_Plane.H
  \ingroup LinearGeometry
  \brief A plane

 Defined by point, normal or distance, normal or Ax + By + Cz + D = 0. <br>
 The vector formed by dA, dB, dC is of unit length

 Files: <br>
  - include/utils/R3_Plane.H <br>
  - include/utils/R3_Plane_i.H <br>
  - src/utils/geometry/R3Plane.cpp
  */
class R3Plane {
protected:
    /// Plane stored as implicit equation, with the vector formed by A,B,C unit length
    double m_dA, m_dB, m_dC, m_dD;

public:
    /**@name Access to data */
    //@{
    ///
    R3Vec Normal() const { return R3Vec( m_dA, m_dB, m_dC ); }

    /// Return any point on the plane
    R3Pt PtOnPlane() const;

    /// How far along the normal is this plane?
    double DistOnNormal() const ;
    /** \brief Returns the tangent plane
       Note that these vectors are not unique.
      @param out_vecDs Returned unit length vector, perpendicular to normal and out_vecDt
      @param out_vecDt Returned unit length vector, perpendicular to normal and out_vecDs
      */
    void TangentPlane( R3Vec &out_vecDs, R3Vec &out_vecDt ) const;
    //@}

    /**@name Intersection */
    //@{
    /// Evaluates implicit equation and checks for (approximate) zero
    WINbool IsPtOnPlane( const R3Pt &in_pt ) const { return RNIsZero( m_dA * in_pt[0] + m_dB * in_pt[1] + m_dC * in_pt[2] + m_dD ); }
    /// Project the point along the normal to the plane
    R3Pt ProjectOnPlane( const R3Pt & ) const;
    /// How far is the point from the plane?
    double DistToPlane( const R3Pt & ) const;
    /// Evaluate implicit equation
    double Evaluate( const R3Pt &in_pt ) const { return m_dA * in_pt[0] + m_dB * in_pt[1] + m_dC * in_pt[2] + m_dD; }
    /** \brief intersect the ray with the plane
     @param in_ray The input ray
     @param out_dT Distance along the ray to the point of intersection (if any)
     @param out_pt Point of intersection (if any)
     @returns TRUE if there was an intersection
    */
    WINbool IntersectRay( const R3Ray &in_ray, double &out_dT, R3Pt &out_pt ) const;
    /// Just return point
    R3Pt IntersectRay( const R3Ray &in_ray ) const;
    /** \brief Intersect the segment with the plane
      If there is an intersection with the line containing the line segment and the plane,
      then returns that intersection, even if it isn't in contained in the line segment. Function
      will still return false in this case.
      @param in_seg The 3D line segment
      @param out_dT Distance along the segment of the intersection (if any)
      @param out_pt Point of intersection (if any)
      @returns TRUE if there was an intersection */
    WINbool Intersect( const R3Line_seg &in_seg, double &out_dT, R3Pt &out_pt ) const;

    /** \brief Intersect the line with the plane
      @param in_line The 3D line 
      @param out_dT Distance along the segment of the intersection (if any)
      @param out_pt Point of intersection (if any)
     @returns TRUE if there was an intersection */
    WINbool Intersect( const R3Line &in_line, double &out_dT, R3Pt &out_pt ) const;
    WINbool IntersectPlane( const R3Plane &in_plane, R3Line &out_line ) const;
    //@}

    /**@name Transformation */
    //@{
    /// Rotate and translate the plane, returning a new one
    R3Plane Transform( const R4Matrix & ) const;
    /// Matrix that takes the plane so that the normal points in the y direction
    R3Matrix ToXYZ() const;
    //@}

    /**@name Constructors, comparitors */
    //@{
    /// Scale the length of the normal to be unit length
    void Normalize();
    /// Flip the normal;
    void Flip();
    /// Approximate comparison (RNEpsilon_d)
    WINbool operator==( const R3Plane & ) const;
    /// Copy plane
    R3Plane & operator=( const R3Plane & );
    ///
    R3Plane( const R3Plane & );
    /// Construction from a point and a normal
    R3Plane( const R3Pt &in_ptOnPlane, const R3Vec &in_vecNormal );
    /// Three non-colinear points form a plane
    R3Plane( const R3Pt &in_pt1, const R3Pt &in_pt2, const R3Pt &in_pt3 );
    /// Distance along normal and normal
    R3Plane( const double in_dDist, const R3Vec &in_vecNormal );
    /// Implicit equation
    R3Plane( const R4Pt &in_ptABCD );
    /// Implicit equation
    R3Plane( const double in_dA, const double in_dB, const double in_dC, const double in_dD ) : m_dA(in_dA), m_dB(in_dB), m_dC(in_dC), m_dD(in_dD) { }
    /// Default plane is x,z plane
    R3Plane(  ) : m_dA(0.0), m_dB(1.0), m_dC(0.0), m_dD(0.0) { }
    ///
    ~R3Plane() {}
    //@}

    /**@name Constructors, comparitors */
    //@{
    /// Writes plane as in implicit equation
    void Write(ofstream &out);
    /// Reads plane as in implicit equation
    void Read(ifstream &in);
    //@}

    ///
    static WINbool Test();
};

#endif


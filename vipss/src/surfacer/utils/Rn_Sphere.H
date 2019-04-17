#ifndef _RN_DEFS_SPHERE_H
#define _RN_DEFS_SPHERE_H

#include "Rn_Defs.H"
#include "Rn_Line.H"
#include "Rn_Line_seg.H"
#include "Rn_Polygon.H"

/** \class R2Sphere Rn_Sphere.H utils/Rn_Sphere.H
 \ingroup ImplicitGeometry
 \brief A circle.

    Defined by a point and a radius.
<br> Files: 
 - include/utils/Rn_Sphere.H <br>
 - include/utils/R2_Sphere_i.H <br>
 - utils/geometry/RNSphere.cpp <br>
 - utils/geometry/R2Sphere.cpp <br>
*/
class R2Sphere {
protected:
    double m_dRadius;
    R2Pt   m_ptCenter;
    
public:
    /**@name Access to data */
    //@{
    /// Constant function returns 2
    inline int     Dim() const         { return 2; }
    /// pi r^2
    inline double  Area() const        { return M_PI * Radius() * Radius(); }
    /// Circumfrance 2 pi r 
    inline double  SurfaceArea() const { return 2.0 * M_PI * Radius(); }
    ///
    inline double  Radius()      const { return m_dRadius; }
    ///
    inline double &Radius()            { return m_dRadius; }
    ///
    inline const R2Pt &Center()  const { return m_ptCenter; }
    ///
    inline R2Pt       &Center()        { return m_ptCenter; }
    //@}

    /**@name Queries */
    //@{
    /// not on bdry
    WINbool Inside ( const R2Pt &in_pt ) const;
    /// On the circumfrance?
    WINbool On     ( const R2Pt &in_pt ) const;
    /// On Boundray or inside
    WINbool OnOrInside ( const R2Pt &in_pt ) const;
    ///
    WINbool operator==(const R2Sphere &p) const;
    /** \brief Return the positive sqrt value.
        $(c_0, c_1)$ is the center point, $r$ is the radius. <br>
        $x = frac{2 c_0 + Sqrt[4 c_0^2  - 4 (c_1^2  + c_0^2  - r^2  - 2 c_1 x + x^2 )}{2}$ */
    double  X_pos(double y) const;
    /** \brief Return the positive sqrt value.
        $(c_0, c_1)$ is the center point, $r$ is the radius. <br>
        $y = frac{2 c_1 + Sqrt[4 c_1^2  - 4 (c_0^2  + c_1^2  - r^2  - 2 c_0 x + x^2 )}{2}$ */
    double  Y_pos(double x) const;
    
    /** \brief Return the negative sqrt value.
        $(c_0, c_1)$ is the center point, $r$ is the radius. <br>
        $x = frac{2 c_0 - Sqrt[4 c_0^2  - 4 (c_1^2  + c_0^2  - r^2  - 2 c_1 x + x^2 )}{2}$ */
    double  X_neg(double y) const;
    /** \brief Return the negative sqrt value.
        $(c_0, c_1)$ is the center point, $r$ is the radius. <br>
        $y = frac{2 c_1 - Sqrt[4 c_1^2  - 4 (c_0^2  + c_1^2  - r^2  - 2 c_0 x + x^2 )}{2}$ */
    double  Y_neg(double x) const;

    /// Generate evenly spaced points, number determined by io_apt.num
    void Points( Array<R2Pt> &io_apt, const double in_dTS = 0.0, const double in_dTE = 2.0 * M_PI ) const;
    //@}

    /**@name Intersect routines */
    //@{
    /** Closest point on the circle
        @param in_p1 Input point
        @return The closest point */
    R2Pt Closest( const R2Pt  & in_p1 ) const;
    /** Intersect the circle with the line.
        @param in_l The line
        @param out_p1 Will be set to first interesection, if any.
        @param out_p2 Will be set to second interesection, if any.
        @return The number of intersections (0,1, or 2). */
    int Intersect( const R2Line  & in_l, 
                   R2Pt          & out_p1, 
                   R2Pt          & out_p2 ) const;
    /** Intersect the circle with the line segment.
     @param in_lseg The segment to intersect with
     @param out_p1 Will be set to first interesection, if any.
     @param out_p2 Will be set to second interesection, if any.
     @return The number of intersections (0,1, or 2). */
    int Intersect( const R2Line_seg & in_lseg, 
                  R2Pt             & out_p1, 
                  R2Pt             & out_p2) const;
    /** Intersect the circle with the circle
     @param in_circ The circle to intersect with
     @param out_p1 Will be set to first interesection, if any.
     @param out_p2 Will be set to second interesection, if any.
     @return The number of intersections (0,1, or 2). */
    int Intersect( const R2Sphere  & in_circ, 
                  R2Pt             & out_p1, 
                  R2Pt             & out_p2) const;
    /** Intersect the circle with the polygon
        @param in_poly The polygon
        @return The polygon with all points outside the circle removed, and any intersections with the circle added */
    R2Polygon Intersect( const R2Polygon & in_poly ) const;
    /** Overlaps. */
    WINbool Overlaps( const R2Polygon & in_poly ) const;
    /** Overlaps. */
    WINbool Overlaps( const R2Sphere & in_circ ) const;
    //@}
   

    /**@name Constructors, assignment */
    //@{
    ///
    R2Sphere &operator=( const R2Sphere &in_sphere );
    ///
    R2Sphere(const R2Pt &in_ptCenter, double in_dR);
    /// Circle from 3 points
    R2Sphere(const R2Pt &in_ptP1, const R2Pt &in_ptP2, const R2Pt &in_ptP3 );
    /// Default center is (0,0)
    R2Sphere(double in_dR = 1.0);
    //@}

    ///
    inline virtual ~R2Sphere() {}

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes center point and radius to file.
    void    Write(ofstream &out) const ;
    /// Reads center point and radius from file.
    WINbool Read(ifstream &in);
    /// Writes center point and radius to file.
    void    WriteBinary(ofstream &out) const ;
    /// Reads center point and radius from file.
    void ReadBinary(ifstream &in);
    //@}

    /// Test intersection routines
    static WINbool Test();
};

/** \class R3Sphere Rn_Sphere.H utils/Rn_Sphere.H
 \ingroup ImplicitGeometry
 \brief A sphere defined by a point and a radius.

<br> Files: 
 - include/utils/Rn_Sphere.H   <br>
 - include/utils/R2_Sphere_i.H <br>
 - utils/geometry/RNSphere.cpp <br>
 - utils/geometry/R3Sphere.cpp <br>
*/
class R3Sphere {
protected:
    double m_dRadius;
    R3Pt   m_ptCenter;
    
public:
    /**@name Access to data */
    //@{
    /** Dimension.
        This function is constant.
        @return 3
        */
    inline int     Dim() const  { return 3; }
    /// $ 4/3 pi r^3 $
    inline double  Area() const;
    /// $ 4 pi r^2 $
    inline double  SurfaceArea() const { return 4.0 * M_PI * Radius() * Radius(); }
    ///
    inline double  Radius() const      { return m_dRadius; }
    ///
    inline double &Radius()            { return m_dRadius; }
    ///
    inline const R3Pt &Center() const { return m_ptCenter; }
    ///
    inline R3Pt       &Center()       { return m_ptCenter; }
    //@}

    /**@name Queries */
    //@{
    ///
    WINbool Inside ( const R3Pt &in_pt ) const;
    /// On the surface?
    WINbool On     ( const R3Pt &in_pt ) const;
    ///
    WINbool operator==( const R3Sphere &in_sphere ) const;
    /** Intersect the sphere with the line segment.
     @param in_lseg The segment to intersect with
     @param out_p1 Will be set to first interesection, if any.
     @param out_p2 Will be set to second interesection, if any.
     @return The number of intersections (0,1, or 2). */
    int Intersect( const R3Line_seg & in_lseg, 
                   R3Pt             & out_p1, 
                   R3Pt             & out_p2) const;
    ///
    Array<R3Pt> IntersectRay( const R3Ray &in_ray, Array<double> &out_adT ) const;
    //@}

    /**@name Constructors, assignment */
    //@{
    ///
    R3Sphere &operator=(const R3Sphere &p);
    ///
    R3Sphere(const R3Pt &in_ptCenter, double in_dR);
    /// Default center is (0,0,0)
    R3Sphere(double in_dR);
    //@}

    ///
    inline virtual ~R3Sphere() {}

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes center point and radius to file.
    void    Write(ofstream &out) const ;
    /// Reads center point and radius from file.
    WINbool Read(ifstream &in);
    //@}

    /// Test intersection routines
    static WINbool Test();
};

#include "R2_Sphere_i.H"
#include "R3_Sphere_i.H"
#endif


#ifndef _RN_DEFS_ELLIPSE_H
#define _RN_DEFS_ELLIPSE_H

#include "Rn_Defs.H"
#include "Rn_Line.H"
#include "Rn_Line_seg.H"

/** \class R2Ellipse Rn_Ellipse.H utils/Rn_Ellipse.H
  \ingroup ImplicitGeometry
  \brief  An ellipse

   Keeps both forms (implicit and explicit) because they're both useful, and expensive to convert between

<p>    Defined by x^2/a^2 + y^2/b^2 -1 = 0 
<p>    Rotated by theta and translated by (x0,y0)

<br> Files: 
 - include/utils/Rn_Ellipse.H <br>
 - include/utils/R2_Ellipse_i.H <br>
 - src/utils/geometry/RNEllipse.cpp <br>
 - src/utils/geometry/R2Ellipse.cpp <br>
*/
class R2Polygon;

class R2Ellipse {
protected:
    /// Implicit equation
    double m_adCoef[6];
    /// x and y axis length, rotation
    double m_dA, m_dB, m_dTheta;
    /// Center of ellipse
    R2Pt   m_ptCenter;

    /// Set center, theta, radius from implicit equation
    void SetCenterTheta();
    /// Set implicit equation from center, theta, radii
    void SetImplicit();
public:
    /**@name Access to data */
    //@{
    /// Constant function returns 2
    inline int     Dim() const         { return 2; }
    /// 
    inline double  Area() const        { return M_PI * XRadius() * YRadius(); }
    ///
    inline double  XRadius()     const { return m_dA; }
    ///
    inline double  YRadius()     const { return m_dB; }
    ///
    inline const R2Pt &Center()  const { return m_ptCenter; }
    /// Clockwise rotation
    inline double  Rotation()    const { return m_dTheta; }
    /// Vector from mid point to x extreme
    R2Vec          XAxis() const;
    /// Vector from mid point to y extreme
    R2Vec          YAxis() const;
    /// Angle of rotation (clockwise)
    double         GetTheta() const { return m_dTheta; }
    /// Implicit equation coefficients
    double         Coefs( const int in_i ) const { return m_adCoef[in_i]; }
    /// Check that the implicit equation defines an ellipse or cirlce (not hyperbola/line)
    WINbool        IsEllipse() const;
    /// Eval point on ellipse, 0 - 2pi
    R2Pt           operator()( const double dT ) const ;
    //@}

    /**@name Queries */
    //@{
    /// To the place where the A axis lies along the x axis and centered at the origin
    R3Matrix ToNormal() const;
    /// From the place where the A axis lies along the x axis and centered at the origin
    R3Matrix FromNormal() const;
    /// Evaluate implicit equation
    double  operator()( const R2Pt & ) const;
    /// Is point inside?
    WINbool Inside ( const R2Pt &in_pt ) const;
    /// On the circumfrance?
    WINbool On     ( const R2Pt &in_pt ) const;

    /// Project onto ellipse, unless point inside
    R2Pt ClosestPt( const R2Pt &in_pt ) const;
    /// Get x as a function of y. Only valid if y is between y extremes of ellipse
    WINbool XFromY( const double in_dY, double &out_dX1, double &out_dX2 ) const;
    /// Get y as a function of x. Only valid if x is between x extremes of ellipse
    WINbool YFromX( const double in_dX, double &out_dY1, double &out_dY2 ) const;
    //@}

    /**@name Intersections */
    //@{
    /** \brief Do they overlap?
      This is an algebraic solution, but it is very numerically unstable so it is better to use the divide and conquer approach below
      @param in_ell Ellipse to intersect with 
      @returns If they overlap y/n */
    WINbool Overlap( const R2Ellipse &in_ell ) const;
    /// Do they overlap? Checks polygon edges for intersection with ellipse
    WINbool Overlap( const R2Polygon &in_poly ) const;
    /** \brief Do they overlap?
      This is a divide and conquer approach. Chops up the boundary into in_iDiv pieces, checks the pieces for inclusion and continues.
      Optional parameters indicate what section of the boundary to check for inclusion in the other ellipse.
      @param in_ell Ellipse to intersect with 
      @param in_iDiv Number of divisions to start with (8-32 are usually good start values)
      @param in_dTS Start point of boundary to check (0 to 2pi)
      @param in_dTE End point of boundary to check (0 to 2pi), must be greater than in_dTS
      @returns If they overlap y/n */
    WINbool Overlap( const R2Ellipse &in_ell, const int in_iDiv, const double in_dTS = 0.0, const double in_dTE = 2.0 * M_PI ) const;

    /** \brief Find the actual intersection points
      This function doesn't exist.
      @param in_ell Ellipse to intersect with 
      @param out_apt Intersection points */
    void Intersect( const R2Ellipse &in_ell, Array<R2Pt> &out_apt ) const;

    /** \brief Find the intersection points where the line crosses the ellipse (if any)
      Analytical solution
      @param in_line Line to intersect with 
      @param out_pt1 If there is one (or more) intersections, returns the first point here
      @param out_pt2 If there are two intersections, returns the second point here
      @returns Number of intersections (0,1,2) */
    int Intersect( const R2Line &in_line, R2Pt &out_pt1, R2Pt &out_pt2 ) const;
    /** \brief Find the intersection points where the line segment crosses the ellipse (if any)
      Analytical solution. Computes line solution and then checks that intersection points are inside line segment. So
      will return the line intersection points, sorted by whether or not they are on the line segment
      @param in_line Line segment to intersect with 
      @param out_pt1 If there is one (or more) intersections, returns the first point here
      @param out_pt2 If there are two intersections, returns the second point here
      @returns Number of intersections (0,1,2) */
    int Intersect( const R2Line_seg &in_line, R2Pt &out_pt1, R2Pt &out_pt2 ) const;
    //@}

    /// Approximate equal with RNEpsilon_d
    WINbool operator==(const R2Ellipse &p) const;
    /// Approximate comparison of two ellipses. Compares rotation angles, center points, and radii
    static WINbool ApproxEqual(const R2Ellipse &in_ell1, const R2Ellipse &in_ell2, const double in_dEps = RNEpsilon_d ) ;
    /** \brief Generate evenly spaced points, number determined by io_apt.num
      @param io_apt Points are returned here. Number of points in io_apt determines number of samples.
      @param in_dTS Starting point (0,2pi) of sample points, default 0.0
      @param in_dTE Ending point (0,2pi) of sample points, default 2pi. in_dTE must be greater than in_dTS */
    void Points( Array<R2Pt> &io_apt, const double in_dTS = 0.0, const double in_dTE = 2.0 * M_PI ) const;
    //@}   

    /**@name Constructors, assignment */
    //@{
    /// Explicit construction
    R2Ellipse(const R2Pt &in_ptCenter, const double in_dRMinor, const double in_dRMajor, const double in_dRotation = 0 );
    /// Default center is (0,0), unit radius (i.e., a sphere)
    R2Ellipse(const double in_dRMinor = 1, const double in_dRMajor = 1);
    /// Implicit equation: A x^2 + B x y + C y^2 + D x + E y + F = 0
    R2Ellipse(const double in_dA, const double in_dB, const double in_dC, const double in_dD, const double in_dE, const double in_dF);
    //@}

    ///
    inline virtual ~R2Ellipse() {}

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes center point, radius, rotation and implicit equation to file.
    void    Write(ofstream &out) const ;
    /// Reads center point, radius, rotation and implicit equation to file.
    WINbool Read(ifstream &in);
    ///  Writes center point, radius, rotation and implicit equation to file. Binary: Will not work cross platform (little indian vs. big indian)
    void    WriteBinary(ofstream &out) const ;
    ///  Reads center point, radius, rotation and implicit equation to file. Binary: Will not work cross platform (little indian vs. big indian)
    WINbool ReadBinary(ifstream &in);
    //@}

    /// Test intersection routines
    static WINbool Test();
};

/** \class R3Ellipse Rn_Ellipse.H utils/Rn_Ellipse.H
  \ingroup ImplicitGeometry
  \brief  A 3D ellipsoid shape. Not a true 3D ellipse.

    Defined by a center point, 3 axes, and two rotations (about x and y axis)

<br> Files: 
 - include/utils/Rn_Ellipse.H   <br>
 - include/utils/R3_Ellipse_i.H <br>
 - utils/geometry/RNEllipse.cpp <br>
 - utils/geometry/R3Ellipse.cpp <br>
*/
class R3Ellipse {
protected:
    /// Three axes and three Euler rotations
    double m_dA, m_dB, m_dC, m_dXRot, m_dYRot, m_dZRot;
    /// Center point
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
    ///
    inline double  XRadius() const    { return m_dA; }
    ///
    inline double &XRadius()          { return m_dA; }
    ///
    inline double  YRadius() const    { return m_dB; }
    ///
    inline double &YRadius()          { return m_dB; }
    ///
    inline double  ZRadius() const    { return m_dC; }
    ///
    inline double &ZRadius()          { return m_dC; }
    ///
    inline double  XRotation() const  { return m_dXRot; }
    ///
    inline double &XRotation()        { return m_dXRot; }
    ///
    inline double  YRotation() const  { return m_dYRot; }
    ///
    inline double &YRotation()        { return m_dYRot; }
    ///
    inline double  ZRotation() const  { return m_dZRot; }
    ///
    inline double &ZRotation()        { return m_dZRot; }
    ///
    inline const R3Pt &Center() const { return m_ptCenter; }
    ///
    inline R3Pt       &Center()       { return m_ptCenter; }
    /// Vector from center point to x extreme
    R3Vec          XAxis() const;
    /// Vector from center point to y extreme
    R3Vec          YAxis() const;
    /// Vector from center point to z extreme
    R3Vec          ZAxis() const;
    //@}

    /**@name Queries */
    //@{
    /// Rotation matrix taking three axes to x,y,z axes
    R3Matrix RotateToNormal() const;
    /// Rotation matrix taking x,y,z axes to three axes of ellipsoid
    R3Matrix RotateFromNormal() const;
    /// Rotation plus translation to take ellipsoid to origin, x,y,z axes
    R4Matrix ToNormal() const;
    /// Rotation plus translation to take origin to ellipsoid center, x,y,z to ellipsoid axes
    R4Matrix FromNormal() const;
    /// Extension of 2D implicit equation to 3D: returns negative number if outside, zero on boundary, positive inside
    double  operator()( const R3Pt & ) const;
    /// Is point inside the ellipsoid?
    WINbool Inside ( const R3Pt &in_pt ) const;
    /// On the surface?
    WINbool On     ( const R3Pt &in_pt ) const;
    /// Approximate comparison
    WINbool operator==( const R3Ellipse &in_Ellipse ) const;
    /// Generate evenly spaced points, number determined by io_apt.num
    void Points( Array<R3Pt> &io_apt, const double in_dTAlongXAxis ) const;
    /// Generate evenly spaced points, number determined by io_apt.num
    void PointsAndNormals( Array<R3Pt> &io_apt, Array<R3Vec> &io_avec, const double in_dTAlongXAxis ) const;
    //@}

    /**@name Constructors, assignment */
    //@{
    ///
    R3Ellipse &operator=(const R3Ellipse &p);
    ///
    R3Ellipse(const R3Pt &in_ptCenter);
    /// Default center is (0,0), radius 1
    R3Ellipse();
    //@}

    ///
    inline virtual ~R3Ellipse() {}

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

#include "R2_Ellipse_i.H"
#include "R3_Ellipse_i.H"
#endif


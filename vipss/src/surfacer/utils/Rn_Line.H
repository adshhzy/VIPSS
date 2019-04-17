#ifndef _RN_DEFS_LINE_H
#define _RN_DEFS_LINE_H

#include "Rn_Polynomial.H"

/** \class R2Line Rn_Line.H utils/Rn_Line.H
  \ingroup LinearGeometry
  \brief   A 2D line.

 Defined by a point and a unit vector. <br>
 Files: <br>
  - include/utils/Rn_Line.H <br>
  - include/utils/R2_Line_i.H <br>
  -  src/utils/geometry/R2Line.cpp <br>
  */
class R2Line {
protected:
    R2Pt   m_pt;
    R2Vec  m_vec;
    
public:
    /**@name Access to data */
    //@{
    /** Dimension.
        This function is constant.
        @return 2
        */
    inline int Dim()     const { return 2; }
    /// Anchor point
    inline const R2Pt  &Pt()   const { return m_pt; }
    /// Direction of line
    inline const R2Vec &Vec()  const { return m_vec; }
    /// Direction perpendicular to line (counter clockwise rotation of vec)
    inline R2Vec Normal()  const { return R2Vec( -m_vec[1], m_vec[0] ); }
    
    /// Slope of line
    inline double Slope() const;
    //@}

    /**@name Set data */
    //@{
    /// Set base point.
    inline void SetPt(const R2Pt &in_pt)    { m_pt = in_pt; }
    /// Set direction. Must be non-zero length
    inline void SetVec(const R2Vec &in_vec) ;
    //@}

    /**@name Project point onto line */
    //@{
    /** Project the point onto the line.
       @param in_pt Point to project
       @param out_ptClosest The closest point on the line
       @param out_t_on_line Parameter of closest point, i.e., operator(t) = closest point
       @param out_d_to_line Distance point is from line
       @return TRUE if point is on the line.
    */
    WINbool FindPtOnLine(const R2Pt & in_pt, 
                         R2Pt       & out_ptClosest, 
                         double     & out_t_on_line, 
                         double     & out_d_to_line ) const;
    /// Is the point within RNEpsilon_d of the line?
    WINbool IsPtOnLine    ( const R2Pt &in_pt ) const ;
    /// Distance from the given point to the line.
    double  Dist_to_line  ( const R2Pt &in_pt ) const ;
    /// Find the closest point on the line
    R2Pt    Project       ( const R2Pt &in_pt ) const ;
    /// Find t such that operator(t) = closest point 
    double  Projected_dist_on_line( const R2Pt &in_pt ) const ;
    //@}

    /**@name Types of line and line comparison queries */
    //@{
    /// Is the line vertical?
    inline WINbool Vertical()     const { return RNIsZero( m_vec[0] ); }
    /// Is the line horizontal?
    inline WINbool Horizontal()   const { return RNIsZero( m_vec[1] ); }
    /// Are the two lines parallel?
    friend WINbool IsParallel(const R2Line &in_line1, const R2Line &in_line2) ;
    /// Are the two lines perpendicular?
    friend WINbool IsPerpendicular(const R2Line &in_line1, const R2Line &in_line2) ;
    /// Are the two lines the same?
    WINbool operator==(const R2Line &in_line1) const ;
    //@}

    /**@name Evaluate line */
    //@{
    /** Y-Intercept.
        Intercept with the y axis, if exists
        @return zero if no intercept
    */
    inline double  Intercept()  const;

    /** Construct point along line.
        @return point + unit_vec * d
    */
    inline R2Pt operator()(double in_d) const { return m_pt + m_vec * in_d; }

    /** Evaluate line function.
        @return $x = (y - b ) / m$
    */
    double X( double in_dY ) const ;
    /** Evaluate line function.
        @return $y = m x + b$
    */
    double Y( double in_dX ) const ;
    /** Return line function as polynomial in $y$.
        @return $x = (y - b ) / m$
    */    
    RNPolynomial X() const ;
    /** Return line function as polynomial in $x$.
        @return $y = m x + b$
    */
    RNPolynomial Y() const ;
    ///
    void Implicit( double &out_adA, double &out_adB, double &out_adC ) const;
    //@}

    /**@name Intersect routines */
    //@{
    /** Intersect this line with another.
        @param in_l The line to intersect with.
        @param out_pt The intersection point.
        @return TRUE if the two lines intersect.
        */
    WINbool Intersect( const R2Line &in_l, R2Pt &out_pt ) const;
    //@}

    /**@name Construct lines */
    //@{
    /// Set this line equal to the input line
    R2Line &operator=(const R2Line &);

    /// Produce a line parallel to this one passing through the given point.
    R2Line Parallel(const R2Pt &in_pt) const ;
    /// Produce a line perpendicular to this one passing through the given point.
    R2Line Perpendicular(const R2Pt &in_pt) const ;

    /// Default constructor.
    R2Line() ;
    /// Construct and set to input line
    R2Line(const R2Line &in_line) ;
    /// Construct a line passing through the given point in the given direction
    R2Line(const R2Pt &in_pt, const R2Vec &in_vec) ;
    /// Construct a line passing through the two points
    R2Line(const R2Pt &in_pt1, const R2Pt  &in_pt2) ;

    /// Construct a line with the given slope and y-intercept.
    R2Line( double slope,  double y_intercept);
    /// Construct a line with the implicit form Ax + By + C
    R2Line( double in_dA,  double in_dB,  double in_dC);
    /// Construct a vertical line with the given x-intercept.
    R2Line( double x_intercept);
    //@}

    ///
    inline ~R2Line() {}

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes point and vector to file.
    void    Write ( ofstream &out ) const ;
    /// Reads point and vector from file.
    WINbool Read  ( ifstream &in );
    //@}

    /// Check closest point and intersection routines
    static WINbool Test();
};

/** \class R3Line Rn_Line.H utils/Rn_Line.H
  \ingroup LinearGeometry
  \brief   A 3D line.

  Defined by a point and a unit vector. <br>
 Files: <br>
  - include/utils/Rn_Line.H <br>
  - include/utils/R3_Line_i.H <br>
  -  utils/geometry/R3Line.cpp <br>
  */
class R3Line {
protected:
    R3Pt   m_pt;
    R3Vec  m_vec;
    
public:
    /**@name Access to data */
    //@{
    /** Dimension.
        This function is constant.
        @return 3
        */
    inline int Dim()     const { return 3; }
    /// Anchor point
    inline const R3Pt  &Pt()   const { return m_pt; }
    /// Direction of line
    inline const R3Vec &Vec()  const { return m_vec; }
    //@}

    /**@name Set data */
    //@{
    /// Set base point.
    inline void SetPt(const R3Pt &in_pt)    { m_pt = in_pt; }
    /// Set direction
    inline void SetVec(const R3Vec &in_vec);
    //@}

    /**@name Project point onto line */
    //@{
    /** Project the point onto the line.
       @param in_pt Point to project
       @param out_ptClosest The closest point
       @param out_t_on_line Parameter of closest point, i.e., operator(t) = closest point
       @param out_d_to_line Distance point is from line
       @return TRUE if point is on the line.
    */
    WINbool FindPtOnLine(const R3Pt & in_pt, 
                         R3Pt       & out_ptClosest, 
                         double     & out_t_on_line, 
                         double     & out_d_to_line ) const;
    /// Is the point within epsilon of the line?
    WINbool IsPtOnLine    ( const R3Pt &in_pt ) const ;
    /// Distance from the given point to the line.
    double  Dist_to_line  ( const R3Pt &in_pt ) const ;
    /// Find the closest point on the line
    R3Pt    Project       ( const R3Pt &in_pt ) const ;
    /// Find t such that operator(t) = closest point 
    double  Projected_dist_on_line( const R3Pt &in_pt ) const ;
    //@}

    /**@name Line comparison queries */
    //@{
    ///
    WINbool Intersect( const R3Line &, R3Pt &out_pt ) const;
    ///
    WINbool Intersect( const R3Line &, R3Pt &out_pt, double &out_dS, double &out_dT ) const;
    /// Are the two lines parallel?
    friend WINbool IsParallel(const R3Line &in_line1, const R3Line &in_line3) ;
    /// Are the two lines perpendicular?
    friend WINbool IsPerpendicular(const R3Line &in_line1, const R3Line &in_line3) ;
    /// Are the two lines the same, within RNEpsilon_d?
    WINbool operator==(const R3Line &in_line1) const ;
    //@}

    /**@name Evaluate line */
    //@{
    /** Construct point along line.
        @return point + unit_vec * d
    */
    inline R3Pt operator()(double in_d) const { return m_pt + m_vec * in_d; }
    //@}

    /**@name Construct lines */
    //@{
    /// Set this line equal to the input line
    inline R3Line &operator=(const R3Line &);
    /// Produce a line parallel to this one passing through the given point.
    R3Line Parallel(const R3Pt &in_pt) const ;
    /// Produce a line perpendicular to this one passing through the given point.
    R3Line Perpendicular(const R3Pt &in_pt) const ;

    /// Default constructor.
    inline R3Line() : m_pt(0,0,0), m_vec(0,1,0) {}
    /// Construct and set to input line
    inline R3Line( const R3Line &in_l ) { (*this) = in_l; }
    /// Construct from point and (unified) vector
    R3Line(const R3Pt &in_pt, const R3Vec &in_vec) ;
    /// Construct a line passing through the two points
    R3Line(const R3Pt &in_pt1, const R3Pt  &in_pt2) ;
    //@}

    ///
    inline ~R3Line() {}
    
    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes point and vector to file.
    void    Write(ofstream &out) const ;
    /// Reads point and vector from file.
    WINbool Read(ifstream &in);
    //@}

    /// Test closest point
    static WINbool Test();
};



#include "R2_Line_i.H"
#include "R3_Line_i.H"



#endif

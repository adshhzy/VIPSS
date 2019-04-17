#ifndef _RN_DEFS_LINE_SEG_H
#define _RN_DEFS_LINE_SEG_H

#include "Rn_Line.H"

/** \class R2Line_segTC Rn_Line_seg.H utils/Rn_Line_seg.H
  \ingroup LinearGeometry
  \brief   A 2D line segment.

    Defined by two points. This class is templated to support differentiation.

<br> Files: 
- include/utils/Rn_Line_seg.H <br>
- include/utils/R2_Line_seg_i.H <br>
- utils/geometry/R2Line_seg.cpp <br>
*/
template <class T>
class R2Line_segTC  {
protected:
    R2PointTC<T> m_pt1, m_pt2;
    
public:
    /**@name Access to data */
    //@{
    /// Constant function returns 2
    inline int Dim()          const { return 2; }
    ///
    inline const R2PointTC<T>   &P1() const { return m_pt1; }
    ///
    inline const R2PointTC<T>   &P2() const { return m_pt2; }
    ///
    inline       R2PointTC<T>   &P1()       { return m_pt1; }
    ///
    inline       R2PointTC<T>   &P2()       { return m_pt2; }
    //@}

    /**@name Queries: Comparison and type */
    //@{
    ///
    WINbool        Vertical() const;
    ///
    WINbool        Horizontal() const;
    /// Are the two line segments parallel?
    static WINbool IsParallelSeg(const R2Line_segTC<T> &, const R2Line_segTC<T> &) ;
  //    friend WINbool IsParallelSeg(const R2Line_segTC<T> &, const R2Line_segTC<T> &) ;
    /// Are the two line segments perpendicular?
    static WINbool IsPerpendicularSeg(const R2Line_segTC<T> &, const R2Line_segTC<T> &);
  //    friend WINbool PerpendicularSeg(const R2Line_segTC<T> &, const R2Line_segTC<T> &);
    ///
    WINbool        operator==(const R2Line_segTC<T> &) const;
    //@}

    /**@name Calculations on the segment */
    //@{
    /** Point on the segment.
        @return $s * p_1 + (1-s) * p_2$
    */
    R2PointTC<T> operator()( T s) const ;
    ///
    T            Length() const ;
    ///
    T            Slope()  const;
    //@}
    
    /**@name Project point onto line segment */
    //@{
    /** Project the point onto the line segement, clipping at the end points.
       @param in_pt Point to project
       @param out_ptClosest The closest point
       @param out_t_on_seg Parameter of closest point, i.e., operator(t) = closest point. Clips to [0,1].
       @param out_d_to_seg Distance point is from line
       @return TRUE if point is on the line.
    */
    WINbool     FindPtOnSeg( const R2PointTC<T> & in_pt, 
                             R2PointTC<T>       & out_ptClosest, 
                             T     & out_t_on_seg, 
                             T     & out_d_to_seg) const;
    /// Is the point within epsilon of the line segment? 
    WINbool      IsPtOnSeg   ( const R2PointTC<T> &in_pt ) const ;
    /// Distance from given point to closest point on line segment
    T            Dist_to_seg ( const R2PointTC<T> &in_pt ) const ;
    /// Closest point on segment
    R2PointTC<T> Project     ( const R2PointTC<T> &in_pt ) const;
    /// t between 0 and 1
    T            Projected_dist_on_seg( const R2PointTC<T> &in_pt ) const;
    //@}

    /**@name Intersection routines */
    //@{
    /** Intersect this line segment with another line segment.
        @param in_lineSeg The other line segment.
        @param out_pt The intersection point, if exists.
        @param out_s  The parameter value of the intersection point for this line segment.
        @param out_t  The parameter value of the intersection point for in_lineSeg
        @return TRUE if the line segments intersect.
        */
    WINbool Intersect( const R2Line_segTC<T> & in_lineSeg, 
                       R2PointTC<T>          & out_pt,
                       T                     & out_s, 
                       T                     & out_t) const;
    //@}
    
    /**@name Construct line segments */
    //@{
    /** Make a line segment perpendicular to this one.
        One end is on this line segment (or line) and the other is at the input point */
    R2Line_segTC<T>   Perpendicular ( const R2PointTC<T> &in_pt ) const ;
    /** Make a line segment parallel to this one.
        The new p1 is the input point, the new p2 is at the new p1 + the old (p2-p1) */
    R2Line_segTC<T>   Parallel      ( const R2PointTC<T> &in_pt ) const ;
    /// 
    R2Line_segTC<T>  &operator=     ( const R2Line_segTC<T> &in_line );
    /// Default constructor. Unit line segment.
    inline R2Line_segTC() : m_pt1(0,0), m_pt2(0,1) {}
    /// Construct from two points.
    inline R2Line_segTC( const R2PointTC<T> &p1, const R2PointTC<T> &p2 ) : m_pt1(p1), m_pt2(p2) {}
    /// Construct from another line.
    inline R2Line_segTC( const R2Line_segTC<T> &in_l ) { (*this) = in_l; }
    //@}

    ///
    inline ~R2Line_segTC() {}

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes 2 points to file.
    void    Write ( ofstream &out ) const ;
    /// Reads 2 points from file.
    WINbool Read  ( ifstream &in );
    /// Writes  2 points to file in binary. Does not work cross-platform (big vs little indian)
    void    WriteBinary ( ofstream &out ) const ;
    /// Reads  2 points to file in binary.  Does not work cross-platform (big vs little indian)
    WINbool ReadBinary  ( ifstream &in );
    //@}

    /// Test routine
    static WINbool Test();
};

/** \class R3Line_seg Rn_Line_seg.H utils/Rn_Line_seg.H
  \ingroup LinearGeometry
  \brief   A 3D line segment.

  Defined by two points.
<br> Files: 
 - include/utils/Rn_Line_seg.H <br>
 - include/utils/R3_Line_seg_i.H <br>
 - utils/geometry/R3Line_seg.cpp <br>
*/
class R3Line_seg  {
protected:
    R3Pt m_pt1, m_pt2;
    
public:
    /**@name Access to data */
    //@{
    /// Constant function returns 3
    inline int Dim()          const { return 3; }
    ///
    inline const R3Pt   &P1() const { return m_pt1; }
    ///
    inline const R3Pt   &P2() const { return m_pt2; }
    ///
    inline       R3Pt   &P1()       { return m_pt1; }
    ///
    inline       R3Pt   &P2()       { return m_pt2; }
    //@}

    /**@name Queries: Comparison and type */
    //@{
    ///
    WINbool        Vertical() const;
    ///
    WINbool        Horizontal() const;
    /// Are the two line segments parallel?
    friend WINbool IsParallelSeg(const R3Line_seg &, const R3Line_seg &) ;
    /// Are the two line segments perpendicular?
    friend WINbool PerpendicularSeg(const R3Line_seg &, const R3Line_seg &);
    ///
    WINbool        operator==(const R3Line_seg &) const;
    //@}

    /**@name Calculations on the segment */
    //@{
    /** Point on the segment.
        @return $s * p_1 + (1-s) * p_2$
    */
    R3Pt        operator()( double s) const ;
    ///
    double      Length() const ;
    //@}
    
    /**@name Project point onto line segment */
    //@{
    /** Project the point onto the line segement, clipping at the end points.
       @param in_pt Point to project
       @param out_ptClosest The closest point
       @param out_t_on_seg Parameter of closest point, i.e., operator(t) = closest point. Clips to [0,1].
       @param out_d_to_seg Distance point is from line
       @return TRUE if point is on the line.
    */
    WINbool     FindPtOnSeg( const R3Pt & in_pt, 
                             R3Pt       & out_ptClosest, 
                             double     & out_t_on_seg, 
                             double     & out_d_to_seg) const;
    /// Is the point within epsilon of the line segment? 
    WINbool     IsPtOnSeg   ( const R3Pt &in_pt ) const ;
    /// Distance from given point to closest point on line segment
    double      Dist_to_seg ( const R3Pt &in_pt ) const ;
    /// Closest point on segment
    R3Pt        Project     ( const R3Pt &in_pt ) const;
    /// t between 0 and 1
    double      Projected_dist_on_seg( const R3Pt &in_pt ) const;
    //@}

    /**@name Intersection routines */
    //@{
    /** Intersect this line segment with another line segment.
        @param in_lineSeg The other line segment.
        @param out_pt The intersection point, if exists.
        @param out_s  The parameter value of the intersection point for this line segment.
        @param out_t  The parameter value of the intersection point for in_lineSeg
        @return TRUE if the line segments intersect.
        */
    WINbool Intersect( const R3Line_seg & in_lineSeg, 
                       R3Pt             & out_pt,
                       double           & out_s, 
                       double           & out_t) const;
    /** Intersect with a plane.
        Plane is given by normal and distance along normal. 
        @param in_vecPlaneNormal The normal of the plane.
        @param in_dDistAlong Where along the normal the plane starts.
        @param out_dAlongSeg Where along the segment the point lies
        @param out_pt The intersection point, if it exists.
        @return TRUE if the segment intersects the plane.
    */
    WINbool Intersect( const R3Vec & in_vecPlaneNormal,
                       const double  in_dDistAlong,
                       double      & out_dAlongSeg,
                       R3Pt        & out_pt ) const;
    //@}
    
    /**@name Construct line segments */
    //@{
    /// Make a line segment perpendicular to this one and passing through the given point.
    R3Line_seg   Perpendicular ( const R3Pt &in_pt ) const ;
    /// Make a line segment parallel to this one and passing through the given point.
    R3Line_seg   Parallel      ( const R3Pt &in_pt ) const ;
    ///
    inline R3Line_seg  &operator=     ( const R3Line_seg &in_line );
    /// Default constructor. Unit line segment in y.
    inline R3Line_seg() : m_pt1(0,0,0), m_pt2(0,1,0) {}
    /// Construct from two points.
    inline R3Line_seg( const R3Pt &p1, const R3Pt &p2 ) : m_pt1(p1), m_pt2(p2) {}
    /// Construct from another line.
    inline R3Line_seg( const R3Line_seg &in_l ) { (*this) = in_l; }
    //@}

    ///
    inline ~R3Line_seg() {}

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes 2 points to file.
    void    Write ( ofstream &out ) const ;
    /// Reads 2 points from file.
    WINbool Read  ( ifstream &in );
    //@}

    /// Test routine
    static WINbool Test();
};


#include "R2_Line_seg_i.H"
#include "R3_Line_seg_i.H"

/** \brief 2D line segment typedef
   \ingroup LinearGeometry */
typedef R2Line_segTC<double> R2Line_seg;

/** \brief Intersect a plane with a line segment.
   \ingroup LinearGeometry

   Plane is given by normal and distance along normal. 

   Line segment is given by two points.

<br> Files: include/utils/Rn_Line_seg.H src/utils/geometry/RNDefs.cpp <br>

   @param in_vecPlaneNormal The normal of the plane.
   @param in_dDistAlong Where along the normal the plane starts.
   @param in_line_seg  line segment.
   @param out_pt The intersection point, if it exists.
   @return TRUE if the segment intersects the plane.
*/
extern WINbool 
RNIntersectPlaneSegment( const R3Vec &in_vecPlaneNormal,
                         double in_dDistAlong,
                         const R3Line_seg &in_line_seg,
                         R3Pt &out_pt);

#endif


#ifndef _RN_DEFS_CONTOUR_H
#define _RN_DEFS_CONTOUR_H

#include "Rn_Polygon.H"

/** \class R2Contour R2_Contour.H utils/R2_Contour.H
  \ingroup PolyGeometry
  \brief A collection of 2D polygons with a thickness and a location in space.
 Each plygon is defined by lists of points which make a polygon in a plane.
 This class isn't really finished and is a bit schitzoid.
 
   There are two uses for this class. The first is as a stack of contours. The second is as
    a collection of polygons in a plane. In both cases, the plane is defined by a plane normal and a distance (z level) 
    along that unit normal. For a stack of contours, the plane contains the first contour and the others
    are spaced dWidth apart from there.

 Files: <br>
  -  include/utils/R2_Contour.H <br>
  -  include/utils/R2_Contour_i.H <br>
  -   utils/geometry/R2Contour.cpp <br>
  */
class R2Contour {
protected:
    /// Planar polygon
    Array<R2Polygon> m_apolygons;        
    /// Vector at each vertex of each polygon
    Array< Array< R3Vec > > m_aavecNorm; 
    /// Spacing between slices
    double           m_dWidth;           
    /// Plane normal of all slices
    R3Vec            m_vecPlaneNormal;   
    /// Where plane starts along normal
    double           m_dZLevel;          
public:
    /**@name Access to data */
    //@{
    /** Dimension.
        This function is constant.
        @return 2
        */
    int     Dim() const { return 2; }
    /// Spacing between planes
    inline double Width() const  { return m_dWidth; }
    /// Where along normal planes start
    inline double ZLevel() const { return m_dZLevel; }
    /// Plane normal of all slices
    inline const R3Vec &PlaneNormal() const { return m_vecPlaneNormal; }
    /// Number of slices
    inline int Num_contours() const { return m_apolygons.num(); }
    /// Polygon for slice
    inline const R2Polygon & operator[] (int in_iPoly) const;
    /// Set normals of slice
    inline const Array<R3Vec> &Normals( const int in_iPoly ) const { return m_aavecNorm[in_iPoly]; }
    /// Do we have slices and are they valid?
    inline WINbool Valid() const;
    //@}

    /**@name \brief Setting data */
    //@{
    /// Set spacing
    inline double &    Width()   { return m_dWidth; }
    /// Set slicing plane direction
    inline void        SetPlane( const R3Vec &in_vec, double in_dZLevel );
    /// Set point in polygon in slice
    inline R2Polygon & operator[] (int in_iPoly) ;
    //@}

    /**@name Queries */
    //@{
    /// Total area of all slices
    double  Area() const ;
    /// Average center point of all slices
    R2Pt    Centroid() const;
    //@}

    /**@name Inside/outside/closest routines 
       Note: These routines assume you have a collection of polygons in the plane. Checks each polygon in turn.
    */
    //@{
    /// Is it on or inside the boundary? Checks all; returns true if true for one polygon.
    WINbool Inside( const R2Pt & ) const;
    /// Is it within epsilon of the boundary? Checks all; returns true if true for one polygon.
    WINbool IsBoundary( const R2Pt & ) const;
    ///  Checks all; returns true if true for all polygons.
    WINbool IsConvex() const;
    /** Use the winding test to see if the point is inside.
        Assumes the polygon is convex. Counts the number of crossings 
        on a ray from the point to infinity.  Checks all; returns true if true for one polygon. */
    WINbool Inside_winding(const R2Pt &) const;
    /// Project onto plane 
    R3Pt    Project( const R2Pt & ) const;
    /// Closest point on the boundary, if not inside already.  Checks all; returns true if true for one polygon. 
    R2Pt    Closest_pt( const R2Pt & ) const;
    //@}

    /**@name Comparisons with other contours */
    //@{
    /// Approx equal on the polygons themselves
    WINbool   operator==( const R2Contour & ) const;
    /// This could be very slow. Intersects all polygons of one contour with the polygons of the other.
    R2Contour Intersect(const R2Contour &in_contour) const;
    //@}
    
    
    /**@name Construction and assignment routines */
    //@{
    ///
    inline R2Contour &operator=(const R2Contour &);
    /// Add a polygon and normals
    inline void Add(const R2Polygon &, const Array<R3Vec> &in_avec);
    ///
    inline R2Contour( const R2Contour &in_poly );
    ///
    inline R2Contour(  );
    //@}

    ///
    virtual inline ~R2Contour() {}
    
    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes points to file.
    void    Write ( ofstream &out ) const ;
    /// Reads points from file.
    WINbool Read  ( ifstream &in );
    //@}

    /// Check closest point and intersection routines
    static WINbool Test();
};

#include "R2_Contour_i.H"

#endif


#ifndef _RN_DEFS_POLYGON_H
#define _RN_DEFS_POLYGON_H

#include "Rn_Line_seg.H"

/** \class R2Polygon Rn_Polygon.H utils/Rn_Polygon.H
 \ingroup PolyGeometry
 \brief A 2D polygon defined by an orderded list of 2D points. 

  Use this class if you're just making a polygon and discarding it, or if you have a concave polygon.
  If you're doing lots of inside/outside/clipping tests, the boxed polygon is faster because it
  includes a bounding box.

 Files: <br>
  - include/utils/Rn_Polygon.H <br>
  - include/utils/R2_Polygon_i.H <br>
  -  src/utils/geometry/R2Polygon.cpp <br>
  -  src/utils/geometry/RNPolygon.cpp (operations that are used for both 2D and 3D polygons)
  */
class R2Polygon {
protected:
    /// The point list
    Array<R2Pt> m_apts; 

    /// Used when intersecting polygons to turn intersection info into a polygon
    R2Polygon Make_polygon() const;

public:
    /** For storing scan line data. */
    class ScanLine {
    public:
        /// y line
        int m_iY;
        /// Start and stop vertices of two edges
        int m_iE1V1, m_iE1V2, m_iE2V1, m_iE2V2;
        /// Percentage along edges in y direction
        double m_dE1Perc, m_dE2Perc;
        /// Where along x
        int m_iX1, m_iX2;
    };

public:
    /**@name Access to data */
    //@{
    /** Dimension.
        This function is constant.
        @return 2
        */
    int     Dim() const { return 2; }
    /// Should be greater than 2
    inline int Num_pts() const { return m_apts.num(); }
    /// Access to array of points
    inline const R2Pt & operator[] (int in_iVert) const;
    /// For setting points
    inline       R2Pt & operator[] (int in_iVert) ;
    /// Accessor that wraps the indices around
    inline const R2Pt & PtWrap ( int in_iVert ) const ;
    /// Accessor that wraps the indices around
    inline R2Pt       & PtWrap ( int in_iVert ) ;
    /// All the points
    inline const Array<R2Pt> &Pts() const { return m_apts; }
    /// All the points
    inline Array<R2Pt> &SetPts() { return m_apts; }
    /// For n > 3, puts point in center and sums triangle areas
    double  Area() const ;
    /// Do we have at least three points? And are they valid?
    WINbool Valid() const;
    //@}

    /**@name Constructed data */
    //@{
    /// Returns unit up (0,0,1) or unit down (0,0,-1), depending upon orientation
    R3Vec       Normal() const;
    /// Average of all points
    R2Pt        Centroid() const;
    /// Line segment from vertex i to vertex i+1
    R2Line_seg  Edge( int inV1ToV2) const;
    /// Line from vertex i to vertex i+1
    R2Line      Line( int inV1ToV2) const;
    /// Returns total length of polygon
    double      SampleEvenly( Array<R2Pt> &out_aptSamples ) const;
    //@}

    /**@name Inside/outside/closest routines */
    //@{
    /// Is it on or inside the boundary?
    WINbool Inside( const R2Pt & ) const;
    /// for convex polygons only
    WINbool InsideQuick( const R2Pt & ) const;
    /// Is it within epsilon of the boundary?
    WINbool IsBoundary( const R2Pt & ) const;
    /// Checks outside angle
    WINbool IsConcave(const int in_iPt) const;
    /// Checks for outside angles less than pi/2
    WINbool IsConvex() const;
    /** \brief Use the winding test to see if the point is inside.

        Counts the number of crossings 
        on a ray from the point to infinity.  */
    WINbool Inside_winding(const R2Pt &) const;
    /// Project onto plane (does nothing)
    R2Pt    Project( const R2Pt & ) const;
    /// Closest point on the boundary, if point not inside already
    R2Pt    Closest_pt( const R2Pt & ) const;
    /// Project onto the boundary
    R2Pt    ProjectBoundary( const R2Pt & ) const;
    /// Return all of the intersections
    WINbool IntersectRay( const R2Pt &, const R2Vec &, Array<R2Pt> &out_apt ) const;
    /// Return all of the intersections plus where they happened
    WINbool IntersectRay( const R2Pt &, const R2Vec &, 
                          Array<R2Pt> &out_apt, Array<double> &out_adTRay,
                          Array<int> &out_aiSide, Array<double> &out_adSide ) const;
    /// Return all of the y-crossings, full data
    void    ScanConvert( Array<ScanLine> &out_aiScan ) const;
    /// Return just the y-crossings
    void    ScanConvert( Array<R2Pt_i> &out_aiptIndex ) const;
    //@}

    /**@name Comparisons with other polygons */
    //@{
    /// Uses approximate equal (RNEpsilon_d)
    WINbool   operator==( const R2Polygon & ) const;
    /// Not sure how well this works
    R2Polygon Intersect(const R2Polygon &in_polygon) const;
    /// Assumes both polygons are convex
    R2Polygon IntersectConvex(const R2Polygon &in_polygon) const;
    /// Assumes both polygons are convex, but will work if a vertex is inside the other
    WINbool Overlaps(const R2Polygon &in_polygon) const;
    /// Find all of the intersections (checks each edge in turn)
    Array<R2Pt> Intersect( const R2Line_seg &in_seg ) const;
    //@}
    
    
    /**@name Construction and assignment routines */
    //@{
    /// Copy points over
    inline R2Polygon &operator=(const R2Polygon &p);
    /// Copy points over
    inline R2Polygon( const R2Polygon &in_poly ) { (*this) = in_poly; }
    /// Default constructor is invalid polygon
    R2Polygon( int in_iNPts = 3);
    /// Initialize with the given point list
    inline R2Polygon( const Array<R2Pt> &in_apt ) { m_apts = in_apt; }

    /// Polygon with unit length sides centered at the origin
    static R2Polygon Unit_polygon( int num_sides ); 
    /// Remove all points
    void Clear() { m_apts.clearcompletely(); }
    /// Add a point in to the end of the list
    void operator +=( const R2Pt &in_pt ) { m_apts += in_pt; }
    /// Take out point in_i
    void Remove( const int in_i ) { m_apts.del( in_i, 1 ); }
    /// Insert a point at the given index
    void Insert( const R2Pt &in_pt, const int in_i ) { m_apts.addItem( in_pt, in_i ); }
    //@}

    ///
    virtual inline ~R2Polygon() {}
    
    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes points to file. Format is Number of points, points
    void    Write ( ofstream &out ) const ;
    /// Writes points to file in binary. Uses point binary write. Format is Number of points, points
    void    WriteBinary( ofstream &out ) const ;
    /// Reads points from file.  Format is Number of points, points
    virtual WINbool Read  ( ifstream &in );
    /// Reads points from file in binary. Uses point binary read. Format is Number of points, points
    virtual WINbool ReadBinary( ifstream &in );
    //@}

    /// Check closest point and intersection routines
    static WINbool Test();
};

/** \class R2BoxedPolygon Rn_Polygon.H utils/Rn_Polygon.H
 \ingroup PolyGeometry
 \brief A 2D polygon with a bounding box.

 Defined by a list of points. Adds max, min information to a 2D polygon. 
  NOTE: You must call Set_dim when the polgon is changed or created in order to update the bounding box. <br>
 Files: <br>
  -  include/utils/Rn_Polygon.H <br>
  -  include/utils/R2_Polygon_i.H <br>
  -  include/utils/R2_BoxedPolygon_i.H <br>
  -   src/utils/geometry/R2Polygon.cpp <br>
  -   src/utils/geometry/RNPolygon.cpp <br>
  */
class R2BoxedPolygon : public R2Polygon {
protected:
    /// Bounding box bounds
    double  m_dHMax, m_dHMin, m_dWMax, m_dWMin;
    /// Ordering of vertices in the plane. Determines normal direction.
    WINbool m_bClockwise; 
    /// Centroid of the bounding box
    R2Pt    m_ptCentroid; 

public:
    /**@name Access to data */
    //@{
    /// Bounding box height
    double Height() const { return m_dHMax - m_dHMin; }
    /// Bounding box width
    double Width()  const { return m_dWMax - m_dWMin; }
    /// Bounding box diagonal
    double Diagonal() const { return Length( R2Pt( m_dHMin, m_dWMin ) - R2Pt( m_dHMax, m_dWMax ) ); }
    /// Bounding box maximum height
    double HMax()   const { return m_dHMax; }
    /// Bounding box minimum height
    double HMin()   const { return m_dHMin; }
    /// Bounding box maximum width
    double WMax()   const { return m_dWMax; }
    /// Bounding box minimum width
    double WMin()   const { return m_dWMin; }
    /// Centroid of polygon
    const R2Pt & Centroid() const { return m_ptCentroid; }
    /// Centroid of bounding box, not polygon
          R2Pt   BBCentroid() const { return R2Pt( (m_dWMin + m_dWMax) * 0.5, (m_dHMin + m_dHMax) * 0.5 ); }
    /// Call when polygon has changed/been set.
    inline void Set_dim();
    //@}

    /**@name Inside/outside routines */
    //@{
    /// Not lying on boundary 
    WINbool IsOutsideBox( double dX, double dY ) const ;
    /// Not lying on boundary
    WINbool IsInsideBox( double dX, double dY ) const ;
    /// Assumes convexity
    WINbool Inside_quick( const R2Pt & ) const;
    /// Assumes both polygons are convex. Returned polygon may be empty size
    R2BoxedPolygon IntersectConvex(const R2BoxedPolygon &in_polygon) const;
    /// Slightly faster than checking for intersection
    WINbool Overlaps(const R2BoxedPolygon &in_polygon) const;
    /** \brief Intersect a ray; assumes convexity.
        Zeros out data if no intersection

       @param in_pt Starting point of ray
       @param in_vec Ray direction
       @param out_apt Points of intersection (if any)
       @param out_adTRay Where on the ray the intersection points happened (positive numbers)
       @param out_aiSide Which edges of the polygon were intersected
       @param out_adSide Where on the edge of the polygon the ray intersected (between [0,1])
       @returns If there were any intersection
       */
    WINbool IntersectRayQuick( const R2Pt &in_pt, const R2Vec &in_vec, 
                               Array<R2Pt> &out_apt, Array<double> &out_adTRay,
                               Array<int> &out_aiSide, Array<double> &out_adSide ) const;
    //@}

    /**@name Constructors and assignments */
    //@{
    /// Copies polygon and bounding box info
    inline R2BoxedPolygon &operator=( const R2BoxedPolygon & );
    /// Calls Set_dim for you after copying points
    inline R2BoxedPolygon &operator=( const R2Polygon & );
    /// Default, triangular polygon
    inline R2BoxedPolygon(int in_iNPts = 3);
    /// Calls Set_dim for you
    inline R2BoxedPolygon( const Array<R2Pt> &in_apt ) : R2Polygon( in_apt ) { Set_dim(); }
    /// Calls Set_dim for you
    inline R2BoxedPolygon( const R2Polygon & );
    /// Copies bounding box data as well
    inline R2BoxedPolygon( const R2BoxedPolygon & );
    //@}

    ///
    virtual inline ~R2BoxedPolygon() {}
    /// Reads points from file and calls Set_dim() (Data read is same as R2Polygon).
    virtual WINbool Read  ( ifstream &in ) { const WINbool bRes = R2Polygon::Read(in); Set_dim(); return bRes; }
};

/** \class R3Polygon Rn_Polygon.H utils/Rn_Polygon.H
 \ingroup PolyGeometry
 \brief Defined by a list of 3D points. 

 Files: <br>
  -  include/utils/Rn_Polygon.H <br>
  -  include/utils/R3_Polygon_i.H <br>
  -  src/utils/geometry/R3Polygon.cpp <br>
  -  src/utils/geometry/RNpolygon.cpp <br> Where operations common to both 2D and 3D polygons live
  */
class R3Polygon {
protected:
    /// The point list
    Array<R3Pt> m_apts; 

public:
    /**@name Access to data */
    //@{
    /** Dimension.
        This function is constant.
        @return 3
        */
    int     Dim() const { return 3; }
    /// Number of points
    inline int Num_pts() const { return m_apts.num(); }
    /// Access to points
    inline const R3Pt & operator[] (int in_iVert) const;
    /// For setting points
    inline       R3Pt & operator[] (int in_iVert) ;
    /// Accessor that wraps the indices around
    inline const R3Pt & PtWrap ( int in_iVert ) const ;
    /// Accessor that wraps the indices around
    inline R3Pt       & PtWrap ( int in_iVert ) ;
    /// All points
    const Array<R3Pt> &Pts() const { return m_apts; }
    
    /// Splits into triangles by adding a center point for Num_pts > 3
    double  Area() const ;
    /// Do we have at least three points? And are they valid?
    WINbool Valid() const;
    //@}

    /**@name Constructed data */
    //@{
    /// Note: not right if polygon is not planar
    R3Vec Normal() const;
    /// Average of points
    R3Pt        Centroid() const;
    /// Line segment from vertex i to vertex i+1
    R3Line_seg  Edge( int inV1ToV2) const;
    /// Line from vertex i to vertex i+1
    R3Line      Line( int inV1ToV2) const;
    //@}

    /**@name Inside/outside/closest routines */
    //@{
    /// Is it on or inside the boundary?
    WINbool Inside( const R3Pt & ) const;
    /// Is it within epsilon of the boundary?
    WINbool IsBoundary( const R3Pt & ) const;
    /// Checks outside angles
    WINbool IsConvex() const;
    /// Project onto plane of polygon
    R3Pt    Project( const R3Pt & ) const;
    /// Closest point on the boundary or inside
    R3Pt    Closest_pt( const R3Pt & ) const;
    /** \brief Project the point onto the triangle (assumes 3 points)
      * First projects the point onto the plane containing the triangle, then onto the triangle.
      @param in_pt Input point
      @param out_dDistSq Distance to the plane containing the triangle (squared)
      @param out_ptPlane The input point projected onto the plane containing the triangle
      @param out_ptTri The input point projected onto the triangle itself (on boundary if out_ptPlane not in triangle)
      @param out_adBary The three barycentric coordinates for the point inside the triangle. Will be in [0,1] if point in triangle, otherwise not
      @returns TRUE if point projects into triangle
      */
    WINbool ProjectTri( const R3Pt &in_pt, double &out_dDistSq, R3Pt &out_ptPlane, R3Pt &out_ptTri, Array<double> &out_adBary ) const;
    /** \brief Intersect the line segment with the triangle (assumes 3 points)
       Assumes polygon is a triangle and sees if the line segment passes through the triangle
       @param in_seg Line segment
       @param out_pt Point of intersection, if any
       @returns TRUE if there was an intersection */
    WINbool IntersectSegment( const R3Line_seg &in_seg, R3Pt &out_pt ) const;
    /// Actual distance (squared) from point to triangle
    double DistSq( const R3Pt &in_pt ) const;
    /// Worst possible distance (squared) from point to triangle
    double LBDistSq( const R3Pt &in_pt ) const;
    /** \brief Clip the polygon with the plane
       If all of the points are on the positive side of the plane, returns the polygon. If no points are on
       the positive side of the plane, returns FALSE and an empty polygon. Otherwise, throws out all points on the
       negative side of the plane and adds in any intersection points with the plane and the polygon edges.
      @param in_pt Point on the plane
      @param in_vec Normal of the plane
      @param out_polyClip The polygon clipped to the plane (i.e., throw out all points on opposite side of plane)
      @returns TRUE if there was an intersection */
    WINbool IntersectHyperplane(const R3Pt & in_pt, const R3Vec &in_vec, R3Polygon &out_polyClip ) const;
    //@}

    /**@name Comparisons with other polygons */
    //@{
    /// Uses approx equal on the points (RNEpsilon_d)
    WINbool   operator==( const R3Polygon & ) const;
    //@}

    /**@name Constructors and assignment */
    //@{
    /// Copy points
    inline R3Polygon &operator=( const R3Polygon & );
    /// Copy points
    inline R3Polygon( const R3Polygon &in_poly ) { (*this) = in_poly; }
    /// Copy points
    inline R3Polygon( const Array<R3Pt> &in_apt ) { m_apts = in_apt; }
    /// Default is size zero polygon
    R3Polygon( int in_iNPts = 3 );
    /// Remove all points (makes invalid)
    void Clear() { m_apts.clearcompletely(); }
    /// Add in a point at the end of the list
    void operator +=( const R3Pt &in_pt ) { m_apts += in_pt; }
    //@}

    ///
    inline  ~R3Polygon() {}
    
    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes points to file. (Number of points, points)
    void    Write ( ofstream &out ) const ;
    /// Reads points from file. (Number of points, points)
    WINbool Read  ( ifstream &in );
    /// Writes points to file. (Number of points, points)
    void    WriteBinary ( ofstream &out ) const ;
    /// Reads points from file. (Number of points, points)
    WINbool ReadBinary  ( ifstream &in );
    //@}

    /// Check closest point and intersection routines
    static WINbool Test();
};

#include "R2_Polygon_i.H"
#include "R2_Boxed_polygon_i.H"
#include "R3_Polygon_i.H"

#endif


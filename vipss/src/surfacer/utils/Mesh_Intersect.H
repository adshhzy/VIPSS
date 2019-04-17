#ifndef _MESH_INTERSECT_DEFS_H
#define _MESH_INTERSECT_DEFS_H


#include <utils/Mesh_PolygonSpatial.H>

/** \class GMeshIntersect Mesh_Intersect.H include/Mesh_Intersect.H
  \ingroup MeshesClosest
  \brief Intersecting a mesh/closest point 

  This class takes as input a GMesh, and builds a spatial partition. This makes
  finding the closest point/intersecting rays much faster. However, it takes
  a bit of time to set up, so if the mesh is constantly changing, building
  the structures is not worth the time.

  Note that if the mesh changes you must call ClearSpatialStructures(). The class
  will automatically build a new spatial partition the next time closest point/intersect
  is called.

  To get GMesh from PMeshLite, call Mesh().
*/
class GMeshIntersect {
protected:

    /// This smattering of variables is what's needed to
    /// compute intersections/closest points. 
    PolygonSpatial                   m_oPsp;
    MMap<const R3Polygon*,Face>      m_oMpf;
    MMap<Face, Face>                 m_oTriMap;
    MMap<Face, int>                  m_oSideMap;
    GMesh                            m_meshTrans;
    R4Matrix                         m_frameTrans;
    R4Matrix                         m_frameTransi;
    Array<R3Polygon>                 m_apoly;
    WINbool                          m_bInitd;

    const GMesh &m_mesh;


public:
    /**@name Accessor */
    //@{
    const GMesh &GetMesh() const { return m_mesh; }
    //@}
    
    /**@name Create and destroy data */
    //@{
    /// This will get called if structures don't exist
    void SetSpatialStructures();
    /// Call when mesh changes
    void ClearSpatialStructures();
    //@}

    /**@name Ray intersection */
    //@{
    /** Intersect the ray with the mesh. Mesh may be non-triangular.
        @param in_pt Starting point of ray
        @param in_vec Direction of ray
        @param out_faces For each intersection, the face intersected
        @param out_barys For each intersection, the barycentric coords of the intersection point (only 3 terms will be non-zero)
        @param out_pts The intersection points
        @return TRUE if intersections
        */
    WINbool IntersectRay( const R3Pt &in_pt, const R3Vec &in_vec, 
                          Array<Face >            & out_faces, 
                          Array< Array<double > > & out_barys,
                          Array< R3Pt > &out_pts ) ;

    /** Intersect the ray with the mesh. Mesh may be non-triangular.
        @param in_pt Starting point of ray
        @param in_vec Direction of ray
        @param out_faces For each intersection, the face intersected
        @param out_barys For each intersection, the barycentric coords of the intersection point (only 3 terms will be non-zero)
        @return TRUE if intersections
        */
    WINbool IntersectRay( const R3Pt &in_pt, const R3Vec &in_vec, 
                          Array<Face >            & out_faces, 
                          Array< Array<double > > & out_barys  ) ;

    /** Intersect the ray with the mesh.
        @param in_pt Starting point of ray
        @param in_vec Direction of ray
        @param out_pts The intersection points
        @return TRUE if intersections
        */
    WINbool IntersectRay( const R3Pt &in_pt, const R3Vec &in_vec, 
                          Array< R3Pt > &out_pts ) ;
    //@}

    /**@name Closest point */
    //@{
    /** Closest point on the mesh. 
        @param in_pt Point
        @param out_face The face containing the closest point
        @param out_adBarys The barycentric coords of the closest point
        @return The closest point
        */
    R3Pt    Closest( const R3Pt &in_pt, Face &out_face,  Array<double> &out_adBarys ) ;
    /** Closest vertex on the mesh. 
     @param in_pt Point
     @return The closest vertex
     */
    Vertex  ClosestVertex( const R3Pt &in_pt ) ;
    /** Closest point on the mesh.
        @param in_pt Point
        @param out_bIsBdry Was the closest point on the boundary of the mesh?
        @return The closest point
        */
    R3Pt    Closest( const R3Pt &in_pt, WINbool &out_bIsBdry ) ;
    /** Closest point on the mesh.
        @param in_pt Point
        @return The closest point
        */
    R3Pt    Closest( const R3Pt &in_pt ) ;

    /** Closest point on the mesh. 
        @param in_pt Point
        @param in_face The face to look in
        @param out_adBarys The barycentric coords of the closest point
        @return The closest point
        */
    R3Pt    ClosestInFace( const R3Pt &in_pt, const Face in_face,  Array<double> &out_adBarys ) const ;
    /** Closest point on the boundary mesh. Mesh better have a boundary.
        @param in_pt Point
        @param out_edge The edge the point was lying on
        @param out_blend The distance along the edge, [0,1]
        @return The closest point lying on the boundary
        */
    R3Pt    ClosestOnBoundary( const R3Pt &in_pt, Edge &out_edge, R2Pt &out_blend ) ;
    //@}

    ///
    GMeshIntersect( const GMesh &in_mesh);
    ///
    virtual ~GMeshIntersect();
};
#endif


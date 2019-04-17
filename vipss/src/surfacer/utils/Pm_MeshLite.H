  #ifndef _PMESHLITE_DEFS_H
#define _PMESHLITE_DEFS_H

#include <map>
#include <utils/Rn_Defs.H>
#include <utils/Utils_Color.H>
#include <utils/R3_Plane.H>
#include <boost/dynamic_bitset.hpp>

#include <utils/Mesh_GMesh.H>
#include <utils/Mesh_Intersect.H>

class R3Plane;

/** \class PMeshLite Pm_MeshLite.H utils/Pm_MeshLite.H
  \ingroup Meshes
  \brief A mesh, wrapper around GMesh

    See OGLMeshes_Texture for making texture maps  <br>

    All adjacency information is found in this class (e.g., an edge's
    two vertices), not on the vertex, edge, and face classes. <br>

    This class uses the GMesh underlying data structures <br>


    If using GMesh, call Compute_normals() if the
    mesh has been changed to update the normals.
    

<br><br> Files:
- include/utils/Pm_MeshLite.H
- utils/PMeshLite/PMeshLite.cpp
- utils/PMeshLite/PMeshLite_IO.cpp
- utils/PMeshLite/PMeshLite_Intersect.cpp
- utils/PMeshLite/PMeshLite_Tools.cpp */
class PMeshLite : private GMesh
{
protected:
    WINbool            m_bHasChanged; 

    /// Color 
    UTILSColor         m_colDefault;

    /// Check that the face can be added to the boundary
    WINbool CheckVs( Face in_of, const Array<Vertex > &in_aopVs ) const;

    /// grow disk subroutine
    WINbool EatWedge( Array<Vertex> &io_avBdry, Array<Edge> &io_aeBdry, boost::dynamic_bitset<> &io_abInDisk, const boost::dynamic_bitset<> &in_abNoCross ) const;

public:
    /**@name Tools */
    //@{
    /// Find a path from vertex 1 to any vertex in the list
    /// No crossing of marked vertices
    /// Use edge weights
    Array<Vertex> FindPath( const Vertex &in_v1, const Array<Vertex> &in_aopBdryVs, const boost::dynamic_bitset<> &in_abEdgeNoCross, const Array< double > &in_adEdgeWeights ) const;
    /// Find a path from vertex 1 to any vertex in the list
    Array<Vertex> FindPath( const Vertex &in_v1, const Array<Vertex> &in_aopBdryVs, bool in_bUseEdgeLen ) const;
    /// Don't cross
    Array<Vertex> FindPath( const Vertex &in_v1, const Array<Vertex> &in_aopBdryVs, const boost::dynamic_bitset<> &in_abEdgeNoCross, bool in_bUseEdgeLen ) const;
    /// Use edge weights
    Array<Vertex> FindPath( const Vertex &in_v1, const Array<Vertex> &in_aopBdryVs, const Array< double > &in_adEdgeWeights ) const;
    /// No edge weights, all cross
    Array<Vertex> FindPath( const Vertex &in_v1, const Vertex &in_v2, bool in_bUseEdgeLen = FALSE ) const;
    /// Use edge weights
    Array<Vertex> FindPath( const Vertex &in_v1, const Vertex &in_v2, const boost::dynamic_bitset<> &in_abEdgeNoCross, const Array< double > &in_adEdgeWeights ) const;
    /// Use no cross
    Array<Vertex> FindPath( const Vertex &in_v1, const Vertex &in_v2, const boost::dynamic_bitset<> &in_abEdgeNoCross, bool in_bUseEdgeLen ) const;
    /// Use edge weights
    Array<Vertex> FindPath( const Vertex &in_v1, const Vertex &in_v2, const Array<double> &in_adEdgeWeights ) const;
    /// run dijkstra's
    void GrowShortestPathTree( const Vertex &in_v1, const boost::dynamic_bitset<> &in_abNoCross, const Array<double> &in_adEdgeWeights, Array<Vertex> &out_avFrom, Array<double> &out_adDist, boost::dynamic_bitset<> &out_abTreeEdges, const double in_dStopGrow = 1e30) const;
    ///
    void DualMST( const Face &in_f1, const boost::dynamic_bitset<> &in_abNoCross, const Array<double> &in_adEdgeWeights, const double in_dStopGrow, Array<double> &out_ad, boost::dynamic_bitset<> & out_abTreeEdges) const;
    int GrowDisk(const Vertex &in_vSeed,
                 const boost::dynamic_bitset<> &in_abVsNoCross,
                 Array<Vertex> &out_aopVs,
                 Array<Face>   &out_aopFs,
                 int in_iMaxNumVs = -1 ) const;
    /// Keep boundary
    WINbool GrowDisk( Array<Vertex> &io_avBdry,
                      Array<Edge> &io_aeBdry, 
                      boost::dynamic_bitset<> &io_abInDisk, 
                      const boost::dynamic_bitset<> &in_abNoCross ) const;
    /// Add input face to disk and eat wedges
    int GrowDiskByFace(  const Face &in_face,
                                        boost::dynamic_bitset<> &io_abInDisk,
                                        Array<Vertex> & io_avBdry, 
                                        Array<Edge>   & io_aeBdry,
                                        Array<int>   & out_wedgeEaten) const;
    /// Grow a disk with input face as center until the disk has targeted number of faces
    bool GrowDiskByDijkstra(  const Face &in_face,
                          const int targetFsInDisk,
                          boost::dynamic_bitset<> &out_abFsInDisk,
                          Array<Vertex> & out_avBdry, 
                          Array<Edge>   & out_aeBdry ) const;

    WINbool AddFaceToDisk(  const Face &in_face,
                                  boost::dynamic_bitset<> &io_abFChart,
                                  Array<Vertex> & io_avBdry, 
                                  Array<Edge>   & io_aeBdry) const ;
    ///check if Vs and Es on the boundary of the disk-chart are in correct order
    void CheckBoundary( const Array<Vertex> &in_av, const Array<Edge> &in_ae) const;
    ///check if the chart will remain a disk after adding this face
    WINbool OkayToAdd( const Face &in_of, const Array<Vertex > &in_aopVs ) const;

    /// Find the shortest path between the two faces,
    /// Not crossing the marked edges
    WINbool FindPath( const Face in_oF1, const Face in_oF2,
                      const boost::dynamic_bitset<> &in_aiEdgeNoCross,
                      Array<Face >   &out_aopFs ) const;
    /// Get the connected components
    int GetComponents( Array< unsigned short > &out_ausLabels, Array<int> &out_aiLabels ) const;
    /// Find the biggest component in the mesh and return it
    void GetOneComponent( PMeshLite &out_mesh ) const;
    /// Return two vertices that are furthest apart in the set of selected vertices
    pair <int, int> ExtremePointsDistance( const boost::dynamic_bitset<> &in_abVs ) const ;
    /// Return two vertices that are furthest apart in x,y, or z dims in the set of selected vertices
    pair <int, int> ExtremePointsXYZ( const boost::dynamic_bitset<> &in_abVs ) const ;
    /// Exponential map. in_abNoCross is vertices not to use
    int ExponentialMap( const Vertex &in_v, const boost::dynamic_bitset<> &in_abNoCross, Array< std::pair<Vertex,R2Pt> > &out_avPtDisk, const double in_dStopGrow = 1e30) const;
    /// Exponential map, version II, in_abNoCross is vertices not to use
    int ExponentialMapII( const Vertex &in_v, const boost::dynamic_bitset<> &in_abNoCross, Array< std::pair<Vertex,R2Pt> > &out_avPtDisk, const double in_dEpsDone, const double in_dStopGrow = 1e30) const;
    /* Exponential map, version II, multiple source
     * @param in_av At least one vertex to serve as the source
     * @param out_avPtDisk For each vertex, the disk of vertices attached to that vertex and its 2d location
     * @param in_abNoCross List of vertices not to cross
     * @param in_dEpsDone what error tolerance to allow when placing vertices in disk
     * @param in_dStopGrow maximum size of disk
     */
    int ExponentialMapIIMultiSource( const Array<Vertex> &in_av,  
                                     Array< Array< std::pair< Vertex,R2Pt> > > &out_aavPtDisk, 
                                     const boost::dynamic_bitset<> &in_abNoCross,
                                     const double in_dEpsDone = 1e-12, const double in_dStopGrow = 1e30) const;
    /* Exponential map, version II, new source(s)
     * @param in_av At least one vertex to serve as the source
     * @param io_abValidDepth For each vertex, whether or not it's been reached
     * @param io_aiLabel For each vertex, which source it belongs to
     * @param io_adDist For each vertex, the distance to the source
     * @param out_avPtDisk For each vertex, the disk of vertices attached to that vertex and its 2d location
     * @param in_abNoCross List of vertices not to cross
     * @param in_dEpsDone what error tolerance to allow when placing vertices in disk
     * @param in_dStopGrow maximum size of disk
     */
    int ExponentialMapIIMultiSource( const Array<Vertex> &in_av,  
                                    boost::dynamic_bitset<> &io_abValidDepth,
                                    Array<int> &io_aiLabel,
                                    Array<double> &io_adDist,
                                    Array< Array< std::pair< Vertex,R2Pt> > > &out_aavPtDisk, 
                                    const boost::dynamic_bitset<> &in_abNoCross,
                                    const double in_dEpsDone = 1e-12, const double in_dStopGrow = 1e30) const;
    /// Get the connected components
    void MakeFlattenedMesh( const Array< std::pair<Vertex,R2Pt> > &in_aptVLoc, PMeshLite &out_mesh ) const;
    //@}


    /**@name Coloring */
    //@{
    ///
    void Color( const UTILSColor &in_col );
    ///
    inline void SetDefaultColor( const UTILSColor &in_col ) { m_colDefault = in_col; }
    ///
    inline const UTILSColor & DefaultColor(  ) const { return m_colDefault; }
    ///
    inline const UTILSColor & GetColor( Vertex in_v ) const { return in_v->GetColor(); }
    ///
    inline void SetColor( Vertex in_v, const UTILSColor &in_col ) { in_v->SetColor( in_col ); }
    ///
    void ColorByDistanceToInput( PMeshLite in_mesh );
    //@}

    /**@name Updating normals */
    //@{
    /// If point locations have changed
    void Compute_normals( const Array<Vertex> &aopVs, const Array<Face> &aopFs );
    /// Do all faces
    void Compute_face_normals(  );
    /// Do all
    void Compute_normals(  );
    //@}

    /**@name Geometry properties */
    //@{
    /// Centroid and scale of entire mesh
    R3Pt         GetCenter( R3Pt &out_ptScale ) const;
    /// Point furthest from centroid
    Vertex       GetExtreme() const;

    /// Point in a face
    R3Pt         GetPointInFace(       Face            in_oFace, 
                                 const Array<double> & in_adPerc ) const;
    R3Pt         GetPointInFace( Face        in_oFace,
                                 const R3Pt &in_ptBary ) const;
    /// Point in a face
    R3Pt         GetFaceCentroid( Face in_oFace ) const;
    /// Area of face
    double       GetFaceArea( const Face in_oFace ) const;
    ///
    void         GetFacePolygon( Face in_oFace, R3Polygon & out_poly ) const { GMesh::polygon(in_oFace, out_poly); }
    /// Only valid for closed, manifold, triangular meshes
    double       GetVolume() const;

    /// Point along an edge
    R3Pt         GetPointOnEdge( Edge in_oEdge, double in_dPerc ) const;
    ///
    double       GetEdgeLength( const Edge in_oEdge ) const;
    /// Returns array size of edge
    const Array<double> EdgeLengths() const;
    /// Center the mesh at the given point with the given raduis
    void Center( const R3Pt &in_ptCenter = R3Pt(0,0,0), const double in_dRadius = 0.5 );
	/// Make anything that seems like it should align, do so
	virtual void Regularize( double in_dPerc = 0.1 );
    ///
    Array<Face> PatchHoles( const bool in_bPatchBigWithVertex, const int in_iMaxHoleToPatch = -1 );
    /// Returns the total area
    double GetVertexAreas( Array<double> &out_adVArea ) const;

    ///migrated from Tim's code
    int AlphaBeta( const Face &in_face, const Vertex &in_vertex,
                             double &out_theta, 
                             double &out_alpha, double &out_beta,
                             double &t, double &a, double &b ) const;

    /// Oriented edge, returns angles
    int AlphaBeta( const Vertex &in_vertex, const Edge &in_edge,
                             double &out_theta1, double &out_theta2, 
                             double &out_alpha, double &out_beta,
                             double &out_delta, double &out_gamma ) const;

    /// Oriented edge, returns sum of cot of alpha beta
    double CotAlphaBeta( const Vertex &v1, const Vertex &v2 ) const;

    /// Oriented edge, returns sum of cot of gamma delta
    double CotGammaDelta( const Vertex &v1, const Vertex &v2 ) const;

    /// Oriented edge, returns sum of cot of alpha beta, gamma delta
    void CotAlphaBetaGammaDelta( const Vertex &v1, const Vertex &v2, double &out_dSumAlphaBeta, double &out_dSumGammaDelta ) const;
    
    ///
    int FaceAngles(const Face &in_face, const Vertex &in_vertex,
                             int &out_ThetaId, double &out_Theta, double &out_t, 
                             int &out_AlphaId, double &out_Alpha, double &out_a, 
                             int &out_BetaId, double &out_Beta, double &out_b ) const;
    
    void FaceAngles( const Face &in_face, Array<double> &out_adAngles ) const;

    /// For a face on the edge, return cot of all of the angles (theta is angle at vertex, alpha is ccw, beta is cw)
    int CotAngles( const Face &in_face, const Vertex &in_vertex,
                             int &out_ThetaId, double &out_dCotTheta, 
                             int &out_AlphaId, double &out_dCotAlpha, 
                             int &out_BetaId, double &out_dCotBeta ) const;

    //@}

    /**@name Changing the topology of the mesh */
    //@{
    ///
    inline Vertex AddVertex();
    ///
    inline Vertex AddVertex(const R3Pt &in_pt);
    ///
    inline Vertex AddVertex(const R3Pt &in_pt, const UTILSColor &);
    ///
    inline Vertex AddVertex(const R3Pt &in_pt, const R3Vec &in_vec, const UTILSColor &);
    /// Orient  all the faces the opposite way (turn it inside out)
    void           Flip();

    /// Make sure the incoming face is oriented correctly before adding
    WINbool IsOkToAdd( const Array<Vertex> &in_aopVs );
    /** Add face.
     @param aiVertId Index indices of vertices to add */
    Face AddFace( const Array<Vertex> &in_aopVs, const int id = -1 );
    /** Add face. Tries reverse order if first order fails
     @param aiVertId Index indices of vertices to add */
    Face AddFaceTryReverse( const Array<Vertex> &in_aopVs, const int id = -1 );
    ///
    void AddFacesFixOrientation( const Array< Array<Vertex> > &in_aavface );

    /// Also deletes any edges which no longer have any faces
    inline virtual void    RemoveFace          ( Face );
    /// Also removes adjacent edges and faces
    inline virtual void    RemoveVertex        ( Vertex  );
    /// Also removes adjacent faces
    inline virtual void    RemoveEdge          ( Edge );
    /// Also deletes any edges which no longer have any faces
    virtual void           RemoveFaces          ( Array<Face> & );
    /// Also removes adjacent edges and faces
    virtual void           RemoveVertices        ( Array<Vertex> & );
    /// Remove any vertices with no edges
    void RemoveDanglingVertices();
    /// Remove any vertices with no edges; return re-mapped vertices
    void RemoveDanglingVertices( Array<int> &out_aiMapVs );
    /// Does not remove vertex
    void    RemoveAdjacentFaces ( Vertex  );
    /// Blow it all away
    virtual void    Clear();

    /// Transform the mesh by the given affine transformation
    virtual void Transform( const R4Matrix &in_oFrame ) ;

	/// Subdivide mesh N times, N > 0
	///		subdivided mesh created in meshOut
	void SubdivideLoop( const int in_iNTimes, 
                              PMeshLite &out_mesh,
                        const WINbool in_bExactPts = FALSE, const WINbool in_bExactNorms = FALSE ) const;

	/// Subdivide mesh N times, N > 0
	///		subdivided mesh created in meshOut
	void SubdivideCC( const int in_iNTimes, 
                            PMeshLite &out_mesh ) const;
    /// Assume this mesh has moved; adjust locations of next level of subdivision
    void AdjustLocationsCC(       PMeshLite     & out_mesh, 
                            const Array<Vertex> & in_avMe, 
                                  Array<Vertex> & out_avNextLevel ) const;
    /// The set of vertices in the next level subdiv that are influenced by the in_av
    void InfluencesCC( const PMeshLite     & out_mesh, 
                       const Array<Vertex> & in_avMe, 
                             Array<Vertex> & out_avNextLevel ) const;
    /// The set of vertices in the next level subdiv that are influenced by the vertex
    void InfluenceCC( const PMeshLite       & in_meshChild, 
                      const Vertex            in_vChild, 
                            Array<Vertex>   & out_avParent, 
                            Array<double>   & out_adInfluence ) const;
    /// Which faces were derived from thisone
    void SubFacesCC( const PMeshLite   & in_meshNextLevel, 
                     const Array<Face> & in_afMe, 
                           Array<Face> & out_afNextLevel ) const;

    /// Make sure we can swap the edge
    WINbool IsOkToSwap( const Edge in_e ) const { return legalEdgeSwap( in_e ) ? TRUE : FALSE; }
    /// Make sure we can collapse the edge
    WINbool IsOkToCollapse( const Edge in_e ) const { return niceEdgeCollapse( in_e ) ? TRUE : FALSE; }

    ///
    void EdgeSwap( const Edge in_e ) { swapEdge(in_e); }
	///
	int SwapEdges( const Array<double> in_adDist );
	///
	void EdgeSplit( const Edge in_e ) { splitEdge(in_e); }
    ///
    void EdgeCollapse( const Edge in_e ) { collapseEdge(in_e); }    
    ///
    void EdgeCollapse( const Edge in_e, const Vertex in_v ) { collapseEdgeVertex(in_e, in_v); }
    ///
    void RemoveDupVertices( double in_dPerc = 0.01);
    //@}


    /**@name Elements of the entire mesh */
    //@{
    ///
    inline int NumVertices()  const;
    ///
    inline int NumEdges()     const;
    ///
    inline int NumFaces()     const;
    ///
    inline int AllMeshVertices ( Array<Vertex > & ) const ;
    ///
    inline int AllMeshEdges    ( Array<Edge >   & ) const ;
    ///
    inline int AllMeshFaces    ( Array<Face >   & ) const ;
    //@}

    /**@name Grow regions of elements.
      Adds to existing list, does not delete. */
    //@{
    ///
    int Grow( const Array<Vertex > &in_aopVs, 
                    Array<Edge   > &out_aopEs,
                    Array<Face   > &out_aopFs ) const;
    ///
    int Grow( const Array<Face   > &in_aopFs, 
                    Array<Vertex > &out_aopVs,
                    Array<Edge   > &out_aopEs ) const;
    ///
    int Grow( const Array<Edge   > &in_aopEs, 
                    Array<Vertex > &out_aopVs,
                    Array<Face   > &out_aopFs ) const;
    ///
    int Grow( boost::dynamic_bitset<> &io_abVs, const int in_iLoops ) const;
    ///
    int GrowRings( const Vertex &in_vertex, 
                   const int nRings,
                   Array<Vertex> &out_aopVs,
				   Array<Face>   &out_aopFs ) const;
	void allAttached(const Array<Vertex> &in_Vs, 
								boost::dynamic_bitset<> &out_Vs, 
								boost::dynamic_bitset<> &out_Es, 
								boost::dynamic_bitset<> &out_Fs);
	void allAttached(const boost::dynamic_bitset<> &in_Vs, 
					boost::dynamic_bitset<> &out_Es, 
					boost::dynamic_bitset<> &out_Fs);


    //Copy a part of the mesh and make a new mesh out of it
    //might be the same as ExtractCopy
    void CopyMesh(PMeshLite &newMesh) const;
    void MakeSubMesh(  const boost::dynamic_bitset<> &in_abFChart,
					    PMeshLite &newMesh ) const;
    ///
    int ExtractCopy( const Array<Face> &in_aopFs,
					 PMeshLite &newMesh ) const;
    ///
    int ExtractCopy( const Array<Face> &in_aopFs,
					 PMeshLite &newMesh, Array<Vertex> &out_aopVs ) const;
    //@}

    /**@name Element from IndexId()
      VertexFromIndexId( IndexId ) is an identity function. */
    //@{
    ///
    int MaxVertexId() const { return vertexnum; }
    ///
    int MaxEdgeId() const { return edgenum; }
    ///
    int MaxFaceId() const { return facenum; }
    ///
    WINbool Valid( Vertex in_v ) const { return valid(in_v) ? TRUE : FALSE; }
    ///
    WINbool Valid( Edge in_e ) const { return valid(in_e) ? TRUE : FALSE; }
    ///
    WINbool Valid( Face in_f ) const { return valid(in_f) ? TRUE : FALSE; }
    ///
    Vertex VertexFromId( int in_iId ) const { return idvertex( in_iId ); }
    ///
    Vertex VertexFromIndexId( int in_iId ) const { return idvertex( in_iId + 1 ); }
    ///
    Edge   EdgeFromId  ( int in_iId ) const { return idedge( in_iId ); }
    ///
    Edge   EdgeFromIndexId  ( int in_iId ) const { return idedge( in_iId + 1 ); }
    ///
    Face   FaceFromId  ( int in_iId ) const { return idface( in_iId ); }
    ///
    Face   FaceFromIndexId  ( int in_iId ) const { return idface( in_iId + 1 ); }
    ///
    Face   FaceFromEdges( const Edge in_e1, const Edge in_e2 ) const;
    ///
    Face   FaceFromVertices( const Vertex in_v1, const Vertex in_v2, const Vertex in_v3 ) const;
    ///
    Face   FaceFromVertices( const Vertex in_v1, const Vertex in_v2 ) const;
    //@}

    /**@name Adjacency information */
    //@{
    /**@name Vertex */
    //@{
    /// Number of adjacent edges
    int                    Degree            ( Vertex in_v ) const { return degree(in_v); }
    /// The vertices are NOT sorted. See VertexVerticesOrdered.
    inline Array< Vertex > VertexVertices    ( Vertex in_v ) const;
    /// The Edges are NOT sorted. See VertexEdgeFace.
    inline Array< Edge >   VertexEdges       ( Vertex in_v ) const;
    /// 
    inline Array< Face >   VertexFaces       ( Vertex in_v ) const;
    /// So they would lie flat if you laid them out. Requires nice vertices.
    inline Array< Vertex > VertexVerticesOrdered( Vertex in_v ) const;
    /// So they would lie flat if you laid them out. Requires nice vertices.
    /// Face i lies between edge i and edge i+1
    inline Array< Face >   VertexFacesOrdered( Vertex in_v ) const;
    /// So they would lie flat if you laid them out. Requires nice vertices.
    inline Array< Edge >   VertexEdgesOrdered( Vertex in_v ) const;
    /** So they would lie flat if you laid them out (counter clockwise orientation).
      Requires nice vertices. */
    void VertexEdgeFace(  Vertex in_v ,
                          Array< Edge > &,
                          Array< Face > & ) const;

    /// returns if there is an edge from v to w, or w to v
    inline Edge  EdgeFromVertices ( Vertex in_v,  Vertex in_w ) const { return queryEdge(in_v, in_w); }
    /// only checks for half edge v to w
    inline Edge  EdgeFromVerticesOrdered ( Vertex in_v,  Vertex in_w ) const { return queryEdgeOrdered(in_v, in_w); }
    ///
    Edge  EdgeFromFaces( Face in_f,  Face in_g ) const; 
    //@}

    /**@name Edge */
    //@{
    ///
    inline Vertex  EdgeVertex1 ( Edge in_e ) const { return vertex1( in_e ); }
    ///
    inline Vertex  EdgeOppositeVertex ( Edge in_e, Vertex in_v ) const { return oppVertex( in_v, in_e ); }
    ///
    inline Vertex  EdgeVertex2 ( Edge in_e ) const { return vertex2( in_e ); }
    /// Assumes triangular mesh; may return NULL if boundary
    inline Vertex  SideVertex1 ( Edge in_e ) const { return sideVertex1( in_e ); }
    /// Assumes triangular mesh; may return NULL if boundary
    inline Vertex  SideVertex2 ( Edge in_e ) const { return sideVertex2( in_e ); }
    ///
    inline Face    EdgeFace1   ( Edge in_e ) const { return face1( in_e ); }
    ///
    inline Face    EdgeFace2   ( Edge in_e ) const { return face2( in_e ); }
    ///
    inline Face    EdgeOtherFace   ( Edge in_e, Face in_f ) const;

    //@}

    /**@name Face */
    //@{
    /// Ordered
    inline Array< Vertex >  FaceVertices  ( Face )   const;
    /// Ordered 
    inline Array<Edge >     FaceEdges     ( Face )   const;
    /// return v0-e0-v1 where e0 goes from v0 to v1
    inline Array<Face >     FaceFaces     ( Face )   const;
    ///
    inline Edge FaceOppositeEdge( Face f, Vertex v ) const;
    /// return f0-f1-f2 where e0 goes from v0 to v1 and f0 shares e0 with in_face
    void  FaceVertexEdge ( Face,
                           Array< Vertex > &,
                           Array< Edge >   & )   const;

    int   NSides( Face in_f ) const { return numVertices(in_f); }
    //@}
    //@}

    /**@name Boundary queries */
    //@{
    ///
    int     NumBoundaries   ( Vertex v ) const { return numBoundaries( v ); }
    ///
    WINbool IsBoundary      ( Vertex v ) const { return isBoundary( v ) ? TRUE : FALSE; }
    ///
    WINbool IsMultipleBoundary( Vertex v ) const { return isNice( v ) ? TRUE : FALSE; }
    ///
    WINbool IsBoundary      ( Edge e )   const { return isBoundary( e ) ? TRUE : FALSE; }
    /// 
    void    FollowBoundary( Vertex v, Array<Vertex > & ) const;
    //@}

    /**@name Silhouette queries */
    //@{
    ///
    WINbool IsSilho(const Edge & in_e, const R3Vec & in_lookVec) const;
    ///
    void SilhouetteEdgesVertices(const R3Vec &in_lookVec, Array< Array<Vertex> > &o_aVertices, Array< Array<Edge> > &o_aEdges, Array < Array<int> > &o_adDist ) const;
    //@}

    /**@name Filtering/smoothing */
    //@{
    ///
    void FilterLocations( Array<R3Pt> &out_apt, const double dPercMove ) const;
    ///
    void Filter( const double in_dPercMove );
    ///
    void LaplacianSmoothing( const int in_iNLoops, const double in_dPercMove, const bool in_bPreserveArea );
    ///
    void LaplacianSmoothing( const int in_iNLoops, const Array<double> &in_adPercMove, const bool in_bPreserveArea );
    /// Smooth the boundary with only the boundary edges
    void SmoothBoundary( const int in_iNLoops, const double in_dPercMove );
    //@}

    /**@name Renumber and compress.
       Returns mapping of ids of elements. Called automatically before writing mesh */
    //@{
    ///
    void Renumber( );
    ///
    void Renumber( Array<int> &out_aiVsMap, Array<int> &out_aiEsMap, Array<int> &out_aiFsMap );
    /// Merge vertices, if safe. 
    void CombineTriangles( PMeshLite         & out_mesh );
    /// Returns true if had to split
    bool MakeTriangular();
    
    /// Returns true if had to split
    bool IsTriangular(){
        Array<Face> af;
        Array<Vertex> avf;
        AllMeshFaces(af);
        bool bRet = true;
        for ( int i = 0; i < af.num(); i++ ) {
            avf = FaceVertices( af[i] );
            if ( avf.num() != 3 ) {
                bRet = false;
                break;
            }
        }
        return bRet;
    };
    //@}

    /**@name Read and write */
    //@{
    ///
    void Write( const char *in_str ) const;
    ///
    void Read( const char *in_str, const WINbool in_bReverse = FALSE ) ;
    ///
    void Write( ofstream & ) const;
    ///
    void Read( ifstream &, const WINbool in_bReverse = FALSE );
    ///
    void ReadOFF( ifstream &, const WINbool in_bReverse = FALSE );
    /// A limited version of that file format (single mesh)
    void ReadObj( ifstream &, const WINbool in_bReverse = FALSE );
    ///
    void WriteBinary( ofstream & ) const;
    ///
    void ReadBinary( ifstream &, const WINbool in_bReverse = FALSE );
    ///
    void ReadInventor( ifstream & );
    /** Write the smf format (v... f...).
     @param io_strm Where to write to */
    void WriteSMF( ofstream &io_strm) const;
    /** Write the inventor format.
     @param io_strm Where to write to
     @param in_strHeader (optional) where to grab the header from */
    void WriteInventor( ofstream &io_strm,
                       const char *in_strHeader = NULL) const;
    /** Write the inventor format.
        @param in_strFName Where to write to
        @param in_strHeader (optional) where to grab the header from */
    void WriteInventor( const char *in_strFName,
                        const char *in_strHeader = NULL) const;
    /** Write the off format */
	void WriteOFF( ofstream &io_strm ) const;
	/** Write out only mesh data (vertices & faces) in obj format */
	void WriteOBJ( ofstream &io_strm ) const;
    /** Write Lu's suf format **/
    void WriteSUF( ofstream &io_strm ) const;
    //@}

    /// 
    const GMesh &Mesh() const { return *this; }
    /// 
    GMesh &Mesh()  { return *this; }

    /**@name Contours
     ** A set of contours along the given vector at the given spacing. 
      * Each array is a consecutive list of edges. There may be more than
      * one contour per z level 
      * @param in_dStartOffset is an optional offset along the normal to start slicing at (use zero for default)
      * @param in_dSpacing is how far to step along the normal between slices
      */
    void  Contours( const R3Vec &in_vecPlaneNormal,
                    const double       in_dStartOffset,
                    const double       in_dSpacing,
                    Array< Array< std::pair<Edge, double> > >  &out_aaRings ) const;
    void  Contours( const R3Vec &in_vecPlaneNormal,
                    const double       in_dStartOffset,
                    const double       in_dSpacing,
                    Array<R3Pt>  &out_apt,
                    Array<R3Vec> &out_avec) const;

   /**@name SlicingPlaneContour
     ** Contour(s) that is the intersection of the mesh with one single plane. */
    //the first contour only
    bool SlicingPlaneContour(const R3Plane &in_plane,
                             Array<R3Pt>  &out_apt) const;
    bool SlicingPlaneContour(const R3Plane &in_plane,
                                   Array < int > &out_apiVsIds,
                                   Array<double> &out_adPos) const;
    //the first contour only with priority edges
    bool SlicingPlaneContour( const R3Plane &in_plane,
                              const boost::dynamic_bitset<> &io_abPriorityEdges,
                             Array<R3Pt>  &out_apt) const;
    bool SlicingPlaneContour( const R3Plane &in_plane,
                             const boost::dynamic_bitset<> &io_abPriorityEdges,
                             Array < int > &out_apiVsIds,
                             Array<double> &out_adPos) const;
    bool SlicingPlaneContour( const R3Plane &in_plane,
                              const boost::dynamic_bitset<> &io_abPriorityEdges,
                              Array<R3Pt>  &out_apt,
                              Array < int > &out_aiEdgeIds,
                              Array<double> &out_adPos) const;
    //all contour(s)
    boost::dynamic_bitset<> SlicingPlaneContour( const R3Plane &in_plane, Array < Array<R3Pt> >  &out_aapt) const;
    boost::dynamic_bitset<> SlicingPlaneContour( const R3Plane &in_plane,
                                                 Array < Array<R3Pt> >  &out_aapt,
                                                 Array < Array < int > >  &out_aaiEdgeIds,
                                                 Array < Array<double> >  &out_aadPos) const;
    boost::dynamic_bitset<> SlicingPlaneContour( const R3Plane &in_plane,
                                                 const boost::dynamic_bitset<> &io_abPriorityEdges,
                                                 const bool in_bCheckAllEdges,
                                                 Array < Array<R3Pt> >  &out_aapt, 
                                                 Array < Array < int > >  &out_aaiEdgeIds,
                                                 Array < Array<double> >  &out_aadPos) const;
    Array<Vertex> SplitFaces(const Array < int > &in_aiEdgeIds,
                             const Array<double> &in_adPos);
    bool ClipMesh( const R3Plane &in_clippingPlane,
                   const Vertex &in_seedVertex);
    map<int, int> CutMesh(Array<Vertex> &in_avCutPath, Array<Vertex> &out_avNewVs);
    Face ReplaceFace(Face &in_faceToRemove, Array<Vertex> &in_aovNew);
    Vertex SplitFace(const int in_iFaceId, const Array<double> &in_adBary);
    
    //@}

    /**@name OpenGL stuff -- See utils/Ogl_Routines.H */
    //@{

    /// Call SetLights before calling draw. No texture, but colors if defined.
    void DrawOpenGL() ;
    /// Draws just the edges
    void DrawEdgesOpenGL() ;
    /// Call SetLights before calling draw. Uses color given.
    void DrawTransparent( const UTILSColor &, const double in_dAlpha, const WINbool in_bDoEdges ) ;
    /// Call SetLights before calling draw. Uses mesh color.
    void DrawTransparent( const float in_fAlpha ) ;

    /// 
    PMeshLite & operator=( const PMeshLite &in_oPMMesh );
    
    ///
    PMeshLite() ;
    ///
    PMeshLite( const PMeshLite &in_mesh ) { *this = in_mesh; }
    ///
    inline virtual ~PMeshLite() ;

    static void TestSubdivCC();
};

#include "Pm_MeshLite_i.H"

#endif


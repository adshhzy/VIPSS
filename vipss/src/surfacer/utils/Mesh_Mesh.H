// This may look like C code, but it is really -*- C++ -*-

// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#ifndef Mesh_h
#define Mesh_h

#include <utils/Rn_Defs.H>
#include <utils/Mesh_Map.H>
#include <utils/Mesh_Stack.H>
#include <utils/Mesh_Pqueue.H>
#include <utils/Mesh_Array.H>
#include <utils/Mesh_Pool.H>
#include <utils/Utils_Color.H>

// template<class T> class Array; // GNUG problems in Recon.cxx and G3dX.cxx
class MeshRandom;
class MHEdge; typedef MHEdge* HEdge; // private

class MVertex; 
/// \ingroup Meshes
typedef MVertex* Vertex;
class MFace; 
/// \ingroup Meshes
typedef MFace* Face;
class MEdge; 
/// \ingroup Meshes
typedef MEdge* Edge;
/// \ingroup Meshes
typedef MHEdge MCorner; 
/// \ingroup Meshes
typedef MCorner* Corner;

/** \class Mesh Mesh_Array.H utils/Mesh_Array.H
    \ingroup Meshes

  \brief Mesh: a set of Vertices, Faces, and Edges and their topological relations.

 Properties:

   - (1a)  vertices appear at most once per face
   - (1b)  an oriented edge (2 consecutive vertices of a face) appears
           in at most one face; hence the mesh is orientable.

 A Mesh must always satisfy (1a) and (1b); such a mesh is called "legal".

   - (2a) vertices are nice (contain at most 1 (possibly partial) face ring)
   - (2b) faces are nice: (a,b,c) implies no (a,c,b)

 A Mesh is "nice" if it also satisfies (2a) and (2b);
  basically, if it is everywhere a 2D manifold.

   - (3a) all faces are triangular

 A Mesh is a "nice triangular mesh" if in addition it satisfies (3a).

 MVertex allocates space for points, normals, and colors, which are used later in GMesh

 MFace allocates space for a normal, which is used later in GMesh.

 MVertex, MFace, MEdge allocate space for string, also used in GMesh
*/
class Mesh {
 public:
     /**@name Construction/copy */
     //@{
     ///
    Mesh();
    ///
    virtual ~Mesh();
    /// Completely clear all vertices, edges, faces
    void clear();
    /// copy carries flags (but not sac fields). Mesh must be empty  (does not take a GMesh!)
    void copy(const Mesh& m); 
    //@}

    /**@name  Raw manipulation functions, may lead to non-nice Meshes */
    //@{
    /// always legal, just makes a vertex and adds it to the list
    Vertex createVertex();
    /// You must remove all of the adjacent faces first. Dies if degree(v)>0
    virtual void destroyVertex(Vertex v);
    /// returns 0 if duplicate vertices or if existing edge. va must be of size 3 or bigger.
    int legalCreateFace(const Array<Vertex>& va) const;
    /// Actually create the face. Dies if !legalCreateFace()
    Face createFace(const Array<Vertex>& va);
    /// Create a triangle, optionally with the given id (which should not already exist). Otherwise,
    /// a new id is created if id == -1
    Face createFace(Vertex v1, Vertex v2, Vertex v3, int id=-1); 
    /// always legal
    virtual void destroyFace(Face f);
    /// Swap out one vertex for another.
    int substituteFaceVertex(Face f, Vertex vold, Vertex vnew);
    //@}

    /**@name Vertex operations */
    //@{
    /// Checks if the faces lie flat, .i.e, the vertex is manifold. May die if butterfly or no faces.
    int isNice(Vertex v) const;
    /// number of adjacent vertices/edges
    int degree(Vertex v) const; 
    /// A nice vertex has 0 (interior) or 1 (edge of mesh) boundaries. Bigger than that indicates a butterfly vertex
    int numBoundaries(Vertex v) const;
    /// Must be nice (calls isNice(v)), and have degree(v)>0
    int isBoundary(Vertex v) const;   
    /// f must be a triangle (isTriangle(f)). Returns the edge opposite the vertex. Slow
    Edge oppEdge(Vertex v, Face f) const; 
    /// Assumes e contains v, i.e., is (v, w); returns w
    Vertex oppVertex(Vertex v, Edge e) const;

    /**@name Elements around a vertex.
       Can get elements clockwise or counter clockwise. Checks that the vertex is nice. 
       If interior vertex, the first face/edge/vertex pulled out is random.  These functions may return NULL
       if no element exists. */
    //@{
    /// If this is an interior vertex it will just return any vertex. Otherwise, return the boundary vertex which is most clockwise. 
    Vertex mostClwVertex(Vertex v) const; 
    /// If this is an interior vertex it will just return any vertex. Otherwise, return the boundary vertex which is most counter clockwise. 
    Vertex mostCcwVertex(Vertex v) const; 

    /// Return the vertex which is clockwise from vext. May return 0 if boundary/notnice. Slow.
    Vertex clwVertex(Vertex v, Vertex vext) const; 
    /// Return the vertex which is counter clockwise from vext. May return 0 if boundary/notnice. Slow.
    Vertex ccwVertex(Vertex v, Vertex vext) const;

    /// If this is an interior vertex it will just return any face. Otherwise, return the boundary face which is most clockwise. 
    Face mostClwFace(Vertex v) const; 

    /// If this is an interior vertex it will just return any face. Otherwise, return the boundary face which is most counter clockwise. 
    Face mostCcwFace(Vertex v) const; 

    /// Return the face which is clockwise from f. May return NULL if boundary/notnice. Slow.
    Face clwFace(Vertex v, Face f) const;
    /// Return the face which is clockwise from f. May return NULL if boundary/notnice. Slow.
    Face ccwFace(Vertex v, Face f) const;

    /// If this is an interior vertex it will just return any edge. Otherwise, return the boundary edge which is most clockwise. 
    Edge mostClwEdge(Vertex v) const;
    /// If this is an interior vertex it will just return any edge. Otherwise, return the boundary edge which is most counter clockwise. 
    Edge mostCcwEdge(Vertex v) const; 

    /// Return the Edge which is clockwise from f. May return NULL if boundary/notnice. Slow.
    Edge clwEdge(Vertex v, Edge e) const;
    /// Return the Edge which is counter clockwise from f. May return NULL if boundary/notnice. Slow.
    Edge ccwEdge(Vertex v, Edge e) const; 

    /// get face which is clockwise from e relative to vertex. May return NULL if e boundary
    Face ccwFace(Vertex v, Edge e) const; 
    /// get face which is counterclockwise from e relative to vertex. May return NULL if e boundary
    Face clwFace(Vertex v, Edge e) const; 
    //@}
    //@}

    /**@name Face operations */
    //@{
    /// Edges/vertices exist
    int isNice(Face f) const;
    /// Should be >= 3
    int numVertices(Face f) const;
    /// Three vertices?
    int isTriangle(Face f) const;
    /// Has a boundary vertex or edge
    int isBoundary(Face f) const; 
    /// Returns NULL if e is a boundary edge
    Face oppFace(Face f, Edge e) const; 
    /// returns vertices in counter clockwise order
    void vertices(Face f, Array<Vertex>& va) const;
    /// Any corners on the face
    void corners(Face f, Array<Corner>& ca) const;
    /**@name move about a face */
    //@{
    /// Get the next edge after e, clockwise order
    Edge clwEdge(Face f, Edge e) const;
    /// Get the next edge after e, counter clockwise order
    Edge ccwEdge(Face f, Edge e) const;
    /// Get the next vertex after v, clockwise order (slow)
    Vertex clwVertex(Face f, Vertex v) const; 
    /// Get the next vertex after v, counter clockwise order (slow)
    Vertex ccwVertex(Face f, Vertex v) const; 
    /// Get the edge v, v1, where v1 is clockwise from v  (slow)
    Edge clwEdge(Face f, Vertex v) const;     
    /// Get the edge v, v1, where v1 is counter clockwise from v  (slow)
    Edge ccwEdge(Face f, Vertex v) const;     
    //@}
    //@}

    /**@name Edge operators */
    //@{
    /// Does this edge only have one adjacent face?
    int isBoundary(Edge e) const;
    /// Endpoint 1 (Ordering has no significance)
    Vertex vertex1(Edge e) const;
    /// Endpoint 2 (Ordering has no significance)
    Vertex vertex2(Edge e) const;
    /// Left hand face.
    Face face1(Edge e) const;
    /// Right hand face. Returns NULL if e is boundary
    Face face2(Edge e) const;
    /// Assumes that the left hand face is triangular. Returns vertex in face that is not in edge
    Vertex sideVertex1(Edge e) const; 
    /// Assumes that the right hand face is triangular. Returns vertex in face that is not in edge, or NULL if boundary
    Vertex sideVertex2(Edge e) const; 
    /// Assumes f is triangular. Returns the vertex in f that is not in e.
    Vertex oppVertex(Edge e, Face f) const;
    /// Assumes e is on the boundary. Get edge on face that is opposite the given vertex
    Edge oppBoundary(Edge e, Vertex v) const; 
    /// Follow the boundary (clockwise)
    Edge clwBoundary(Edge e) const;  
    /// Follow the boundary (counter clockwise)
    Edge ccwBoundary(Edge e) const;  
    //@}

    /**@name Corner operators
      Corners are vertices that have different normals, one for each adjacent face. Not used
      very often, but useful for tagging sharp edges/corners */
    //@{
    /// Get the corner at this vertex that corresponds to the face (if any)
    Corner corner(Vertex v, Face f) const;
    /// Get the vertex that corresponds to this corner
    Vertex corner_vertex(Corner c) const;
    /// Get the Face that corresponds to this corner
    Face corner_face(Corner c) const;
    /// around vertex
    Corner ccwCorner(Corner c) const; 
    /// around vertex
    Corner clwCorner(Corner c) const; 
    /// around face
    Corner ccwFaceCorner(Corner c) const;
    /// around face
    Corner clwFaceCorner(Corner c) const;
    /// 
    Corner ccwCorner(Vertex v, Edge e) const;
    ///
    Corner clwCorner(Vertex v, Edge e) const;
    //@} 

    /**@name Other associations */
    //@{
    /// obtain edge from vertices. Returns NULL if no such edge
    Edge queryEdge( Vertex v, Vertex w) const;
    /// obtain edge with v1 == v, v2 == w
    Edge queryEdgeOrdered( Vertex v, Vertex w) const;
    /// asserts exists
    Edge edge(Vertex v, Vertex w) const; 
    /// asserts exists, oriented
    Edge orderedEdge(Vertex v1, Vertex v2) const; 
    /// get face from 2 consecutive vertices (ccw order)
    Face face(Vertex v, Vertex vccw) const; // may ret 0
    //@}

    /**@name Counting routines (fast) */
    //@{
    ///
    int numVertices() const;
    ///
    int numFaces() const;
    ///
    int numEdges() const;
    //@}

    /**@name MeshRandom access (fast), assert there exist at least one */
    //@{
    /// Get a random vertex
    Vertex randomVertex(MeshRandom& r) const;
    /// Get a random face
    Face randomFace(MeshRandom& r) const;
    /// Get a random edge,unbiased for a closed triangle mesh
    Edge randomEdge(MeshRandom& r) const; 
    //@}

    /**@name Flags */
    //@{
    /// all flags
    enum { ALL=~0 };
    /// Get a unique id for a new flag
    static int allocateFlag() { MESHassertx(sflag<32); return 1<<(sflag++); }
    /// get selected flag(s) from the flags for the whole mesh
    int gflag(int flagmask) const; 
    /// Is this flag set on the vertex?
    int flag(Vertex v, int flagmask) const;
    /// Is this flag set on the face?
    int flag(Face f, int flagmask) const;
    /// Is this flag set on the edge?
    int flag(Edge e, int flagmask) const;
    /// set or reset selected flag(s), return previous state of flag(s)
    int gmodFlag(int flagmask, int pbool);

    /// set or reset selected flag(s), return previous state of flag(s)
    int modFlag(Vertex v, int flagmask, int pbool);
    /// set or reset selected flag(s), return previous state of flag(s)
    int modFlag(Face f, int flagmask, int pbool);
    /// set or reset selected flag(s), return previous state of flag(s)
    int modFlag(Edge e, int flagmask, int pbool);
    //@}

    /**@name Triangular mesh operations (die if not triangular!) */
    //@{
    /// would collapse be legal?
    int legalEdgeCollapse(Edge e) const;
    /// would collapse preserve a nice mesh?
    int niceEdgeCollapse(Edge e) const;
    /// would edge swap be legal?  (legal implies nice here)
    int legalEdgeSwap(Edge e) const;
    /// Do the collapse
    virtual void collapseEdgeVertex(Edge e, Vertex vs);
    /** die if !legalEdgeCollapse(e)
        remove f1,[f2], v2, (v2,{*})
       add (v1,{**})  where {**}={*}-{v1,vo1,vo2} */
    virtual void collapseEdge(Edge e);
    /** splitEdge(e) always legal
        remove f1,[f2], (v1,v2)
       add 2/4 faces, vnew, (vnew,v1), (vnew,v2), (vnew,vo1), [(vnew,vo2)] */
    virtual Vertex splitEdge(Edge e, int vid=0);
    /** die if !legalEdgeSwap(e)
       remove f1,f2, (v1,v2)
       add 2 faces, (vo1,vo2) */
    virtual Edge swapEdge(Edge e);
    //@}

    /**@name More mesh operations */
    //@{
    /// vs2 can be 0, returns v2, leaves hole
    virtual Vertex splitVertex(Vertex v1, Vertex vs1, Vertex vs2, int v2i);
    /// introduce one center vertex and triangulate face. Always legal
    virtual Vertex center_split_face(Face f); 
    /// introduce one center vertex and triangulate face. Always legal
    virtual void center_split_face(Face f, Vertex v); // always legal
    /// introduce an edge to split face on (v1,v2)
    virtual Edge split_face(Face f, Vertex v1, Vertex v2);
    /** remove the consecutive set of edges separating two faces
        (may destroy some vertices if >1 edge shared by two faces) */
    int legal_coalesce_faces(Edge e);
    /// die if !legal
    virtual Face coalesce_faces(Edge e); 
    ///
    virtual Vertex insert_vertex_on_edge(Edge e);
    ///
    virtual Edge remove_vertex_between_edges(Vertex vr); 
    //@}

    /**@name Element ids. Number from 1. Not guaranteed to be sequential */
    //@{
    /// May return Null if id does not exist
    Vertex idvertex(int id) const;
    /// Same as v->Id()
    int vertexid(Vertex v) const;
    /// May return Null if id does not exist
    Face idface(int i) const;
    /// Same as f->Id()
    int faceid(Face f) const;
    /// Same as e->Id()
    int edgeid(Edge e) const;
    /// May return Null if id does not exist
    Edge idedge( int i ) const;
    /// May return Null if id does not exist
    Vertex id_retrieve_vertex(int i) const;
    /// May return Null if id does not exist
    Face id_retrieve_face(int i) const;
    /// Are all the vertices in the mesh nice?
    int isNice() const;
    /// Call this if you want the ids to be sequential.
    void renumber();    
    /// Returns new ids for vertices, -1 if vertex not there anymore
    void renumber( Array<int> &out_aiVsMap, Array<int> &out_aiEsMap, Array<int> &out_aiFsMap );
    //@}

    /**@name  Validation. 
      Because elements are pointers, it is possible for them to "exist" but not be valid. 
      These methods check that the data pointed to exists and belongs to this mesh.
    */
    //@{
    ///
    void OK() const;
    ///
    WINbool valid(Vertex v) const;
    ///
    WINbool valid(Face f) const;
    ///
    WINbool valid(Edge e) const;
    ///
    void valid(Corner c) const;
    //@}

 public:
     /**@name Create vertices/faces with specific ids.
        Do not use this unless you really, really know what you are doing */
     //@{
    /// die if index is already used
    virtual Vertex createVertexI(int id);
    /// die if index is already used
    virtual Face createFaceI(int id, const Array<Vertex>& va);
    ///
    int nextVertexId() const { return vertexnum; }
    ///
    int nextEdgeId() const { return edgenum; }
    ///
    int nextFaceId() const { return facenum; }
    //@}

 public:       
     /**@name Should be private */
     //@{
     ///
    const MMap<int,Vertex>& fastiter_vertices() const { return id2vertex; }
     ///
    const MMap<int,Face>& fastiter_faces() const { return id2face; }
    /// isNice(v), may ret 0
    HEdge mostClwHEdge(Vertex v) const; 
    /// isNice(v), may ret 0
    HEdge mostCcwHEdge(Vertex v) const; 
    /// May return NULL
    HEdge clwHEdge(HEdge e) const;       
    /// May return NULL
    HEdge ccwHEdge(HEdge e) const;       
    //@}

 protected:
     /// 0=no, 1=min, 2=max
    static int sdebug;          
     /// id to assign to next new edge
    int edgenum;                
    /// id to assign to next new vertex
    int vertexnum;              
    /// id to assign to next new face
    int facenum;                

 private:
    int iflags;
    MMap<int,Vertex> id2vertex;  // also act as set of vertices and faces
    MMap<int,Face> id2face;
    MMap<int,Edge> id2edge;
    int nedges;
    static int sflag;
    HEdge herep(Vertex v) const;
    HEdge herep(Face v) const;
    HEdge herep(Edge e) const;
    int isBoundary(HEdge he) const;
    HEdge hefromev1(Edge e, Vertex v) const; // may ret 0
    HEdge hefromev2(Edge e, Vertex v) const; // may ret 0
    HEdge hefromef(Edge e, Face f) const;    // may ret 0
    HEdge getHEdge(Vertex v, Face f) const; // slow; on f pointing to v
    HEdge queryHEdge(Vertex v1, Vertex v2) const;
    void enterHEdge(HEdge he, Vertex v1);
    void removeHEdge(HEdge he, Vertex v1);
    void create_bogus_hedges(Array<HEdge>& ar_he);
    void remove_bogus_hedges(const Array<HEdge>& ar_he);
    void gatheredgecoalescevertices(Edge e, Array<Vertex>& va) const;
friend class MeshHEdgeIter;
friend class FaceHEdgeIter;
friend class EdgeVertexIter;
friend class EdgeFaceIter;
friend class VertexCcwHEdgeIter;
    DISABLECOPY(Mesh);
};

/** \defgroup MeshesIterators Iterators over mesh elements
  \ingroup Meshes

 Iterators can crash if continued after any change in the Mesh.
 HEdge iterators should not be used by the general public.

  Iterators are used as follows:
  
  For[Mesh/Element][Element]( Mesh, v/e/f, .. ) { <br>
       use v/e/f <br>
  }EndFor; <br>

  Don't forget the EndFor;
 */
//@{

/// private use only
#define ForVertexHEdge(m,vv,he) \
ForStack((vv)->fastiter_hedges(),HEdge,zzz) HEdge he=zzz->fastiter_prev();
/// private use only
#define ForFaceHEdge(m,ff,he) \
{ register HEdge zz_hef=(ff)->fastiter_herep(), he=zz_hef; \
  for (;he;he=he->fastiter_next(),he=he==zz_hef?0:he) {
#define DummyEndFor }}
/// private use only
#define ForVertexCcwHEdge(m,vv,he) \
{ register HEdge zz_hef=(m).mostClwHEdge(vv), he=zz_hef; \
  for (;he;he=(m).ccwHEdge(he),he=he==zz_hef?0:he) {
#define DummyEndFor }}
  
/// These mesh iterators do not define an order.
#define ForMeshVertex(m,v) ForMapValue((m).fastiter_vertices(),int,Vertex,v)
/// These mesh iterators do not define an order.
#define ForMeshFace(m,f) ForMapValue((m).fastiter_faces(),int,Face,f)
/// These mesh iterators do not define an order.
#define ForMeshEdge(m,e) \
{ Edge e; for (MeshEdgeIter zz(m);(e=zz.next());) {
#define DummyEndFor }}

/// These mesh iterators go in order (and are slower).
#define ForMeshOrderedVertex(m,v) \
{ Vertex v; for (MeshOrderedVertexIter zz(m);(v=zz.next());) {
#define DummyEndFor }}
/// These mesh iterators go in order (and are slower).
#define ForMeshOrderedFace(m,f) \
{ Face f; for (MeshOrderedFaceIter zz(m);(f=zz.next());) {
#define DummyEndFor }}
/// These mesh iterators go in order (and are slower).
#define ForMeshOrderedEdge(m,e) \
{ Edge e; for (MeshOrderedEdgeIter zz(m);(e=zz.next());) {
#define DummyEndFor }}

///  Vertex iterators do not specify order! Work correctly on non-nice vertices.
#define ForVertexVertex(m,vv,v) \
{ Vertex v; for (VertexVertexIter zz(m,vv);(v=zz.next());) {
#define DummyEndFor }}
///  Vertex iterators do not specify order! Work correctly on non-nice vertices.
#define ForVertexFace(m,vv,f) \
{ Face f; for (VertexFaceIter zz(m,vv);(f=zz.next());) {
#define DummyEndFor }}
///  Vertex iterators do not specify order! Work correctly on non-nice vertices.
#define ForVertexEdge(m,vv,e) \
{ Edge e; for (VertexEdgeIter zz(m,vv);(e=zz.next());) {
#define DummyEndFor }}
///  Vertex iterators do not specify order! Work correctly on non-nice vertices.
#define ForVertexCorner(m,vv,c) ForVertexHEdge(m,vv,c)

/// Face iterators all go CCW 
#define ForFaceVertex(m,ff,v) \
ForFaceHEdge(m,ff,zz_he) Vertex v=zz_he->fastiter_vert();
/// Face iterators all go CCW 
#define ForFaceFace(m,ff,f) \
{ Face f; for (FaceFaceIter zz(m,ff);(f=zz.next());) {
#define DummyEndFor }}
/// Face iterators all go CCW 
#define ForFaceEdge(m,ff,e) \
ForFaceHEdge(m,ff,zz_he) Edge e=zz_he->fastiter_edge();
/// Face iterators all go CCW 
#define ForFaceCorner(m,ff,c) ForFaceHEdge(m,ff,c)

/// Edge iterators do not define an order.
#define ForEdgeVertex(m,ee,v) \
{ Vertex v; for (EdgeVertexIter zz(m,ee);(v=zz.next());) {
#define DummyEndFor }}
/// Edge iterators do not define an order.
#define ForEdgeFace(m,ee,f) \
{ Face f; for (EdgeFaceIter zz(m,ee);(f=zz.next());) {
#define DummyEndFor }}

///  These vertex iterators go CCW, but require nice vertices. 
#define ForVertexCcwVertex(m,vv,v) \
{ Vertex v; for (VertexCcwVertexIter zz(m,vv);(v=zz.next());) {
#define DummyEndFor }}
///  These vertex iterators go CCW, but require nice vertices. 
#define ForVertexCcwFace(m,vv,f) \
ForVertexCcwHEdge(m,vv,zz_he) Face f=zz_he->fastiter_face();
///  These vertex iterators go CCW, but require nice vertices. 
#define ForVertexCcwCorner(m,vv,c) ForVertexCcwHEdge(m,vv,c)

//@}

//--------------------------------------------------------------------------
/** \class MVertex Mesh_mesh.H utils/Mesh_Mesh.H
    \ingroup Meshes
  \brief Vertex: position, id, normal, color

  Note: These are pool allocated, which means that the data may still look
  valid even after the vertex has been deleted.

  You should never actually use this class; use the Vertex typedef, which is
  a pointer to one of these. The mesh is responsible for allocating these.
*/
class MVertex {
 public:
     /**@name accessors */
     //@{
     ///
     const R3Pt       & GetLoc()   const { return point; }
     ///
     const UTILSColor & GetColor() const { return col; }
     ///
     const R3Vec      & GetNorm()  const { return norm; }
     ///
     void SetLoc   (const R3Pt &in_pt)        { point = in_pt; }
     ///
     void SetNorm  (const R3Vec &in_vec)      { norm = in_vec; }
     ///
     void SetColor (const UTILSColor &in_col) { col = in_col; }
     //@}

     /**@name Ids 
        Id() is always IndexId+1(); use IndexId() to get indices numbered from zero */
     //@{
     ///
     int Id() const      { return id; }
     ///
     int IndexId() const { return id-1; }
     //@}

     /// Only use if you REALLY know what you're doing
     void SetIndexId( const int in_iId ) { id = in_iId + 1; }

     ///
    static int allocateFlag() { MESHassertx(sflag<32); return 1<<(sflag++); }
 private:
    Stack<HEdge> sthe;           // hedges he such that he->prev->vert==this
    int id;
    int flags;
    const char* string;
    R3Pt  point;
    R3Vec norm;
    UTILSColor col;
    static int sflag;
friend class Mesh;
friend class GMesh;
friend class MeshHEdgeIter;
friend class VertexHEdgeIter;
 public:                        
     /// should be private
    const Stack<HEdge>& fastiter_hedges() const { return sthe; }
 public:
     /**@name Pool allocation */
     //@{
     ///
    DISABLECOPY(MVertex);
     ///
    POOLALLOCATION(MVertex);   // must be last entry of class!
     ///
    MVertex() : col(0.5,0.5,0.5) { string=0; }
     ///
    ~MVertex() { delete[] (char*)string; }
    //@}
};


/** \class MFace Mesh_mesh.H utils/Mesh_Mesh.H
    \ingroup Meshes
  \brief Face: id, normal, vertices

  Note: These are pool allocated, which means that the data may still look
  valid even after the Face has been deleted.

  You should never actually use this class; use the Face typedef, which is
  a pointer to one of these. The mesh is responsible for allocating these.
*/
class MFace {
 public:
     /// Always one more than IndexId()
     int Id() const { return id; }
     /// Use for numbering from zero
     int IndexId() const { return id-1; }
     ///
     const R3Vec &GetNorm() const { return norm; }
     ///
     void SetNorm(const R3Vec &in_vec)  { norm = in_vec; }
     /// Default flag
   static int allocateFlag() { MESHassertx(sflag<32); return 1<<(sflag++); }
 private:
    HEdge herep;
    int id;
    int flags;
    R3Vec norm;
    const char* string;
    MFace() { string=0; }
    ~MFace() { delete[] (char*)string; }
    static int sflag;
friend class Mesh;
friend class GMesh;
 public:                        
     /// should be private
    HEdge fastiter_herep() const { return herep; }
 public:
     /// must be last entry of class!
    POOLALLOCATION(MFace);     
    DISABLECOPY(MFace);
};


/** \class MEdge Mesh_mesh.H utils/Mesh_Mesh.H
    \ingroup Meshes
  \brief Edge: id, end points, faces

  Note: These are pool allocated, which means that the data may still look
  valid even after the edge has been deleted.

  You should never actually use this class; use the Edge typedef, which is
  a pointer to one of these. The mesh is responsible for allocating these.
*/
class MEdge {
 public:
     /// Always one more than IndexId()
     int Id() const { return id; }
     /// For indexing from zero
     int IndexId() const { return id-1; }

     /// Default flag
    static int allocateFlag() { MESHassertx(sflag<32); return 1<<(sflag++); }
 private:
    HEdge herep;
    int flags;
    int id;
    const char* string;
    MEdge() { string=0; }
    ~MEdge() { delete[] (char*)string; }
    static int sflag;
friend class Mesh;
friend class GMesh;
friend class MeshEdgeIter; // for efficient access to herep
 public:
     /// Pool allocation
    DISABLECOPY(MEdge);
    POOLALLOCATION(MEdge);     // must be last entry of class!
};


/** \class MHEdge Mesh_mesh.H utils/Mesh_Mesh.H
    \ingroup Meshes
  \brief Half edge

  Note: These are pool allocated, which means that the data may still look
  valid even after the edge has been deleted.

  You should never actually use this class; it is only an internal class.
*/
class MHEdge {
 private:
    HEdge prev;                 // previous HEdge in ring around face
    HEdge next;                 // next HEdge in ring around face
    HEdge sym;                  // pointer to symmetric HEdge (or 0)
    Vertex vert;                // Vertex to which this HEdge is pointing
    Face face;                  // Face on which this HEdge belongs
    Edge edge;                  // Edge to which this HEdge belongs
    const char* string;
    MHEdge() { string=0; }
    ~MHEdge() { delete[] (char*)string; }
friend class Mesh;
friend class GMesh;
friend class MeshEdgeIter;
friend class VertexHEdgeIter;
friend class VertexVertexIter;
friend class VertexFaceIter;
friend class VertexEdgeIter;
friend class FaceHEdgeIter;
friend class FaceVertexIter;
friend class FaceFaceIter;
friend class FaceEdgeIter;
friend class EdgeVertexIter;
friend class EdgeFaceIter;
friend class VertexCcwVertexIter;
friend class VertexCcwFaceIter;
 public:                        // should be private
    HEdge fastiter_next() { return next; }
    HEdge fastiter_prev() { return prev; }
    Vertex fastiter_vert() { return vert; }
    Face fastiter_face() { return face; }
    Edge fastiter_edge() { return edge; }
 public:
    DISABLECOPY(MHEdge);
    POOLALLOCATION(MHEdge);    // must be last entry of class!
};


//*** inlines

inline HEdge Mesh::herep(Vertex v) const
{
    return v->sthe.empty() ? 0 : v->sthe.top()->prev;
}

inline HEdge Mesh::herep(Face f) const { return f->herep; }
inline HEdge Mesh::herep(Edge e) const { return e->herep; }

inline int Mesh::isBoundary(Edge e) const { return herep(e)->sym ? 0 : 1; }

inline Vertex Mesh::vertex1(Edge e) const { return herep(e)->prev->vert; }
inline Vertex Mesh::vertex2(Edge e) const { return herep(e)->vert; }

inline Face Mesh::face1(Edge e) const { return herep(e)->face; }
inline Face Mesh::face2(Edge e) const {
    HEdge he=herep(e);
    return he->sym ? he->sym->face : 0;
}

inline Vertex Mesh::corner_vertex(Corner c) const { return c->vert; }
inline Face Mesh::corner_face(Corner c) const { return c->face; }

inline Corner Mesh::ccwCorner(Corner c) const { return c->sym?c->sym->prev:0; }
inline Corner Mesh::clwCorner(Corner c) const { return c->next->sym; }
inline Corner Mesh::ccwFaceCorner(Corner c) const { return c->next; }
inline Corner Mesh::clwFaceCorner(Corner c) const { return c->prev; }

inline Vertex Mesh::idvertex(int i) const { return id2vertex.get(i); }
inline int Mesh::vertexid(Vertex v) const { return v->id; }
inline Face Mesh::idface(int i) const { return id2face.get(i); }
inline int Mesh::faceid(Face f) const { return f->id; }
inline int Mesh::edgeid(Edge e) const { return e->id; }
inline Edge Mesh::idedge(int i) const { return id2edge.get(i); }

inline Vertex Mesh::id_retrieve_vertex(int i) const
{ return id2vertex.retrieve(i); }
inline Face Mesh::id_retrieve_face(int i) const
{ return id2face.retrieve(i); }

inline int Mesh::gflag(int flm) const { return iflags&flm; }
inline int Mesh::flag(Vertex v, int flm) const { return v->flags&flm; }
inline int Mesh::flag(Face f, int flm) const { return f->flags&flm; }
inline int Mesh::flag(Edge e, int flm) const { return e->flags&flm; }

inline int Mesh::gmodFlag(int flm, int pbool) {
    int st=gflag(flm);
    if (pbool) { if (st!=flm) iflags|=flm; }
    else { if (st) iflags&=~flm; }
    return st;
}
inline int Mesh::modFlag(Vertex v, int flm, int pbool) {
    int st=flag(v,flm);
    if (pbool) { if (st!=flm) v->flags|=flm; }
    else { if (st) v->flags&=~flm; }
    return st;
}
inline int Mesh::modFlag(Face f, int flm, int pbool) {
    int st=flag(f,flm);
    if (pbool) { if (st!=flm) f->flags|=flm; }
    else { if (st) f->flags&=~flm; }
    return st;
}
inline int Mesh::modFlag(Edge e, int flm, int pbool) {
    int st=flag(e,flm);
    if (pbool) { if (st!=flm) e->flags|=flm; }
    else { if (st) e->flags&=~flm; }
    return st;
}

/** \class MeshHEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Half edge iterator */
class MeshHEdgeIter {
 public:
    MeshHEdgeIter(const Mesh& m);
    HEdge next();
 private:
    MapIter<int,Vertex> mi;
    StackIter<HEdge> si;
};

/** \class MeshEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Edge iterator */
class MeshEdgeIter {
 public:
    MeshEdgeIter(const Mesh& m);
    Edge next();
 private:
    MeshHEdgeIter it;
};

/** \class MeshOrderedVertexIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex iterator */
class MeshOrderedVertexIter {
 public:
    MeshOrderedVertexIter(const Mesh& m);
    Vertex next();
 private:
    Pqueue<Vertex> pq;
};

/** \class MeshOrderedFaceIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Face iterator */
class MeshOrderedFaceIter {
 public:
    MeshOrderedFaceIter(const Mesh& m);
    Face next();
 private:
    Pqueue<Face> pq;
};

/** \class MeshOrderedEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Edge iterator */
class MeshOrderedEdgeIter {
 public:
    MeshOrderedEdgeIter(const Mesh& m);
    Edge next();
 private:
    Pqueue<Edge> pq;
};

//* Vertex

/** \class VertexHEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex iterator */
class VertexHEdgeIter {         // return HEdges pointing to v
 public:
    VertexHEdgeIter(const Mesh& m, Vertex v);
    HEdge next();
 private:
    StackIter<HEdge> si;
};

/** \class VertexVertexIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex iterator */
class VertexVertexIter {
 public:
    VertexVertexIter(const Mesh& m, Vertex v);
    Vertex next();
 private:
    VertexHEdgeIter it;
    Vertex extrav;
};

/** \class VertexFaceIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex face iterator */
class VertexFaceIter {
 public:
    VertexFaceIter(const Mesh& m, Vertex v);
    Face next();
 private:
    VertexHEdgeIter it;
};

/** \class VertexEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex edge iterator */
class VertexEdgeIter {
 public:
    VertexEdgeIter(const Mesh& m, Vertex v);
    Edge next();
 private:
    VertexHEdgeIter it;
    Edge extrae;
};

//* Face

/** \class FaceHEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Face edge iterator */
class FaceHEdgeIter {
 public:
    FaceHEdgeIter(const Mesh& m, Face f);
    HEdge next();
 private:
    HEdge hef;
    HEdge he;
};

/** \class FaceVertexIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Face vertex iterator */
class FaceVertexIter {
 public:
    FaceVertexIter(const Mesh& m, Face f);
    Vertex next();
 private:
    FaceHEdgeIter it;
};

/** \class FaceFaceIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Face faceiterator */
class FaceFaceIter {
 public:
    FaceFaceIter(const Mesh& m, Face f);
    Face next();
 private:
    FaceHEdgeIter it;
};

/** \class FaceEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Face edge iterator */
class FaceEdgeIter {
 public:
    FaceEdgeIter(const Mesh& m, Face f);
    Edge next();
 private:
    FaceHEdgeIter it;
};

//* Edge

/** \class EdgeVertexIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Edge vertex iterator */
class EdgeVertexIter {
 public:
    EdgeVertexIter(const Mesh& m, Edge e);
    Vertex next();
 private:
    HEdge he;
    int i;
};

/** \class EdgeFaceIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Edge face iterator */
class EdgeFaceIter {
 public:
    EdgeFaceIter(const Mesh& m, Edge e);
    Face next();
 private:
    HEdge he;
    int i;
};

//* VertexCcw

/** \class VertexCcwHEdgeIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex edge iterator */
class VertexCcwHEdgeIter {      // return HEdges pointing to v
 public:
    VertexCcwHEdgeIter(const Mesh& m, Vertex vp);
    HEdge next();
 private:
    const Mesh& mesh;
    Vertex v;
    HEdge hec, hef;
};

/** \class VertexCcwVertexIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex vertex iterator */
class VertexCcwVertexIter {
 public:
    VertexCcwVertexIter(const Mesh& m, Vertex vp);
    Vertex next();
 private:
    VertexCcwHEdgeIter it;
    Vertex extrav;
};

/** \class VertexCcwFaceIter Mesh_Mesh.H utils/Mesh_Mesh.H
  \ingroup MeshesIterators
  \brief Vertex face ccw iterator */
class VertexCcwFaceIter {
 public:
    VertexCcwFaceIter(const Mesh& m, Vertex v);
    Face next();
 private:
    VertexCcwHEdgeIter it;
};

//** Mesh iter

inline HEdge MeshHEdgeIter::next()
{
    if (si) return si.next();
    if (!mi) return 0;
    for (;;) {
        mi.next();
        if (!mi) return 0;
        si.reinit(mi.value()->sthe);
        if (si) return si.next();
    }
}

inline Edge MeshEdgeIter::next()
{
    for (;;) {
        HEdge he=it.next();
        if (!he) return 0;
        Edge e=he->edge;
        if (e->herep==he) return e;
    }
}

inline Vertex MeshOrderedVertexIter::next()
{
    return pq.empty() ? 0 : pq.removemin();
}

inline Face MeshOrderedFaceIter::next()
{
    return pq.empty() ? 0 : pq.removemin();
}

inline Edge MeshOrderedEdgeIter::next()
{
    return pq.empty() ? 0 : pq.removemin();
}

//** Vertex iter

inline VertexHEdgeIter::VertexHEdgeIter(const Mesh&, Vertex v)
: si(v->sthe) { }

inline HEdge VertexHEdgeIter::next() { return si ? si.next()->prev : 0; }

//** Face iter

inline FaceHEdgeIter::FaceHEdgeIter(const Mesh& m, Face f)
{ he=hef=m.herep(f); }

inline HEdge FaceHEdgeIter::next()
{
    if (!he) return 0;
    HEdge her=he;
    he=he->next;
    if (he==hef) he=0;
    return her;
}

inline FaceVertexIter::FaceVertexIter(const Mesh& m, Face f)
: it(m,f) { }

inline Vertex FaceVertexIter::next() {
    HEdge he=it.next();
    return he ? he->vert : 0;
}

#endif

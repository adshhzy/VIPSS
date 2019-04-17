// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Gmesh_h
#define Gmesh_h

#include <utils/Rn_Defs.H>
#include <utils/Rn_Polygon.H>
#include <utils/Mesh_Mesh.H>

/** \class GMesh Mesh_GMesh.H utils/Mesh_GMesh.H
  \ingroup Meshes
  \brief Add geometry to the Mesh class

 Corner data is currently not handled.

  Actual point location is stored in the MVertex class */
class GMesh : public Mesh {
 public:
     ///
    GMesh();
    ///
    ~GMesh();
    /**@name Construction/deletetion */
    //@{
    /// copy carries flags, strings, and geometry (but not sac fields). Must be empty.
    void copy(const GMesh& mo); 
    /// Merge the input mesh with this one; mvvn will return the vertex mapping if non-null
    void merge(const GMesh& mo, MMap<Vertex,Vertex>* mvvn=0);
    /// Make sure to remove adjacent faces first
    void destroyVertex(Vertex v);
    ///
    void destroyFace(Face f);
    /// do appropriate actions with geometry, SHARPE, and face strings. Keep vs (should be one of e's vertices)
    void collapseEdgeVertex(Edge e, Vertex vs);
    /// Just collapse
    void collapseEdge(Edge e);
    /// id is (optional) id of new vertex
    Vertex splitEdge(Edge e, int id=0);
    /// Assumes triangular mesh. Swap the diamond's diagonal.
    Edge swapEdge(Edge e);
    /// Split Vertex v1, creating vs1 and vs2. New vertex has id v2i
    Vertex splitVertex(Vertex v1, Vertex vs1, Vertex vs2, int v2i);
    /// Add a vertex to the center of the face; creates new vertex
    Vertex center_split_face(Face f);
    /// Add a vertex to the center of the face; use v as new vertex
    void  center_split_face(Face f, Vertex v);
    /// Split one edge of face, given by v1, v2
    Edge split_face(Face f, Vertex v1, Vertex v2);
    /// Collapse edge e
    Face coalesce_faces(Edge e);
    /// Add a new vertex on the edge
    Vertex insert_vertex_on_edge(Edge e);
    /// Undo edge split
    Edge remove_vertex_between_edges(Vertex vr);
    /// Swap vold out for vnew
    int substituteFaceVertex(Face f, Vertex vold, Vertex vnew);
    //@}

    /**@name Geometry. Can also use methods on the MVertex, MFace classes */
    //@{
    ///
    const R3Pt & point(Vertex v) const;
    ///
    void setPoint(Vertex v, const R3Pt& p);
    ///
    const R3Vec & norm(Vertex v) const;
    ///
    void setNorm(Vertex v, const R3Vec& n);
    ///
    const R3Vec & norm(Face v) const;
    ///
    void setNorm(Face v, const R3Vec& n);
    /// Pass in the polygon
    void polygon(Face f, class R3Polygon& polygon) const;
    /// Length squared
    double length2(Edge e) const;
    ///
    double length(Edge e) const;
    /// Transform the entire mesh by frame
    void transform(const R4Matrix& frame);
    //@}

    /**@name  Strings. This is one way to add additional data to the elements. 
       The data must be converted to a string, then concatenated onto the current string.
       The strings are stored in a hash table-like structure, indexed by a key. */
    //@{
    ///
    const char* string(Vertex v) const;
    ///
    const char* string(Face f) const;
    ///
    const char* string(Edge e) const;
    ///
    const char* string(Corner c) const;
    ///
    static const char* string_key(const char* str, const char* key);
    ///
    const char* corner_key(Corner c, const char* key) const; // Corner|Vertex
    /// copies string
    void setString(Vertex v, const char* s);
    ///
    void setString(Face f, const char* s);
    ///
    void setString(Edge e, const char* s);
    ///
    void setString(Corner c, const char* s);
    ///
    static const char* string_update(const char* str, const char* key,
                                     const char* val);
    ///
    void updateString(Vertex v, const char* key, const char* val);
    ///
    void updateString(Face f, const char* key, const char* val);
    ///
    void updateString(Edge e, const char* key, const char* val);
    ///
    void updateString(Corner c, const char* key, const char* val);
    //@}

    /**@name standard I/O for my meshes (see format below) */
    //@{
    /// read a whole mesh, discard comments. 
    virtual void read(istream& is);     
    /// no '\n' required
    void readline(char* s);      
    ///
    static int recognizeLine(const char* s);
    ///
    virtual void write(ostream& os) const;
    /// pos may be 0, ret old
    ostream* recordChanges(ostream* pos); 
    //@}

    /**@name Predefined Vertex,Face,Edge flag bits. First 2 are parsed when reading  */
    //@{
    /// "cusp" on Vertex
    static const int& CUSPV;    
    /// "sharp" on Edge
    static const int& SHARPE;  
    //@}

    /**@name  Discouraged; Don't use unless you know how to update the index lists */
    //@{
    ///
    Vertex createVertexI(int id);
    /// May die if face bad
    Face createFaceI(int id, const Array<Vertex>& va); 
    //@}
 private:
    ostream* ios;               // for recordChanges
    int nextOBJVert;            //the ID for the next OBJ-format vertex
    int nextOBJFace;            //the ID for the next OBJ-format face
    
    // Disable these
    GMesh &operator=( const GMesh & );
    GMesh( const GMesh & );
};

/** \class StringKeyIter Mesh_GMesh.H utils/Mesh_GMesh.H
 \ingroup MeshesHelper
 \brief For reading and writing strings from mesh file format
 
   I/O Mesh Format (Vertices and Faces must fit on one line)
   (vertex numbers begin with 1)
   Vertex vi  x y z [{other_info}]
   Face fi  vi1 vi2 ... vin [{other_info}]
   MVertex vi newx newy newz
   Ecol v1 v2
   Eswa v1 v2
   Espl v1 v2 vnew
 Example:
   Vertex 1  1.5e2 0 1.5 {normal=(0,1,0)}
   Vertex 2  0 1.5 0
   Face 1  1 2 3
   Face 2  2 3 4 5 {color=red, phong=2}
  fi may be zero, in which case a number is assigned
*/
class StringKeyIter {
public:
    StringKeyIter(const char* pstr);
    ~StringKeyIter();
    int next(const char*& kb, int& kl, const char*& vb, int& vl);
 private:
    const char* str;
    const char* s;
    // shallow copy is safe
};

/// \ingroup MeshesHelper
#define ForStringKeyValue(S,KS,KL,VS,VL) \
{ StringKeyIter zz(S); \
  const char* KS; const char* VS; \
  int KL; int VL; \
  while (zz.next(KS,KL,VS,VL)) {
#define DummyEndFor }}

//----------------------------------------------------------------------

inline const R3Pt& GMesh::point(Vertex v) const { return v->point; }
inline const R3Vec& GMesh::norm(Vertex v) const { return v->norm; }
inline const R3Vec& GMesh::norm(Face f) const { return f->norm; }

inline const char* GMesh::string(Vertex v) const { return v->string; }
inline const char* GMesh::string(Face f) const { return f->string; }
inline const char* GMesh::string(Edge e) const { return e->string; }
inline const char* GMesh::string(Corner c) const { return c->string; }
inline void GMesh::setString(Vertex v, const char* s)
{ delete[] (char*)v->string; v->string=newString(s); }
inline void GMesh::setString(Face f, const char* s)
{ delete[] (char*)f->string; f->string=newString(s); }
inline void GMesh::setString(Edge e, const char* s)
{ delete[] (char*)e->string; e->string=newString(s); }
inline void GMesh::setString(Corner c, const char* s)
{ delete[] (char*)c->string; c->string=newString(s); }

#endif

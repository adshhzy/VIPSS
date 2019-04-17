
inline
int PMeshLite::AllMeshVertices ( Array<Vertex > &out_vs ) const
{
    out_vs.need( NumVertices() );
    int i = 0;
    ForMeshOrderedVertex( *this, v ) {
        out_vs[i++] = v;
    }EndFor;

    out_vs.need(i);
    return i;
}


inline
int PMeshLite::AllMeshEdges    ( Array<Edge> &out_es ) const
{
    out_es.need( NumEdges() );
    int i = 0;
    ForMeshOrderedEdge( *this, e ) {
        out_es[i++] = e;
    }EndFor;

    out_es.need(i);
    return i;
}


inline
int PMeshLite::AllMeshFaces    ( Array<Face> &out_fs) const
{
    out_fs.need( NumFaces() );
    int i = 0;
    ForMeshOrderedFace( *this, f ) {
        out_fs[i++] = f;
    }EndFor;

    out_fs.need(i);
    return i;
}

inline 
Face PMeshLite::EdgeOtherFace   ( Edge in_e, Face in_f ) const 
{ 
    if ( EdgeFace1(in_e) == in_f ) return EdgeFace2(in_e);
    if ( EdgeFace2(in_e) == in_f ) return EdgeFace1(in_e);

    ASSERT(FALSE);
    return EdgeFace1(in_e); 
}



/*************
 * Changing topology routines
 *************/
inline
void PMeshLite::RemoveVertex( Vertex io_pm )
{
    // Remove the vertex from the Gmesh/bmesh
    destroyVertex( io_pm );

    m_bHasChanged = TRUE;
}

inline
Vertex PMeshLite::AddVertex()
{
    // Fill in with underlying mesh data
    Vertex nv = createVertex(  );
    setPoint( nv, R3Pt(0,0,0) );
    nv->SetColor( m_colDefault );

    m_bHasChanged = TRUE;

    return nv;
}

inline
Vertex PMeshLite::AddVertex(const R3Pt &in_pt)
{
    // Fill in with underlying mesh data
    Vertex nv = createVertex(  );

    setPoint( nv, in_pt );
    nv->SetColor( m_colDefault );

    m_bHasChanged = TRUE;

    return nv;  
}

inline
Vertex PMeshLite::AddVertex(const R3Pt &in_pt, const UTILSColor &in_col)
{
    // Fill in with underlying mesh data
    Vertex nv = createVertex(  );

    setPoint( nv, in_pt );
    nv->SetColor( in_col );

    m_bHasChanged = TRUE;

    return nv;  
}

inline
Vertex PMeshLite::AddVertex(const R3Pt &in_pt, const R3Vec &in_vec, const UTILSColor &in_col)
{
    // Fill in with underlying mesh data
    Vertex nv = createVertex(  );

    setPoint( nv, in_pt );
    setNorm( nv, in_vec );
    nv->SetColor( in_col );

    m_bHasChanged = TRUE;
    return nv;  
}

inline
void PMeshLite::RemoveFace( Face io_f )
{
    destroyFace( io_f );

    m_bHasChanged = TRUE;
}

inline
void PMeshLite::RemoveEdge( Edge io_e )
{
    Face opF1 = face1( io_e );
    Face opF2 = face2( io_e );

    if ( opF1 != NULL )
        RemoveFace( opF1 );

    if ( opF2 != NULL )
        RemoveFace( opF2 );
}


/**************
 * Constructor/destructors
 **************/

inline PMeshLite::~PMeshLite() 
{
    Clear();
}

/*************
 * All mesh elements access routines
 *************/
inline
int PMeshLite::NumVertices() const
{
    return numVertices();
}

inline
int PMeshLite::NumEdges()    const
{
	return numEdges();
}

inline
int PMeshLite::NumFaces()    const
{
    return numFaces();
}


/*************
 * All mesh element element access routines
 *************/

inline
Array<Vertex>   PMeshLite::VertexVertices  ( Vertex in_v ) const
{
    Array<Vertex> out_vs;
    ForVertexVertex( *this, in_v, ov ) {
        out_vs += ov;
    } EndFor;

    return out_vs;
}


inline
Array<Edge >     PMeshLite::VertexEdges  ( Vertex in_v) const
{
    Array<Edge> out_es;

    ForVertexEdge( *this, in_v, oe ) {
        out_es += oe;
    } EndFor;

    return out_es;
}

inline
Array<Edge >     PMeshLite::VertexEdgesOrdered  ( Vertex in_v) const
{
    const Array<Vertex> av = VertexVerticesOrdered(in_v);
    Array<Edge> out_es(av.num());

    for ( int i = 0; i < av.num(); i++ )
        out_es[i] = EdgeFromVertices( av[i], in_v);

    return out_es;
}

inline
Array<Face>     PMeshLite::VertexFaces     ( Vertex in_v ) const
{
    Array<Face > out_fs;
    ForVertexFace( *this, in_v, of ) {
        out_fs += of;
    } EndFor;

    return out_fs;
}


inline Edge PMeshLite::FaceOppositeEdge( Face f, Vertex v ) const
{
    ForFaceEdge( *this, f, e ) {
        if ( EdgeVertex1(e) != v && EdgeVertex2(e) != v ) {
            return e;
        }
    }EndFor;
    ASSERT(FALSE);
    return FaceEdges(f)[0];
}

inline
Array<Vertex> PMeshLite::VertexVerticesOrdered ( Vertex in_v ) const
{
    Array<Vertex> out_vs;
    ForVertexCcwVertex( *this, in_v, ov ) {
        out_vs += ov;
    } EndFor;

    return out_vs;
}



inline
Array<Face>   PMeshLite::VertexFacesOrdered     ( Vertex in_v ) const
{
    Array<Face> out_fs;
    ForVertexCcwFace( *this, in_v, of ) {
        out_fs += of;
    } EndFor;

    return out_fs;
}



inline
Array<Vertex>   PMeshLite::FaceVertices    ( Face in_f )   const
{
    Array<Vertex> out_vs;
    ForFaceVertex( *this, in_f, ov ) {
        out_vs +=  ov;
    } EndFor;

    return out_vs;
}


inline
Array<Edge>     PMeshLite::FaceEdges       ( Face in_f )   const
{
    Array<Edge> out_es;
    ForFaceEdge( *this, in_f, oe ) {
        out_es +=  oe;
    } EndFor;

    return out_es;
}

inline
Array<Face>     PMeshLite::FaceFaces       ( Face in_f )   const
{    
    Array<Face> out_fs;
    ForFaceEdge( *this, in_f, oe ) {
        if(EdgeFace1(oe) == in_f)
            out_fs += EdgeFace2(oe);
        else
            out_fs += EdgeFace1(oe);
    } EndFor;

    return out_fs;
}

inline
double PMeshLite::GetEdgeLength( const Edge in_oEdge ) const
{
    return ::Length( EdgeVertex1( in_oEdge )->GetLoc() - EdgeVertex2( in_oEdge )->GetLoc() );
}

inline
double PMeshLite::GetFaceArea( const Face in_f ) const
{
    Array<Vertex> aopV = FaceVertices(in_f);
    if ( aopV.num() != 3 ) {
        R3Polygon poly;
        GetFacePolygon( in_f, poly );
        return poly.Area();
    } else {
        const double dAng90 = cos( M_PI / 2.0 );
        const R3Vec v1v2 = UnitSafe( aopV[0]->GetLoc() - aopV[1]->GetLoc() );
        const R3Vec v2v3 = UnitSafe( aopV[1]->GetLoc() - aopV[2]->GetLoc() );
        const R3Vec v3v1 = UnitSafe( aopV[2]->GetLoc() - aopV[0]->GetLoc() );
        const double dDot1 = fabs( Dot( v3v1, v1v2 ) - dAng90 );
        const double dDot2 = fabs( Dot( v1v2, v2v3 ) - dAng90 );
        const double dDot3 = fabs( Dot( v2v3, v3v1 ) - dAng90 );

        if ( dDot1 < dDot2 && dDot1 < dDot3 ) {
            return Length( Cross( aopV[0]->GetLoc() - aopV[1]->GetLoc(), aopV[2]->GetLoc() - aopV[0]->GetLoc() ) );
        } else if ( dDot1 < dDot3 ) {
            return Length( Cross( aopV[1]->GetLoc() - aopV[2]->GetLoc(), aopV[0]->GetLoc() - aopV[1]->GetLoc() ) );
        } else {
            return Length( Cross( aopV[2]->GetLoc() - aopV[0]->GetLoc(), aopV[1]->GetLoc() - aopV[2]->GetLoc() ) );
        }
    }
    return 0.0;
}

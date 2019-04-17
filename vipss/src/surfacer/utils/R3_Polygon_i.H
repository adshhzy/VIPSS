/* -------------------------------------------------------------------------
 * DESCR   :	Accessors
 * ------------------------------------------------------------------------- */
inline const R3Pt &R3Polygon::operator[] (int in_iVert) const 
{ 
    return m_apts[in_iVert]; 
}

inline R3Pt &R3Polygon::operator[] (int in_iVert) 
{ 
    return m_apts[in_iVert]; 
}
    
// Accessor that wraps the indices around
inline const R3Pt &R3Polygon::PtWrap (int in_iVert) const 
{ 
    return m_apts.wrap(in_iVert); 
}

inline R3Pt &R3Polygon::PtWrap ( int in_iVert) 
{ 
    return m_apts.wrap(in_iVert); 
}
    
/* -------------------------------------------------------------------------
 * DESCR   :	See if the polygon is valid by checking it\'s points
 * ------------------------------------------------------------------------- */
inline
WINbool R3Polygon::Valid() const 
{
   if (Num_pts() == 0) return FALSE;
   
   for (int i = 0; i < Num_pts(); i++)
     if (RNIsValid(m_apts[i]) == FALSE) return FALSE;

   return TRUE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Constructor
 * ------------------------------------------------------------------------- */
inline R3Polygon::R3Polygon( int in_iNPts) : m_apts(in_iNPts)
{
    for (int i = 0; i < in_iNPts; i++)
        RNZero(m_apts[i]);
}

/* -------------------------------------------------------------------------
 * DESCR   :	Copy a polygon
 * ------------------------------------------------------------------------- */
inline
R3Polygon &R3Polygon::operator=(const R3Polygon &in_polygon)
{
   m_apts = in_polygon.m_apts;

   return (*this);
}


/* -------------------------------------------------------------------------
 * DESCR   :	Edges and lines
 * ------------------------------------------------------------------------- */
inline
R3Line_seg R3Polygon::Edge(int in_iV1ToV2) const
{
    return R3Line_seg( PtWrap( in_iV1ToV2 ), PtWrap( in_iV1ToV2 + 1 ) );
}

inline
R3Line R3Polygon::Line(int in_iV1ToV2) const
{
    return R3Line( PtWrap( in_iV1ToV2 ), PtWrap( in_iV1ToV2 + 1 ) );
}


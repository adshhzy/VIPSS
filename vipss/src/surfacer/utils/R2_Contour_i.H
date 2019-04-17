inline const R2Polygon & R2Contour::operator[] (int in_iPoly) const
{
    return m_apolygons[in_iPoly];
}

///
inline R2Polygon & R2Contour::operator[] (int in_iPoly) 
{
    return m_apolygons[in_iPoly];
}
///
inline WINbool R2Contour::Valid() const
{
    for (int i = 0; i < m_apolygons.num(); i++)
        if ( m_apolygons[i].Valid() == FALSE ) return FALSE;

    if ( m_apolygons.num() == 0 )
        return FALSE;

    return TRUE;
}

inline void R2Contour::SetPlane( const R3Vec & in_vec, 
                                 double        in_dZLevel )
{
    m_vecPlaneNormal = UnitSafe( in_vec );
    m_dZLevel = in_dZLevel;
}


inline R2Contour &R2Contour::operator=(const R2Contour &in_oContour)
{
    m_apolygons = in_oContour.m_apolygons;
    m_aavecNorm = in_oContour.m_aavecNorm;
    m_dZLevel = in_oContour.m_dZLevel;
    m_dWidth = in_oContour.m_dWidth;
    m_vecPlaneNormal = in_oContour.m_vecPlaneNormal;

    return *this;
}

inline R2Contour::R2Contour( )
{
    SetPlane( R3Vec(0,1,0), 0 );
    m_dWidth = 0;
}

///
inline void R2Contour::Add( const R2Polygon &in_oPoly, const Array<R3Vec> &in_avec )
{
    m_apolygons += in_oPoly;
    m_aavecNorm += in_avec;

    ASSERT( in_oPoly.Num_pts() == in_avec.num() );
}
///
inline R2Contour::R2Contour( const R2Contour &in_contour ) 
{ 
    (*this) = in_contour; 
}

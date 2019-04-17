inline WINbool R2BoxedPolygon::IsOutsideBox( double dX, double dY ) const 
{
    if ( dX < m_dWMin || dX > m_dWMax ||
         dY < m_dHMin || dY > m_dHMax ) 
        return TRUE;
    return FALSE;
}

inline WINbool R2BoxedPolygon::IsInsideBox( double dX, double dY ) const 
{
    if ( dX > m_dWMin && dX < m_dWMax && dY > m_dHMin && dY < m_dHMax ) 
        return TRUE;
    return FALSE;
}

inline
R2BoxedPolygon::R2BoxedPolygon(int in_iNPts) : R2Polygon(in_iNPts) 
{
    Set_dim();
}

inline void R2BoxedPolygon::Set_dim()
{
    if (Num_pts() <= 0) {
        m_dHMax = m_dHMin = m_dWMax = m_dWMin = 1e30;
        return;
    }

    m_dHMax = (*this)[0][1];
    m_dHMin = (*this)[0][1];
    m_dWMax = (*this)[0][0];
    m_dWMin = (*this)[0][0];
    m_ptCentroid = (*this)[0];
    for (int i = 1; i < Num_pts(); i++) {
        m_dHMax = WINmax( m_dHMax, (*this)[i][1] );
        m_dHMin = WINmin( m_dHMin, (*this)[i][1] );
        m_dWMax = WINmax( m_dWMax, (*this)[i][0] );
        m_dWMin = WINmin( m_dWMin, (*this)[i][0] );
        m_ptCentroid[0] += (*this)[i][0];
        m_ptCentroid[1] += (*this)[i][1];
    }
    m_ptCentroid[0] /= (double) Num_pts();
    m_ptCentroid[1] /= (double) Num_pts();

    R2Vec v1( (*this)[1] - (*this)[0] ), v2 ( (*this)[2] - (*this)[1] );
    R3Vec cross( Cross( R3Vec(v1[0], v1[1], 0), R3Vec(v2[0], v2[1], 0) ) );
    if ( cross[2] > 0 )
        m_bClockwise = TRUE;
    else
        m_bClockwise = FALSE;
}

inline
R2BoxedPolygon::R2BoxedPolygon( const R2Polygon &in_oPoly) : R2Polygon(in_oPoly)
{
    Set_dim();
}

inline R2BoxedPolygon &R2BoxedPolygon::operator=( const R2BoxedPolygon &in_oPoly)
{
    (*this).R2Polygon::operator=(in_oPoly);
    m_dHMax = in_oPoly.m_dHMax;
    m_dHMin = in_oPoly.m_dHMin;
    m_dWMax = in_oPoly.m_dWMax;
    m_dWMin = in_oPoly.m_dWMin;
    m_ptCentroid = in_oPoly.m_ptCentroid;
    m_bClockwise = in_oPoly.m_bClockwise;

    return *this;
}
inline R2BoxedPolygon &R2BoxedPolygon::operator=( const R2Polygon &in_oPoly)
{
    (*this).R2Polygon::operator=(in_oPoly);
    Set_dim();
    return *this;
}

inline 
R2BoxedPolygon::R2BoxedPolygon( const R2BoxedPolygon &in_oPoly) 
{
    (*this) = in_oPoly;
}

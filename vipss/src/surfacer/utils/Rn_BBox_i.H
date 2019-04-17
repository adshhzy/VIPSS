template <class Point, class Vector, class Matrix>
inline int RNBBoxTC<Point, Vector, Matrix>::Dim() const 
{ 
    return m_ptCenter.Dim(); 
}

template <class Point, class Vector, class Matrix>
const Point & RNBBoxTC<Point, Vector, Matrix>::Corner ( const int in_i ) const 
{ 
    ASSERT(in_i >= 0 && in_i < m_aptCorners.num()); 
    return m_aptCorners[in_i]; 
}

template <class Point, class Vector, class Matrix>
Point RNBBoxTC<Point, Vector, Matrix>::Scale (  ) const 
{ 
    Point pt;
    for ( int iD = 0; iD < Dim(); iD++ )
        pt[iD] = m_adExtent[iD];

    return pt; 
}

template <class Point, class Vector, class Matrix>
inline const Point & RNBBoxTC<Point, Vector, Matrix>::Center ( ) const 
{ 
    return m_ptCenter; 
}

template <class Point, class Vector, class Matrix>
inline double RNBBoxTC<Point, Vector, Matrix>::Extent ( const int in_i ) const 
{ 
    ASSERT(in_i >= 0 && in_i < Dim()); 
    return m_adExtent[in_i]; 
}

template <class Point, class Vector, class Matrix>
inline void RNBBoxTC<Point, Vector, Matrix>::ApplyTransform( const Matrix &in_mat )
{
    m_ptCenter = in_mat * m_ptCenter;

    for ( int i = 0; i < Dim(); i++ ) 
        m_adExtent[i] = 0;

    for ( FORINT i = 0; i < m_aptCorners.num(); i++ ) {
        m_aptCorners[i] = in_mat * m_aptCorners[i];
        for ( int iD = 0; iD < Dim(); iD++ ) {
            m_adExtent[iD] = WINmax( m_adExtent[iD], fabs( m_aptCorners[i][iD] - m_ptCenter[iD] ) );
        }
    }

    m_mat = in_mat * m_mat;

    WINbool bInv = FALSE;
    m_matInv = m_mat.Inverse( bInv );

    ASSERT( bInv == TRUE );
}

template <class Point, class Vector, class Matrix>
inline void RNBBoxTC<Point, Vector, Matrix>::SetExtents( const Point & in_ptExtents )
{
    SArray<int> aiIndx( Dim() );
    for ( int i = 0; i < Dim(); i++ ) {
        m_adExtent[i] = in_ptExtents[i];
        aiIndx[i] = 0;
    }

    Point ptCorner;
    for ( FORINT i = 0; i < m_aptCorners.num(); i++ ) {
        for ( int iD = 0; iD < Dim(); iD++ ) {
            if ( aiIndx[iD] == 0 )
                ptCorner[iD] = -m_adExtent[iD] / 2.0;
            else 
                ptCorner[iD] = m_adExtent[iD] / 2.0;
        }

        m_aptCorners[i] = m_mat * ptCorner;

        for ( FORINT iD = 0; iD < Dim(); iD++ ) {
            if ( aiIndx[iD] == 0 ) {
                aiIndx[iD] = 1;
                break;
            } else {
                aiIndx[iD] = 0;
            }
        }
    }
}

template <class Point, class Vector, class Matrix>
inline void RNBBoxTC<Point, Vector, Matrix>::ResetMatrix( )
{
    ApplyTransform( m_matInv );

    m_mat.MakeIdentity();
    m_matInv.MakeIdentity();
}


template <class Point, class Vector, class Matrix>
inline WINbool RNBBoxTC<Point, Vector, Matrix>::Inside( const Point &in_pt ) const
{
    const Vector vec = m_matInv * (in_pt - m_ptCenter);

    for ( int iD = 0; iD < Dim(); iD++ ) {
        if ( fabs(vec[iD]) > m_adExtent[iD] )
            return FALSE;
    }

    return TRUE;
}

template <class Point, class Vector, class Matrix>
inline WINbool RNBBoxTC<Point, Vector, Matrix>::IsBoundary( const Point &in_pt ) const
{
    const Vector vec = m_matInv * (in_pt - m_ptCenter);

    for ( int iD = 0; iD < Dim(); iD++ ) {
        if ( RNIsZero( fabs(vec[iD]), m_adExtent[iD] ) )
            return FALSE;
    }

    return TRUE;
}

template <class Point, class Vector, class Matrix>
inline Point   RNBBoxTC<Point, Vector, Matrix>::Project( const Point &in_pt ) const
{
    const Vector vec = m_matInv * (in_pt - m_ptCenter);

    for ( int iD = 0; iD < Dim(); iD++ ) {
        vec[iD] = WINminmax( vec[iD], -m_adExtent[iD], m_adExtent[iD] );
    }

    return m_ptCenter + m_mat * vec;
}

template <class Point, class Vector, class Matrix>
inline WINbool RNBBoxTC<Point, Vector, Matrix>::IntersectRay( const Point &in_pt, const Vector &in_vec, Array<Point> &out_apt ) const
{
    const Vector vecStart = m_matInv * (in_pt - m_ptCenter);
    const Vector vecDir = UnitSafe( m_matInv * in_vec );

    out_apt.need(0);

    for ( int iD = 0; iD < Dim(); iD++ ) {
        if ( RNIsZero( vecDir[iD] ) )
            continue;

        if ( vecStart[iD] < -m_adExtent[iD] && vecDir[iD] > 0 ) {
            const double dT = ( -m_adExtent[iD] - vecStart[iD] ) / vecDir[iD];
            const Point p = vecStart + vecDir * dT;

            WINbool bIsInside = TRUE;
            for ( int j = 0; j < Dim(); j++ ) {
                if ( fabs( p[j] ) > m_adExtent[j] )
                    bIsInside = FALSE;
            }

            if ( bIsInside == TRUE )
                out_apt += p;
        }
        if ( vecStart[iD] > -m_adExtent[iD] && vecDir[iD] < 0 ) {
            const double dT = ( vecStart[iD] - m_adExtent[iD] ) / vecDir[iD];
            const Point p = vecStart + vecDir * dT;

            WINbool bIsInside = TRUE;
            for ( int j = 0; j < Dim(); j++ ) {
                if ( fabs( p[j] ) > m_adExtent[j] )
                    bIsInside = FALSE;
            }

            if ( bIsInside == TRUE )
                out_apt += p;
        }
    }

    if ( out_apt.num() > 0 )
        return TRUE;

    return FALSE;
}

template <class Point, class Vector, class Matrix>
inline WINbool   RNBBoxTC<Point, Vector, Matrix>::operator==( const RNBBoxTC<Point, Vector, Matrix> &in_bbox ) const
{
    if ( Dim() != in_bbox.Dim() )
        return FALSE;

    for ( int iD = 0; iD < Dim(); iD++ )
        if ( !RNApproxEqual( m_adExtent[iD], in_bbox.m_adExtent[iD] ) )
            return FALSE;

    if ( !ApproxEqual( m_mat, in_bbox.m_mat ) )
        return FALSE;

    return TRUE;
}


template <class Point, class Vector, class Matrix>
inline RNBBoxTC<Point, Vector, Matrix> &
RNBBoxTC<Point, Vector, Matrix>::operator=(const RNBBoxTC<Point, Vector, Matrix> &in_bbox)
{
    m_mat = in_bbox.m_mat;
    m_matInv = in_bbox.m_matInv;

    m_ptCenter = in_bbox.m_ptCenter;

    for ( int iD = 0; iD < Dim(); iD++ )
        m_adExtent[iD] = in_bbox.m_adExtent[iD];

    for ( int i = 0; i < m_aptCorners.num(); i++ )
        m_aptCorners[i] = in_bbox.m_aptCorners[i];

    return *this;
}

template <class Point, class Vector, class Matrix>
inline RNBBoxTC<Point, Vector, Matrix>::RNBBoxTC( const RNBBoxTC<Point, Vector, Matrix> &in_bbox ) 
{ 
    m_adExtent.init( Dim() );

    int iD = 1;
    for ( int i = 0; i < Dim(); i++ )
        iD *= 2;

    m_aptCorners.init( iD );

    (*this) = in_bbox; 
}

template <class Point, class Vector, class Matrix>
inline RNBBoxTC<Point, Vector, Matrix>::RNBBoxTC( )
{
    RNZero(m_ptCenter);

    m_adExtent.init( Dim() );
    Point ptScl;

    int iD = 1;
    for ( int i = 0; i < Dim(); i++ ) {
        iD *= 2;
        ptScl[i] = m_adExtent[i] = 0.5;
    }

    m_mat.MakeIdentity();

    m_aptCorners.init( iD );

    SetExtents( ptScl );
}

template <class Point, class Vector, class Matrix>
inline RNBBoxTC<Point, Vector, Matrix>::RNBBoxTC( const Point &in_ptCenter, const Point &in_ptScl )
{
    m_adExtent.init( Dim() );

    int iD = 1;
    for ( int i = 0; i < Dim(); i++ ) {
        iD *= 2;
        m_adExtent[i] = 0.5;
    }

    Point pt;
    RNZero(pt);

    m_ptCenter = in_ptCenter;
    m_mat = Matrix::Translation( in_ptCenter - pt );
    m_matInv = Matrix::Translation( pt - in_ptCenter );

    m_aptCorners.init( iD );

    SetExtents( in_ptScl );
}
    //@}

    ///
    /**@name Read, write, print */
    //@{
template <class Point, class Vector, class Matrix>
inline void    RNBBoxTC<Point, Vector, Matrix>::Print() const 
{
    for ( int iD = 0; iD < m_adExtent.num(); iD++ )
        TRACE("%f ", m_adExtent[iD] );

    TRACE("\n");
    m_mat.Print();
}

template <class Point, class Vector, class Matrix>
inline void    RNBBoxTC<Point, Vector, Matrix>::Write ( ofstream &out ) const 
{
    for ( int iD = 0; iD < m_adExtent.num(); iD++ )
        out << m_adExtent[iD] << " ";
    out << m_mat;
}

template <class Point, class Vector, class Matrix>
inline WINbool RNBBoxTC<Point, Vector, Matrix>::Read  ( ifstream &in )
{
    Point ptScl;
    for ( int iD = 0; iD < m_adExtent.num(); iD++ ) {
        in >> ptScl[iD];
    }
    in >> m_mat;

    WINbool bSuc = FALSE;
    m_mat.Inverse( bSuc );
    m_matInv = m_mat.Inverse( bSuc );

    SetExtents( ptScl );

    return bSuc;
}


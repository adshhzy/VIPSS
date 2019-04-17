template <class T>
WINbool C2LFTTC<T>::IsIdentity( const T in_dEps ) const
{
    if ( RNApproxEqual( m_mat[0][0].real(), m_mat[0][0].imag(), in_dEps ) == FALSE )
        return FALSE;

    if ( RNApproxEqual( m_mat[1][1].real(), m_mat[1][1].imag(), in_dEps ) == FALSE )
        return FALSE;

    if ( RNApproxEqual( m_mat[0][0].real(), m_mat[1][1].imag(), in_dEps ) == FALSE )
        return FALSE;

    if ( RNIsZero( m_mat[0][1], in_dEps ) == FALSE )
        return FALSE;

    if ( RNIsZero( m_mat[1][0], in_dEps ) == FALSE )
        return FALSE;

    return TRUE;
}


template <class T>
void C2LFTTC<T>::SetIdentity(  ) 
{
    m_mat[0][0] = m_mat[1][1] = C2PointTC<T>(1,0);
    m_mat[0][1] = m_mat[1][0] = C2PointTC<T>(0,0);
}
/// Scale the point up (s is the scale factor)
///  s 0
///  0 1
template <class T>
void C2LFTTC<T>::SetScale( const C2PointTC<T> &in_pt ) 
{
    m_mat[0][0] = in_pt;
    m_mat[0][1] = m_mat[1][0] = C2PointTC<T>(0,0);
    m_mat[1][1] = C2PointTC<T>(1,0);
}
/// translate 
///  1 t
///  0 1
template <class T>
void C2LFTTC<T>::SetTrans( const C2PointTC<T> &in_pt ) 
{
    m_mat[0][0] = C2PointTC<T>(1,0);
    m_mat[0][1] = in_pt;
    m_mat[1][0] = C2PointTC<T>(0,0);
    m_mat[1][1] = C2PointTC<T>(1,0);
}
/// Invert matrix
///   a b ^-1 = d -b
//    c d       -c a
template <class T>
void C2LFTTC<T>::SetInvert( const C2LFTTC<T> &in_mat ) 
{
    m_mat[0][0] = in_mat.m_mat[1][1];
    m_mat[1][1] = in_mat.m_mat[0][0];
    m_mat[1][0] = in_mat.m_mat[1][0] * C2PointTC<T>(-1.0,0.0);
    m_mat[0][1] = in_mat.m_mat[0][1] * C2PointTC<T>(-1.0,0.0);
}
/// Rotation
///   1  r
///   r* 1
template <class T>
void C2LFTTC<T>::SetFlip( const C2PointTC<T> &in_pt ) 
{
    m_mat[0][0] = C2PointTC<T>(1,0);
    m_mat[0][1] = in_pt;
    m_mat[1][0] = conjg(in_pt);
    m_mat[1][1] = C2PointTC<T>(1,0);
}

/// Multiplication by a point (like multiplication of a homogeneous point by a matrix)
///   a b   z  = (a z + b) / (c z + d)
///   c d   1
template <class T>
C2PointTC<T> C2LFTTC<T>::operator *( const C2PointTC<T> &in_pt ) const 
{
    const C2PointTC<T> ptNumer = (m_mat[0][0] * in_pt + m_mat[0][1]);
    const C2PointTC<T> ptDenom = (m_mat[1][0] * in_pt + m_mat[1][1]);
    return ptNumer / ptDenom ;
}

/// Multiplication by a matrix (like a regular matrix)
template <class T>
C2LFTTC<T> C2LFTTC<T>::operator *( const C2LFTTC<T> &in_mat ) const 
{
    C2LFTTC<T> mat;
    mat.m_mat[0][0] = m_mat[0][0] * in_mat.m_mat[0][0] +
                      m_mat[0][1] * in_mat.m_mat[1][0];

    mat.m_mat[0][1] = m_mat[0][0] * in_mat.m_mat[0][1] +
                      m_mat[0][1] * in_mat.m_mat[1][1];

    mat.m_mat[1][0] = m_mat[1][0] * in_mat.m_mat[0][0] +
                      m_mat[1][1] * in_mat.m_mat[1][0];

    mat.m_mat[1][1] = m_mat[1][0] * in_mat.m_mat[0][1] +
                      m_mat[1][1] * in_mat.m_mat[1][1];
    return mat;
}

/// Comparison
template <class T>
WINbool C2LFTTC<T>::operator==( const C2LFTTC<T> &in_mat ) const 
{
    if (! ( in_mat.m_mat[0][0] == m_mat[0][0] ) ) return FALSE;
    if (! ( in_mat.m_mat[1][0] == m_mat[1][0] ) ) return FALSE;
    if (! ( in_mat.m_mat[1][1] == m_mat[1][1] ) ) return FALSE;
    if (! ( in_mat.m_mat[0][1] == m_mat[0][1] ) ) return FALSE;
    return TRUE;
}

template <class T>
WINbool C2LFTTC<T>::operator!=( const C2LFTTC<T> &in_mat ) const 
{
    if ( in_mat.m_mat[0][0] != m_mat[0][0] ) return TRUE;
    if ( in_mat.m_mat[1][0] != m_mat[1][0] ) return TRUE;
    if ( in_mat.m_mat[1][1] != m_mat[1][1] ) return TRUE;
    if ( in_mat.m_mat[0][1] != m_mat[0][1] ) return TRUE;
    return FALSE;
}


template <class T>
C2LFTTC<T>::C2LFTTC( const C2PointTC<T> &in_a00,
                     const C2PointTC<T> &in_a01,
                     const C2PointTC<T> &in_a10,
                     const C2PointTC<T> &in_a11 ) 
{
    m_mat[0][0] = in_a00;
    m_mat[0][1] = in_a01;
    m_mat[1][0] = in_a10;
    m_mat[1][1] = in_a11;
}

///
template<class T>
void C2LFTTC<T>::Write(ofstream &out) const 
{ 
    m_mat[0][0].Write(out); 
    m_mat[0][1].Write(out); 
    m_mat[1][0].Write(out); 
    m_mat[1][1].Write(out); 
}

///
template<class T>
void C2LFTTC<T>::WriteBinary(ofstream &out) const 
{ 
    for ( int i = 0; i < 2; i++ ) {
        for ( int j = 0; j < 2; j++ ) {
            ::WriteBinary(m_mat[i][j].real(), out); 
            ::WriteBinary(m_mat[i][j].imag(), out); 
        }
    }
}

///
template<class T>
void C2LFTTC<T>::Read(ifstream &in) 
{ 
    m_mat[0][0].Read(in); 
    m_mat[0][1].Read(in); 
    m_mat[1][0].Read(in); 
    m_mat[1][1].Read(in);  
}

///
template<class T>
void C2LFTTC<T>::ReadBinary(ifstream &in) 
{ 
    for ( int i = 0; i < 2; i++ ) {
        for ( int j = 0; j < 2; j++ ) {
            ::ReadBinary(m_mat[i][j].real(), in); 
            ::ReadBinary(m_mat[i][j].imag(), in); 
        }
    }
}

///
template<class T>
inline ostream& operator<<( ostream &out, const C2LFTTC<T> &in_lft )
{
    out << in_lft(0,0) << in_lft(0,1) << in_lft(1,0) << in_lft(1,1);
    return out;
}

///
template<class T>
inline istream& operator>>( istream &in, C2LFTTC<T> &out_lft )
{
    in >> out_lft(0,0) >> out_lft(0,1) >> out_lft(1,0) >> out_lft(1,1);
    return in;
}


inline WINbool
ApproxEqual(const C2LFTTC<double>& v1, const C2LFTTC<double>& v2, const double dEps = RNEpsilon_d )
{
    if ( ApproxEqual( v1(0,0), v2(0,0), dEps ) == FALSE ) return FALSE;
    if ( ApproxEqual( v1(1,0), v2(1,0), dEps ) == FALSE ) return FALSE;
    if ( ApproxEqual( v1(1,1), v2(1,1), dEps ) == FALSE ) return FALSE;
    if ( ApproxEqual( v1(0,1), v2(0,1), dEps ) == FALSE ) return FALSE;

    return TRUE;
}



template <class T>
inline
void C2PointTC<T>::Set( const T in_dRe, const T in_dIm ) 
{
    m_dRe = in_dRe;
    m_dIm = in_dIm; 
}

template <class T>
inline
void C2PointTC<T>::SetRadiusTheta( const T in_dRadius, const T in_dTheta ) 
{
    m_dRe = in_dRadius * cos( in_dTheta );
    m_dIm = in_dRadius * sin( in_dTheta ); 
}

template <class T>
inline
C2PointTC<T>& C2PointTC<T>::operator += ( const C2PointTC<T> &in_pt )
{
    (*this) = (*this) + in_pt;
    return *this;
}
///
template <class T>
inline
C2PointTC<T>& C2PointTC<T>::operator -= ( const C2PointTC<T> &in_pt )
{
    (*this) = (*this) - in_pt;
    return *this;
}
///
template <class T>
inline
C2PointTC<T>& C2PointTC<T>::operator *= ( const C2PointTC<T> &in_pt )
{
    (*this) = (*this) * in_pt;
    return *this;
}
///
template <class T>
inline
C2PointTC<T>& C2PointTC<T>::operator /= ( const C2PointTC<T> &in_pt )
{
    (*this) = (*this) / in_pt;
    return *this;
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator+( const C2PointTC<T> &in_pt ) const 
{
    return C2PointTC<T>( m_dRe + in_pt.m_dRe, m_dIm + in_pt.m_dIm );
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator+( const T &in_pt ) const 
{
    return C2PointTC<T>( m_dRe + in_pt, m_dIm );
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator-( const C2PointTC<T> &in_pt ) const 
{
    return C2PointTC<T>( m_dRe - in_pt.m_dRe, m_dIm - in_pt.m_dIm );
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator-( const T &in_pt ) const 
{
    return C2PointTC<T>( m_dRe - in_pt, m_dIm );
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator-() const 
{
    return C2PointTC<T>( -m_dRe, -m_dIm );
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator*( const C2PointTC<T> &in_pt ) const 
{
    return C2PointTC<T>( m_dRe * in_pt.m_dRe - m_dIm * in_pt.m_dIm, 
                         m_dRe * in_pt.m_dIm + m_dIm * in_pt.m_dRe );
}

template <class T>
inline
C2PointTC<T> C2PointTC<T>::operator/( const C2PointTC<T> &in_pt ) const 
{
    const T dDenom = in_pt.real() * in_pt.real() + in_pt.imag() * in_pt.imag();

    return C2PointTC<T>( (m_dRe * in_pt.m_dRe + m_dIm * in_pt.m_dIm) / dDenom, 
                         (m_dIm * in_pt.m_dRe - m_dRe * in_pt.m_dIm ) / dDenom );
}

template <class T>
inline
WINbool C2PointTC<T>::operator==( const C2PointTC<T> &in_pt ) const 
{
    if ( ( in_pt.real() == real() ) && ( in_pt.imag() == imag() ) )
        return TRUE;
    return FALSE;
}

template <class T>
inline
WINbool C2PointTC<T>::operator!=( const C2PointTC<T> &in_pt ) const 
{
    if ( ( in_pt.real() == real() ) || ( in_pt.imag() == imag() ) )
        return FALSE;
    return TRUE;
}

template <class T>
inline C2PointTC<T> conjg( const C2PointTC<T> &in_pt )
{
    return C2PointTC<T>( in_pt.real(), -in_pt.imag() );
}

template <class T>
inline T norm( const C2PointTC<T> &in_pt )
{
    return in_pt.real() * in_pt.real() + in_pt.imag() * in_pt.imag();
}

template <class T>
inline C2PointTC<T> pow( const C2PointTC<T> &in_pt, const T in_d )
{
    const T dAng = atan2( in_pt.imag(), in_pt.real() );
    const T dRad = in_pt.Length();

    const T dRadNew = pow( dRad, in_d );
    const T dAngNew = dAng * in_d;

    return C2PointTC<T>( dRadNew * cos( dAngNew ), dRadNew * sin( dAngNew ) );
}

template <class T>
inline C2PointTC<T> sqrt( const C2PointTC<T> &in_pt )
{
    return pow(in_pt, (T) 0.5);
}


extern double RNEpsilon_d;
extern float RNEpsilon_f;

template <class T>
inline WINbool
ApproxEqual(const C2PointTC<T>& v1, const C2PointTC<T>& v2, const C2PointTC<T>& in_dEps )
{
    return ( fabs( v1.real() - v2.real() ) < in_dEps.real()  &&
             fabs( v1.imag() - v2.imag() ) < in_dEps.imag() ) ? TRUE : FALSE; 
}

inline WINbool
ApproxEqual(const C2PointTC<double>& v1, const C2PointTC<double>& v2, const double dEps = RNEpsilon_d )
{
    return ( fabs( v1.real() - v2.real() ) < dEps  &&
             fabs( v1.imag() - v2.imag() ) < dEps) ? TRUE : FALSE; 
}

inline WINbool
RNIsZero(const C2PointTC<double>& pt, const double dEps = RNEpsilon_d )
{
    return ( fabs( pt.real() ) < dEps  &&
             fabs( pt.imag() ) < dEps) ? TRUE : FALSE; 
}

inline WINbool
ApproxEqual(const C2PointTC<float>& v1, const C2PointTC<float>& v2, const double dEps = RNEpsilon_f )
{
    return ( fabs( v1.real() - v2.real() ) < dEps  &&
             fabs( v1.imag() - v2.imag() ) < dEps) ? TRUE : FALSE; 
}

inline WINbool
RNIsZero(const C2PointTC<float>& pt, const float dEps = RNEpsilon_f )
{
    return ( fabs( pt.real() ) < dEps  &&
             fabs( pt.imag() ) < dEps) ? TRUE : FALSE; 
}

///
template<class T>
inline ostream& operator<<( ostream &out, const C2PointTC<T> &in_pt )
{
    out << in_pt.real() << " " << in_pt.imag() << " ";
    return out;
}

///
template<class T>
inline istream& operator>>( istream &in, C2PointTC<T> &out_pt )
{
    in >> out_pt.real() >> out_pt.imag();
    return in;
}

// =========================================================
//
//	Binary operators on affine1 classes
//
// =========================================================

/** \ingroup SimpleGeometryAlgebraic */
//@{

/// Add two vectors
template<class Coord>
inline R1VectorTC<Coord>
operator + ( const R1VectorTC<Coord>& v1, const R1VectorTC<Coord>& v2 )
{
    return R1VectorTC<Coord>( v1[0] + v2[0] );
}

/// Subtract two vectors
template<class Coord>
inline R1VectorTC<Coord>
operator - ( const R1VectorTC<Coord>& v1, const R1VectorTC<Coord>& v2 )
{
    return R1VectorTC<Coord>( v1[0] - v2[0] );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R1VectorTC<Coord>
operator * ( const Coord& s, const R1VectorTC<Coord>& v )
{
    return R1VectorTC<Coord>( v[0] * s );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R1VectorTC<Coord>
operator * ( const R1VectorTC<Coord>& v, const Coord& s )
{
    return R1VectorTC<Coord>( v[0] * s );
}

/// Divide a vector by a scaler
template<class Coord>
inline R1VectorTC<Coord>
operator / ( const R1VectorTC<Coord>& v, const Coord& s )
{
     ASSERT( !RNIsZero( s ) );

   return R1VectorTC<Coord>( v[0] / s );
}

/// Dot product of two vectors
template<class Coord>
inline Coord
Dot( const R1VectorTC<Coord>& v1, const R1VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0];
}

/// Dot product of two co vectors
template<class Coord>
inline Coord
Dot( const R1VectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0];
}

/// Dot product of a vector and a co vector
template<class Coord>
inline Coord
Dot( const R1CoVectorTC<Coord>& v1, const R1VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0];
}

/// Multiplication of a vector and a co vector
template<class Coord>
inline Coord
CCW( const R1VectorTC<Coord>& v1, const R1VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0];
}

/// Subtract two points to get a vector
template<class Coord>
inline R1VectorTC<Coord>
operator-( const R1PointTC<Coord>& p1, const R1PointTC<Coord>& p2)
{
    return R1VectorTC<Coord>( p1[0] - p2[0] );
}

/// Add two covectors
template<class Coord>
inline R1CoVectorTC<Coord>
operator + ( const R1CoVectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{
    return R1CoVectorTC<Coord>( v1[0] + v2[0] );
}

/// Subtract two covectors
template<class Coord>
inline R1CoVectorTC<Coord>
operator - ( const R1CoVectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{
    
    return R1CoVectorTC<Coord>( v1[0] - v2[0] );
}

/// Scale a covector
template<class Coord>
inline R1CoVectorTC<Coord>
operator * ( const Coord& s, const R1CoVectorTC<Coord>& v )
{
    
    return R1CoVectorTC<Coord>( v[0] * s );
}


/// Scale a covector
template<class Coord>
inline R1CoVectorTC<Coord>
operator * ( const R1CoVectorTC<Coord>& v, const Coord& s )
{
    
    return R1CoVectorTC<Coord>( v[0] * s );
}


/// Scale a covector
template<class Coord>
inline R1CoVectorTC<Coord>
operator / ( const R1CoVectorTC<Coord>& v, const Coord& s )
{
    
    return R1CoVectorTC<Coord>(v[0] / s );
}

/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R1CoVectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{    
    return v1[0] * v2[0];
}


/// Multiplication of two co vectors
template<class Coord>
inline Coord
CCW( const R1CoVectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{
    
    return v1[0] * v2[0];
}

/// Point plus vector
template<class Coord>
inline R1PointTC<Coord>
operator + ( const R1PointTC<Coord>& p, const R1VectorTC<Coord>& vDelta )
{
    
    return R1PointTC<Coord>( p[0] + vDelta[0] );
}


/// Vector plus point
template<class Coord>
inline R1PointTC<Coord>
operator + ( const R1VectorTC<Coord>& vDelta, const R1PointTC<Coord>& p )
{
    
    return R1PointTC<Coord>( p[0] + vDelta[0] );
}


/// Point plus a (negative) vector
template<class Coord>
inline R1PointTC<Coord>
operator - ( const R1PointTC<Coord>& p, const R1VectorTC<Coord>& vDelta )
{
    
    return R1PointTC<Coord>( p[0] - vDelta[0] );
}


/** \brief Linear interpolation of two points
 * If s is zero, get first point. If s is one, get second point.
 * @param p1 First point
 * @param p2 Second point
 * @param s Percentage along, zero gives first point, one gives second
 * @returns Point along line p1 + (p2 - p1) s
*/
template<class Coord>
inline R1PointTC<Coord>
Lerp( const R1PointTC<Coord>&p1, const R1PointTC<Coord>&p2, const Coord& s )
{
    
    R1VectorTC<Coord> p1p2 = p2 - p1;
    
    return R1PointTC<Coord>( p1 + s*p1p2 );
}


//@}

/** \ingroup SimpleGeometryCompare */
//@{
/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator==(const R1VectorTC<Coord>& v1, const R1VectorTC<Coord>& v2)
{
    return (v1[0] == v2[0]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator!=(const R1VectorTC<Coord>& v1, const R1VectorTC<Coord>& v2)
{
    return (v1[0] != v2[0]) ? TRUE : FALSE;
}

/// Checks if within RNEpsilon_d of each other
inline
WINbool
ApproxEqual(const R1VectorTC<double>& v1, const R1VectorTC<double>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d ) ? TRUE : FALSE; 
}

/// Checks if within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R1VectorTC<float>& v1, const R1VectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f ) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R1CoVectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{
    
    return (v1[0] == v2[0]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R1CoVectorTC<Coord>& v1, const R1CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0]) ? FALSE : TRUE;
}

/// Checks if within RNEpsilon_d of each other
inline WINbool
ApproxEqual( const R1CoVectorTC<double> & v1, const R1CoVectorTC<double> & v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d ) ? TRUE : FALSE; 
}

/// Checks if within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R1CoVectorTC<float>& v1, const R1CoVectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f ) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R1PointTC<Coord>& p1, const R1PointTC<Coord>& p2 )
{
    
    return (p1[0] == p2[0]) ? TRUE : FALSE;
}


/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R1PointTC<Coord>& p1, const R1PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0]) ? FALSE : TRUE;
}

/// Checks if within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R1PointTC<double>& p1, const R1PointTC<double>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_d ) ? TRUE : FALSE; 
}

/// Checks if within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R1PointTC<float>& p1, const R1PointTC<float>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_f ) ? TRUE : FALSE; 
}

//@}


// =========================================================
//
//	Unary operators on affine1 classes
//
// =========================================================

/** \ingroup SimpleGeometryUnary */
//@{

/// Returns non-negative number
template<class Coord>
inline Coord
Length( const R1VectorTC<Coord>& v)
{
    return (Coord) fabs(v[0]);
}

/// Length squared
template<class Coord>
inline Coord
LengthSq( const R1VectorTC<Coord>& v)
{
    return v[0]*v[0];
}

/// Return a unit length vector. Will fail if length zero
template<class Coord>
inline R1VectorTC<Coord>
Unit( const R1VectorTC<Coord>& v)
{
    return v / Length(v);
}

/** \brief "Safe" version of Unit

@param v input vector
@param out_vIfBad vector to return if v is small
@returns Unit length vector
 Returns unit length vector, or out_vIfBad  if length is small (1e-12)
*/
template<class Coord>
inline R1VectorTC<Coord>
UnitSafe( const R1VectorTC<Coord>& v, const R1VectorTC<Coord>& out_vIfBad )
{
    Coord dLength = Length( v );

    if ( fabs( dLength ) < (Coord) 1e-12 )
        return out_vIfBad;

    return v / dLength;
}

/// Return unit length vector, or vector if length is small (1e-12)
template<class Coord>
inline R1VectorTC<Coord>
UnitSafe( const R1VectorTC<Coord>& v )
{
    Coord dLength = Length( v );

    if ( fabs( dLength ) < (Coord) 1e-12 )
        return v;

    return v / dLength;
}

/// Maximum of any coordinate
template<class Coord>
inline Coord
MaximumNorm( const R1VectorTC<Coord>&v)
{
    // L-infinity norm
    return v[0];
}

/// Returns non-negative number
template<class Coord>
inline Coord
Length( const R1CoVectorTC<Coord>& v)
{
    
    return (Coord) fabs(v[0]);
}


/// Length Squared
template<class Coord>
inline Coord
LengthSq( const R1CoVectorTC<Coord>& v)
{
    
    return v[0]*v[0];
}


/// Return a unit length vector. Will fail if length zero
template<class Coord>
inline R1CoVectorTC<Coord>
Unit( const R1CoVectorTC<Coord>& v)
{
    
    return v / Length(v);
}


/// Maximum of any coordinate
template<class Coord>
inline Coord
MaximumNorm( const R1CoVectorTC<Coord>& v )
{
    /// L-infinity norm
    return v[0];
}
//@}

/** \ingroup SimpleGeometryConversion*/
//@{
/// Safe conversion
template<class Coord>
inline R1CoVectorTC<Coord>
VecToCov( const R1VectorTC<Coord>& v )
{
    return R1CoVectorTC<Coord>( v[0] );
}

/// Safe conversion
template<class Coord>
inline R1VectorTC<Coord>
CovToVec( const R1CoVectorTC<Coord>& cv )
{
    return R1VectorTC<Coord>( cv[0] );
}


//@}


// =========================================================
//
//	Unary operators on affine3 classes
//
// =========================================================

/** \ingroup SimpleGeometryUnary */
//@{

/// Returns non-negative number
template<class Coord>
inline Coord
Length( const R3VectorTC<Coord>& v )
{
    return (Coord) sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

/// Length squared
template<class Coord>
inline Coord
LengthSq( const R3VectorTC<Coord>& v )
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/// Return a unit length vector. Will fail if length zero
template<class Coord>
inline R3VectorTC<Coord>
Unit( const R3VectorTC<Coord>& v )
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
inline R3VectorTC<Coord>
UnitSafe( const R3VectorTC<Coord>& v, const R3VectorTC<Coord>& out_vIfBad )
{
    Coord dLength = Length( v );

    if ( fabs( dLength ) < (Coord) 1e-12 )
        return out_vIfBad;

    return v / dLength;
}

/// Return unit length vector, or vector if length is small (1e-12)
template<class Coord>
inline R3VectorTC<Coord>
UnitSafe( const R3VectorTC<Coord>& v )
{
    Coord dLength = Length( v );

    if ( fabs( dLength ) < (Coord) 1e-12 )
        return v;

    return v / dLength;
}

/// Maximum of any coordinate
template<class Coord>
inline Coord
MaximumNorm( const R3VectorTC<Coord>& v )
{
    // L-infinity norm
    Coord max = (Coord) fabs(v[0]);
    if ((Coord) fabs(v[1]) > max)
        max = (Coord) fabs(v[1]);
    if ((Coord) fabs(v[2]) > max)
        max = (Coord) fabs(v[2]);
    return max;
}


/// Returns non-negative number
template<class Coord>
inline Coord
Length( const R3CoVectorTC<Coord>& v )
{
    return (Coord) sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

/// Length Squared
template<class Coord>
inline Coord
LengthSq( const R3CoVectorTC<Coord>& v )
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/// Return a unit length vector. Will fail if length zero
template<class Coord>
inline R3CoVectorTC<Coord>
Unit( const R3CoVectorTC<Coord>& v )
{
    return v / Length(v);
}

/// Maximum of any coordinate
template<class Coord>
inline Coord
MaximumNorm( const R3CoVectorTC<Coord>& v )
{
    // L-infinity norm
    Coord max = (Coord) fabs(v[0]);
    if ((Coord) fabs(v[1]) > max)
        max = (Coord) fabs(v[1]);
    if ((Coord) fabs(v[2]) > max)
        max = (Coord) fabs(v[2]);
    return max;
}


//@}

/** \ingroup SimpleGeometryConversion*/
//@{

/// Safe conversion
template<class Coord>
inline R3CoVectorTC<Coord>
VecToCov( const R3VectorTC<Coord>& v )
{
    return R3CoVectorTC<Coord>( v[0], v[1], v[2] );
}

/// Safe conversion
template<class Coord>
inline R3VectorTC<Coord>
CovToVec( const R3CoVectorTC<Coord>& cv )
{
    return R3VectorTC<Coord>( cv[0], cv[1], cv[2] );
}

//@}

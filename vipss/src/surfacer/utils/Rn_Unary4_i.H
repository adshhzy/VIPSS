// =========================================================
//
//	Unary operators on affine4 classes
//
// =========================================================

/** \ingroup SimpleGeometryUnary */
//@{

/// Returns non-negative number
template<class Coord>
inline Coord
Length( const R4VectorTC<Coord>& v )
{
    return (Coord) sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3] * v[3] );
}

/// Length squared
template<class Coord>
inline Coord
LengthSq( const R4VectorTC<Coord>& v )
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3] * v[3];
}

/// Return a unit length vector. Will fail if length zero
template<class Coord>
inline R4VectorTC<Coord>
Unit( const R4VectorTC<Coord>& v )
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
inline R4VectorTC<Coord>
UnitSafe( const R4VectorTC<Coord>& v, const R4VectorTC<Coord>& out_vIfBad )
{
    Coord dLength = Length( v );

    if ( fabs( dLength ) < (Coord) 1e-12 )
        return out_vIfBad;

    return v / dLength;
}

/// Return unit length vector, or vector if length is small (1e-12)
template<class Coord>
inline R4VectorTC<Coord>
UnitSafe( const R4VectorTC<Coord>& v )
{
    Coord dLength = Length( v );

    if ( fabs( dLength ) < (Coord) 1e-12 )
        return v;

    return v / dLength;
}

/// Maximum of any coordinate
template<class Coord>
inline Coord
MaximumNorm( const R4VectorTC<Coord>& v )
{
    // L-infinity norm
    Coord max = (Coord) fabs(v[0]);
    if ((Coord) fabs(v[1]) > max)
        max = (Coord) fabs(v[1]);
    if ((Coord) fabs(v[2]) > max)
        max = (Coord) fabs(v[2]);
    if ((Coord) fabs(v[3]) > max)
        max = (Coord) fabs(v[3]);
    return max;
}

// -------------------------------------
// friend functions
// -------------------------------------

/// Returns non-negative number
template<class Coord>
tinline Coord
Length( const R4CoVectorTC<Coord>& v )
{
    return (Coord) sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3] );
}

/// Length Squared
template<class Coord>
inline Coord
LengthSq( const R4CoVectorTC<Coord>& v )
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
}

/// Return a unit length vector. Will fail if length zero
template<class Coord>
inline R4CoVectorTC<Coord>
Unit( const R4CoVectorTC<Coord>& v )
{
    return v / Length(v);
}

/// Maximum of any coordinate
template<class Coord>
inline Coord
MaximumNorm( const R4CoVectorTC<Coord>& v )
{
    // L-infinity norm
    Coord max = (Coord) fabs(v[0]);
    if ((Coord) fabs(v[1]) > max)
        max = (Coord) fabs(v[1]);
    if ((Coord) fabs(v[2]) > max)
        max = (Coord) fabs(v[2]);
    if ((Coord) fabs(v[3]) > max)
        max = (Coord) fabs(v[3]);
    return max;
}
//@}

/** \ingroup SimpleGeometryConversion*/
//@{

/// Safe conversion
template<class Coord>
inline R4CoVectorTC<Coord>
VecToCov( const R4VectorTC<Coord>& v )
{
    return R4CoVectorTC<Coord>( v[0], v[1], v[2], v[3] );
}

/// Safe conversion
template<class Coord>
inline R4VectorTC<Coord>
CovToVec( const R4CoVectorTC<Coord>& cv )
{
    return R4VectorTC<Coord>( cv[0], cv[1], cv[2], cv[3] );
}


//@}


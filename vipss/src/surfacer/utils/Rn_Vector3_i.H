// =========================================================
//
//	Methods for R3VectorTC template class
//
// =========================================================


// -------------------------------------
// constructors
// -------------------------------------

template<class Coord>
inline
R3VectorTC<Coord>::R3VectorTC( const Coord& _dx, const Coord& _dy, const Coord& _dz )
{
    dx = _dx; dy = _dy; dz = _dz;
}


// -------------------------------------
// assignment operators
// -------------------------------------

template<class Coord>
inline R3VectorTC<Coord>&
R3VectorTC<Coord>::operator += ( const R3VectorTC<Coord>& v)
{
    dx += v.dx;
    dy += v.dy;
    dz += v.dz;
    return *this;
}

template<class Coord>
inline R3VectorTC<Coord>&
R3VectorTC<Coord>::operator -= ( const R3VectorTC<Coord>& v )
{
    dx -= v.dx;
    dy -= v.dy;
    dz -= v.dz;
    return *this;
}

template<class Coord>
inline R3VectorTC<Coord>&
R3VectorTC<Coord>::operator *= ( const Coord& s )
{
    dx *= s;
    dy *= s;
    dz *= s;
    return *this;
}

template<class Coord>
inline R3VectorTC<Coord>&
R3VectorTC<Coord>::operator /= ( const Coord& s )
{
    const Coord s1 = (Coord) 1.0 / s;
    dx *= s1;
    dy *= s1;
    dz *= s1;
    return *this;
}

// -------------------------------------
// unary operators
// -------------------------------------

template<class Coord>
inline R3VectorTC<Coord>
R3VectorTC<Coord>::operator + () const
{
    return *this;
}

template<class Coord>
inline R3VectorTC<Coord>
R3VectorTC<Coord>::operator - () const
{
    return R3VectorTC<Coord>( -dx, -dy, -dz );
}


// -------------------------------------
// miscellaneous methods
// -------------------------------------
/// Returns a vector which is perpendicular (not unique)
template<class Coord>
inline R3VectorTC<Coord>  R3VectorTC<Coord>::Perpendicular() const
{
    if ( fabs( (*this)[0] ) >= fabs( (*this)[1] ) && fabs( (*this)[0] ) >= fabs( (*this)[2] ) ) {
        if ( fabs( (*this)[1] ) > fabs( (*this)[2] ) ) {
            return R3VectorTC<Coord>( -(*this)[1], (*this)[0], 0 );
        } else {
            return R3VectorTC<Coord>( -(*this)[2], 0, (*this)[0] );
        }
    } else if ( fabs( (*this)[1] ) >= fabs( (*this)[2] ) ) {
        if ( fabs( (*this)[0] ) > fabs( (*this)[2] ) ) {
            return R3VectorTC<Coord>( -(*this)[1], (*this)[0], 0 );
        } else {
            return R3VectorTC<Coord>( 0, -(*this)[1], (*this)[2] );
        }
    }
    if ( fabs( (*this)[1] ) > fabs( (*this)[0] ) ) {
        return R3VectorTC<Coord>( 0, (*this)[2], -(*this)[1] );
    } else {
        return R3VectorTC<Coord>( (*this)[2], 0, -(*this)[0] );
    }

    ASSERT(FALSE);

    return R3VectorTC<Coord>(0,0,0);
}


// -------------------------------------
// Self-editing functions
// -------------------------------------
template<class Coord>
inline Coord R3VectorTC<Coord>::Normalize()
{
    const Coord dLen = sqrt(dx * dx + dy * dy + dz * dz);
    if ( ! RNIsZero( dLen ) ) {
        dx = dx / dLen;
        dy = dy / dLen;
        dz = dz / dLen;
    }

    return dLen;
}


// -------------------------------------
// Read/write/print functions
// -------------------------------------
template<class Coord>
inline void R3VectorTC<Coord>::Write(ofstream &out) const
{
    out << dx << " " << dy << " " << dz << " ";
}

template<class Coord>
inline void R3VectorTC<Coord>::WriteBinary(ofstream &out) const
{
    out.write( (const char *) &dx, Dim() * sizeof(Coord) );
}

template<class Coord>
inline WINbool R3VectorTC<Coord>::Read(ifstream &in) 
{

    in >> dx >> dy >> dz;

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline WINbool R3VectorTC<Coord>::ReadBinary(ifstream &in) 
{
    in.read( (char *) &dx, Dim() * sizeof(Coord) );

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline void R3VectorTC<Coord>::Print( const WINbool in_bDoReturn ) const
{
    cout << dx << " " << dy << " " << dz;
    if ( in_bDoReturn == TRUE )
        cout << "\n";
    else
        cout << " ";
}




// =========================================================
//
//	Methods for R3PointTC template class
//
// =========================================================


// -------------------------------------
// constructors
// -------------------------------------

template<class Coord>
inline
R3PointTC<Coord>::R3PointTC ( const Coord& _x, const Coord& _y, const Coord& _z )
{
    x = _x; y = _y; z = _z;
}



// -------------------------------------
// constructors
// -------------------------------------

template<class Coord>
inline R3PointTC<Coord>&
R3PointTC<Coord>::operator += ( const R3VectorTC<Coord>& vDelta )
{
    x += vDelta[0];
    y += vDelta[1];
    z += vDelta[2];
    return *this;
}

template<class Coord>
inline R3PointTC<Coord>&
R3PointTC<Coord>::operator -= (const R3VectorTC<Coord>& vDelta )
{
    x -= vDelta[0];
    y -= vDelta[1];
    z -= vDelta[2];
    return *this;
}

// -------------------------------------
// binary operators
// -------------------------------------



// -------------------------------------
// point dominance
// -------------------------------------

template<class Coord>
inline WINbool
R3PointTC<Coord>::operator < ( const R3PointTC<Coord>& p ) const
{
    return (x < p.x && y < p.y && z < p.z) ? TRUE : FALSE;
}

template<class Coord>
inline WINbool
R3PointTC<Coord>::operator<= ( const R3PointTC<Coord>& p ) const
{
    return (x <= p.x && y <= p.y && z <= p.z) ? TRUE : FALSE;
}


// -------------------------------------
// Read/write/print functions
// -------------------------------------
template<class Coord>
inline void R3PointTC<Coord>::Write(ofstream &out) const
{
    out << x << " " << y << " " << z << " ";
}

template<class Coord>
inline void R3PointTC<Coord>::WriteBinary(ofstream &out) const
{
    out.write( (const char *) &x, Dim() * sizeof(Coord) );
}

template<class Coord>
inline WINbool R3PointTC<Coord>::Read(ifstream &in) 
{
    in >> x >> y >> z;

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline WINbool R3PointTC<Coord>::ReadBinary(ifstream &in) 
{
    in.read( (char *) &x, Dim() * sizeof(Coord) );

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline void R3PointTC<Coord>::Print( const WINbool in_bDoReturn ) const
{
    TRACE("(%f, %f, %f)", x,y,z);

    if ( in_bDoReturn == TRUE )
        TRACE("\n");
    else
        TRACE("\n");
}

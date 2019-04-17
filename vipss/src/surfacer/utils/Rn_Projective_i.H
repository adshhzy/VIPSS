
// =========================================================
// inline function implementations
// =========================================================

/** \ingroup ProjectiveGeometry 
 \defgroup ProjectiveGeometryAlgebraic Algebraic operators

Adding, subtracting, multiplying, dot products
*/
//@{
///
template<class Coord>
inline Coord
Dot( const R3HPointTC<Coord>& p1, const R3HPlaneTC<Coord>& p2 )
{
    return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z + p1.w*p2.w;
}

///
template<class Coord>
inline Coord
Dot( const R3HPlaneTC<Coord>& p1, const R3HPointTC<Coord>& p2 )
{
    return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z + p1.w*p2.w;
}

///
template<class Coord>
inline
Coord
SignedDistance( const R3HLineTC<Coord>& a, const R3HLineTC<Coord>& b )
{
    return  a.pl[0]*b.pl[4] +
            a.pl[1]*b.pl[5] +
            a.pl[2]*b.pl[3] +
            a.pl[3]*b.pl[0] +
            a.pl[4]*b.pl[1] +
            a.pl[5]*b.pl[2];
}

//@}

template<class Coord>
inline
R2HPointTC<Coord>::R2HPointTC( const R2PointTC<Coord>& p )
{
    u = p.u; v = p.v; w = 1.0f;
}

template<class Coord>
inline
R2HPointTC<Coord>::operator R2PointTC<Coord> ( )
{
    return R2PointTC<Coord>( u / w, v / w );
}

template<class Coord>
inline
R3HPointTC<Coord>::R3HPointTC( const Coord& _x, const Coord& _y,
                               const Coord& _z, const Coord& _w)
{
    x = _x; y = _y; z = _z; w = _w;
}

template<class Coord>
inline
R3HPointTC<Coord>::R3HPointTC( const R3PointTC<Coord>& p )
{
    x = p.x; y = p.y; z = p.z; w = 1.0f;
}

template<class Coord>
inline
R3HPointTC<Coord>::operator R3PointTC<Coord> ( )
{
    return R3PointTC<Coord>( x / w, y / w, z / w );
}


template<class Coord>
inline
R3HPlaneTC<Coord>::R3HPlaneTC( const Coord& _x, const Coord& _y,
                               const Coord& _z, const Coord& _w)
{
    x = _x; y = _y; z = _z; w = _w;
}

template<class Coord>
inline
R3HPlaneTC<Coord>::R3HPlaneTC( const R3PointTC<Coord>& p1,
                               const R3PointTC<Coord>& p2,
                               const R3PointTC<Coord>& p3)
{
    R3VectorTC<Coord> v2 = p2 - p1;
    R3VectorTC<Coord> v3 = p3 - p1;
    x = v2.dy*v3.dz-v2.dz*v3.dy;
    y = v2.dz*v3.dx-v2.dx*v3.dz;
    z = v2.dx*v3.dy-v2.dy*v3.dx;
    w = -p1.x*x - p1.y*y - p1.z*z;
}

template<class Coord>
inline
R3HLineTC<Coord>::R3HLineTC( const R3HPointTC<Coord>& p, const R3HPointTC<Coord>& q )
{
    pl[0] = p.x*q.y - q.x*p.y;
    pl[1] = p.x*q.z - q.x*p.z;
    pl[2] = p.x*q.w - q.x*p.w;
    pl[3] = p.y*q.z - q.y*p.z;
    pl[4] = p.z*q.w - q.z*p.w;
    pl[5] = q.y*q.w - p.y*p.w;
}

///
template<class Coord>
inline
R3HLineTC<Coord>::R3HLineTC( const R3HPlaneTC<Coord>& p, const R3HPlaneTC<Coord>& q )
{
    pl[0] = p.x*q.y - q.x*p.y;
    pl[1] = p.x*q.z - q.x*p.z;
    pl[2] = p.x*q.w - q.x*p.w;
    pl[3] = p.y*q.z - q.y*p.z;
    pl[4] = p.z*q.w - q.z*p.w;
    pl[5] = q.y*q.w - p.y*p.w;
}


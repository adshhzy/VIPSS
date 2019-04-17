// =========================================================
//
//	Binary operators on affine4 classes
//
// =========================================================

/** \ingroup SimpleGeometryAlgebraic */
//@{

/// Add two vectors
template<class Coord>
inline R4VectorTC<Coord>
operator + ( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    return R4VectorTC<Coord>( v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2], v1[3] + v2[3] );
}

/// Subtract two vectors
template<class Coord>
inline R4VectorTC<Coord>
operator - ( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    return R4VectorTC<Coord>( v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2], v1[3] - v2[3] );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R4VectorTC<Coord>
operator * ( const Coord& s, const R4VectorTC<Coord>& v )
{
    return R4VectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s, v[3] * s );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R4VectorTC<Coord>
operator * ( const R4VectorTC<Coord>& v, const Coord& s )
{
    return R4VectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s, v[3] * s );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R4VectorTC<Coord>
operator / ( const R4VectorTC<Coord>& v, const Coord& s )
{
    ASSERT( !RNIsZero( s ) );

    Coord s1 = 1.0f / s;
    return R4VectorTC<Coord>( v[0] * s1, v[1] * s1, v[2] * s1, v[3] * s1 );
}


/// Dot product of two vectors
template<class Coord>
inline Coord
Dot( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}

/// Dot product of a vector and a co vector
template<class Coord>
inline Coord
Dot( const R4VectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}

/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R4CoVectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}

/// Outer product of two vectors
template<class Coord>
inline R4MatrixTC<Coord>
OuterProduct( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of a covector and a vector
template<class Coord>
inline R4MatrixTC<Coord>
OuterProduct( const R4CoVectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of a covector and a vector
template<class Coord>
inline R4MatrixTC<Coord>
OuterProduct( const R4VectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of two covectors
template<class Coord>
inline R4MatrixTC<Coord>
OuterProduct( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Alternating product of two vectors
template<class Coord>
inline R4MatrixTC<Coord>
AlternatingProduct( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   ret(1,2) = -ret(2,1);

    return ret;
}

/// Alternating product of a covector and a vector
template<class Coord>
inline R4MatrixTC<Coord>
AlternatingProduct( const R4VectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = ret(3,3) = 0.0f;

    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; 
    ret(3,0) = v1[3]*v2[0]-v1[0]*v2[3]; 

    ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(3,1) = v1[3]*v2[1]-v1[1]*v2[3];

    ret(3,2) = v1[3]*v2[2]-v1[2]*v2[3];

    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   
    ret(0,3) = -ret(3,0);
    
    ret(1,2) = -ret(2,1);
    ret(1,3) = -ret(3,1);

    ret(3,2) = -ret(2,3);

    return ret;
}

/// Alternating product of a covector and a vector
template<class Coord>
inline R4MatrixTC<Coord>
AlternatingProduct( const R4CoVectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = ret(3,3) = 0.0f;

    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; 
    ret(3,0) = v1[3]*v2[0]-v1[0]*v2[3]; 

    ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(3,1) = v1[3]*v2[1]-v1[1]*v2[3];

    ret(3,2) = v1[3]*v2[2]-v1[2]*v2[3];

    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   
    ret(0,3) = -ret(3,0);
    
    ret(1,2) = -ret(2,1);
    ret(1,3) = -ret(3,1);

    ret(3,2) = -ret(2,3);

    return ret;
}

/// Alternating product of two covectors
template<class Coord>
inline R4MatrixTC<Coord>
AlternatingProduct( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    R4MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = ret(3,3) = 0.0f;

    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; 
    ret(3,0) = v1[3]*v2[0]-v1[0]*v2[3]; 

    ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(3,1) = v1[3]*v2[1]-v1[1]*v2[3];

    ret(3,2) = v1[3]*v2[2]-v1[2]*v2[3];

    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   
    ret(0,3) = -ret(3,0);
    
    ret(1,2) = -ret(2,1);
    ret(1,3) = -ret(3,1);

    ret(3,2) = -ret(2,3);

    return ret;
}


// -------------------------------------
// binary operators
// -------------------------------------

/// Add two covectors
template<class Coord>
inline R4CoVectorTC<Coord>
operator + ( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    return R4CoVectorTC<Coord>( v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2], v1[3] + v2[3] );
}

/// Subtract two covectors
template<class Coord>
inline R4CoVectorTC<Coord>
operator - ( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    return R4CoVectorTC<Coord>( v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2], v1[3] - v2[3] );
}

/// Scale a covector
template<class Coord>
inline R4CoVectorTC<Coord>
operator * ( const Coord& s, const R4CoVectorTC<Coord>& v )
{
    return R4CoVectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s, v[3] * s );
}

/// Scale a covector
template<class Coord>
inline R4CoVectorTC<Coord>
operator * ( const R4CoVectorTC<Coord>& v, const Coord& s )
{
    return R4CoVectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s, v[3] * s );
}

/// Scale a covector
template<class Coord>
inline R4CoVectorTC<Coord>
operator / ( const R4CoVectorTC<Coord>& v, const Coord& s )
{
    ASSERT( !RNIsZero( s ) );

    Coord s1 = Coord(1.0f) / s;
    return R4CoVectorTC<Coord>( v[0] * s1, v[1] * s1, v[2] * s1, v[3] * s1 );
}

/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3];
}


/// Subtract two points
template<class Coord>
inline R4VectorTC<Coord>
operator - ( const R4PointTC<Coord>& p1, const R4PointTC<Coord>& p2 )
{
    return R4VectorTC<Coord>( p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2], p1[3] - p2[3] );
}

/// Add a vector to a point
template<class Coord>
inline R4PointTC<Coord>
operator + ( const R4PointTC<Coord>& p, const R4VectorTC<Coord>& vDelta )
{
    return R4PointTC<Coord>( p[0] + vDelta[0], p[1] + vDelta[1], p[2] + vDelta[2], p[3] + vDelta[3] );
}

/// Add a covector to a point
template<class Coord>
inline R4PointTC<Coord>
operator + ( const R4CoVectorTC<Coord>& vDelta, const R4PointTC<Coord>& p )
{
    return R4PointTC<Coord>( p[0] + vDelta[0], p[1] + vDelta[1], p[2] + vDelta[2], p[3] + vDelta[3] );
}

/// Add a negative vector to a point
template<class Coord>
inline R4PointTC<Coord>
operator - ( const R4PointTC<Coord>& p, const R4VectorTC<Coord>& vDelta )
{
    return R4PointTC<Coord>( p[0] - vDelta[0], p[1] - vDelta[1], p[2] - vDelta[2], p[3] - vDelta[3] );
}

// -------------------------------------
// friend functions
// -------------------------------------

/** \brief Linear interpolation of two points
 * If s is zero, get first point. If s is one, get second point.
 * @param p1 First point
 * @param p2 Second point
 * @param s Percentage along, zero gives first point, one gives second
 * @returns Point along line p1 + (p2 - p1) s
*/
template<class Coord>
inline R4PointTC<Coord>
Lerp( const R4PointTC<Coord>&p1, const R4PointTC<Coord>&p2, Coord s )
{
    R4VectorTC<Coord> p1p2 = p2 - p1;
    return R4PointTC<Coord>( p1 + s*p1p2 );
}


/// Add two matrices
template<class Coord>
inline R4MatrixTC<Coord>
operator + ( const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2 )
{
    R4MatrixTC<Coord> ret;

    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = m1(i, j) + m2(i, j);
    return ret;
}

/// Subtract two matrics
template<class Coord>
inline R4MatrixTC<Coord>
operator - ( const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2 )
{
    R4MatrixTC<Coord> ret;

    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = m1(i, j) - m2(i, j);
    return ret;
}

/// Multiply a matrix by a scaler
template<class Coord>
inline R4MatrixTC<Coord>
operator * ( const R4MatrixTC<Coord>& m, const Coord s )
{
    R4MatrixTC<Coord> ret;

    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = m(i, j) * s;
    return ret;
}

/// Multiply a matrix by a scaler
template<class Coord>
inline R4MatrixTC<Coord>
operator * ( const Coord s, const R4MatrixTC<Coord>& m )
{
    R4MatrixTC<Coord> ret;

    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = m(i, j) * s;
    return ret;
}

/// Multiply a matrix by a matrix
template<class Coord>
inline R4MatrixTC<Coord> 
operator * ( const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2 )
{
    R4MatrixTC<Coord> ret;

    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = m1(i, 0) * m2(0,j) + m1(i,1) * m2(1, j) + m1(i,2) * m2(2,j) + m1(i,3) * m2(3,j);
    return ret;
}

/// Multiply a matrix by a scaler
template<class Coord>
inline R4MatrixTC<Coord>
operator / ( const R4MatrixTC<Coord>& m, const Coord s )
{
    R4MatrixTC<Coord> ret;

    ASSERT( s != 0.0 );

    Coord s1 = 1.0 / s;

    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = m(i, j) * s1;
    return ret;
}

/// Multiply a 4D vector by a 4D matrix (rhs multiplication)
template<class Coord>
inline R4VectorTC<Coord>
operator * ( const R4MatrixTC<Coord>& m, const R4VectorTC<Coord>& v)
{
    return R4VectorTC<Coord>(v[0]*m(0,0) + v[1]*m(0,1) + v[2]*m(0,2) + v[3]*m(0,3),
                             v[0]*m(1,0) + v[1]*m(1,1) + v[2]*m(1,2) + v[3]*m(1,3),
                             v[0]*m(2,0) + v[1]*m(2,1) + v[2]*m(2,2) + v[3]*m(2,3),
                             v[0]*m(3,0) + v[1]*m(3,1) + v[2]*m(3,2) + v[3]*m(3,3)
                             );
}

/// Multiply a 3D vector by a 4D matrix (rhs multiplication). Makes 4th coordinate 0, doesn't divide
template<class Coord>
inline R3VectorTC<Coord>
operator * ( const R4MatrixTC<Coord>& m, const R3VectorTC<Coord>& v)
{
    R4VectorTC<Coord>  vec(v[0]*m(0,0) + v[1]*m(0,1) + v[2]*m(0,2),
                           v[0]*m(1,0) + v[1]*m(1,1) + v[2]*m(1,2),
                           v[0]*m(2,0) + v[1]*m(2,1) + v[2]*m(2,2),
                           v[0]*m(3,0) + v[1]*m(3,1) + v[2]*m(3,2)
                           );

    return R3VectorTC<Coord>( vec[0], vec[1], vec[2] );
}

/// Multiply a 4D covector by a 4D matrix (lhs multiplication).
template<class Coord>
inline R4CoVectorTC<Coord>
operator * ( const R4CoVectorTC<Coord>& v, const R4MatrixTC<Coord>& m )
{
    return R4CoVectorTC<Coord>(v[0]*m(0,0) + v[1]*m(0,1) + v[2]*m(0,2) + v[3]*m(0,3),
                               v[0]*m(1,0) + v[1]*m(1,1) + v[2]*m(1,2) + v[3]*m(1,3),
                               v[0]*m(2,0) + v[1]*m(2,1) + v[2]*m(2,2) + v[3]*m(2,3),
                               v[0]*m(3,0) + v[1]*m(3,1) + v[2]*m(3,2) + v[3]*m(3,3)
                               );
}

/// Multiply a 4D point by a 4D matrix (rhs multiplication).
template<class Coord>
inline R4PointTC<Coord>
operator * ( const R4MatrixTC<Coord>& m, const R4PointTC<Coord>& p )
{
    return R4PointTC<Coord>(p[0]*m(0,0) + p[1]*m(0,1) + p[2]*m(0,2) + p[3]*m(0,3),
                            p[0]*m(1,0) + p[1]*m(1,1) + p[2]*m(1,2) + p[3]*m(1,3),
                            p[0]*m(2,0) + p[1]*m(2,1) + p[2]*m(2,2) + p[3]*m(2,3),
                            p[0]*m(3,0) + p[1]*m(3,1) + p[2]*m(3,2) + p[3]*m(3,3)
                            );
}

/// Multiply a 3D point by a 4D matrix (rhs multiplication). Makes 4th coordinate 1, divides by w. Fails if zero divide
template<class Coord>
inline R3PointTC<Coord>
operator * ( const R4MatrixTC<Coord>& m, const R3PointTC<Coord>& p )
{
    R4PointTC<Coord>  pt(p[0]*m(0,0) + p[1]*m(0,1) + p[2]*m(0,2) + m(0,3),
                         p[0]*m(1,0) + p[1]*m(1,1) + p[2]*m(1,2) + m(1,3),
                         p[0]*m(2,0) + p[1]*m(2,1) + p[2]*m(2,2) + m(2,3),
                         p[0]*m(3,0) + p[1]*m(3,1) + p[2]*m(3,2) + m(3,3)
                        );

    ASSERT( !RNIsZero( pt[3] ) );
    return R3PointTC<Coord>( pt[0] / pt[3], pt[1] / pt[3], pt[2] / pt[3] );
}


//@}

/** \ingroup SimpleGeometryCompare */
//@{

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R4VectorTC<Coord>& v1, const R4VectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3]) ? FALSE : TRUE;
}

/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R4VectorTC<double>& v1, const R4VectorTC<double>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_d  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_d  &&
             fabs( v1[3] - v2[3] ) < RNEpsilon_d) ? TRUE : FALSE; 
}


/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R4VectorTC<float>& v1, const R4VectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_f  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_f  &&
             fabs( v1[3] - v2[3] ) < RNEpsilon_f) ? TRUE : FALSE; 
}
/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R4CoVectorTC<Coord>& v1, const R4CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3]) ? FALSE : TRUE;
}

/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R4CoVectorTC<double>& v1, const R4CoVectorTC<double>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_d  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_d  &&
             fabs( v1[3] - v2[3] ) < RNEpsilon_d) ? TRUE : FALSE; 
}


/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R4CoVectorTC<float>& v1, const R4CoVectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_f  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_f  &&
             fabs( v1[3] - v2[3] ) < RNEpsilon_f) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R4PointTC<Coord>& p1, const R4PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1[3] == p2[3]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R4PointTC<Coord>& p1, const R4PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2] && p1[3] == p2[3]) ? FALSE : TRUE;
}

/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R4PointTC<double>& p1, const R4PointTC<double>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_d &&
             fabs( p1[1] - p2[1] ) < RNEpsilon_d &&
             fabs( p1[2] - p2[2] ) < RNEpsilon_d &&
             fabs( p1[3] - p2[3] ) < RNEpsilon_d ) ? TRUE : FALSE; 
}

/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R4PointTC<float>& p1, const R4PointTC<float>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_f &&
             fabs( p1[1] - p2[1] ) < RNEpsilon_f &&
             fabs( p1[2] - p2[2] ) < RNEpsilon_f &&
             fabs( p1[3] - p2[3] ) < RNEpsilon_f ) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2)
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            if ( m1(i,j) != m2(i,j) )
                return FALSE;
    return TRUE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2)
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            if ( m1(i,j) != m2(i,j) )
                return TRUE;
    return FALSE;
}


/// Checks if elements within in_dEps of each other
inline
WINbool 
ApproxEqual ( const R4MatrixTC<double> &m1, const R4MatrixTC<double> &m2, const double in_dEps = RNEpsilon_d )
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            if ( fabs( m1(i,j) - m2(i,j) ) > in_dEps )
                return FALSE;
    return TRUE;
}


/// Checks if elements within in_fEps of each other
inline WINbool
ApproxEqual ( const R4MatrixTC<float>& m1, const R4MatrixTC<float>& m2, const float in_fEps = RNEpsilon_f)
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            if ( fabs( m1(i,j) - m2(i,j) ) > in_fEps )
                return FALSE;
    return TRUE;
}

//@}

/** \ingroup SimpleGeometryLinInterp Matrix*/
//@{


/** \brief Take a percentage of a matrix transformation

  * E.g., for a translation matrix, translate s percentage along the way. If you took s to be 0.5, and applied the matrix you got back twice, it should be the same as applying m1 once.
  * @param s Percentage (zero gives identity matrix, one gives m1). Can use negative percentages.
  * @param m1 Matrix
  * @returns A matrix that will do the percentage transformation.
  */
template<class Coord>
inline R4MatrixTC<Coord>
ScalarMult ( Coord s, const R4MatrixTC<Coord>& m1) {
   return (s * m1.Log()).Exp();
}

/** \brief Add together two matrix transformations
  * Used in Lerp, below
  * @param m1 Matrix 1
  * @param m2 Matrix 2
  * @returns m1 + m2
  */
template<class Coord>
inline R4MatrixTC<Coord>
LinearComb (const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2) {
   return (m1.Log() + m2.Log()).Exp();
}

/** \brief Take a linear combination of two matrix transformations

  * @param s Percentage (zero gives m1, one gives m2). Can use negative percentages.
  * @param m1 Matrix
  * @param m2 Matrix
  * @returns A matrix that will do the percentage transformation.
  */
template<class Coord>
inline R4MatrixTC<Coord>
Lerp( const R4MatrixTC<Coord>& m1, const R4MatrixTC<Coord>& m2, Coord s )
{
   return LinearComb( ScalarMult(1.0-s, m1), ScalarMult(s, m2) );
}

/** \brief Take a weighted combination of n matrix transformations

  * @param weights Blend values. Should sum to one and be non-zero
  * @param mats Input matrices
  * @returns A matrix that will do the blended transformation.
  */
template<class Coord>
inline R4MatrixTC<Coord>
Blend( Array< R4MatrixTC<Coord> >& mats, const Array<double>& weights )
{
    ASSERT(mats.num() == weights.num());
    
    R4MatrixTC<Coord> out;
    if (weights.num() < 1) {
        out.MakeIdentity();
        return out;
    } else out = weights[0] * mats[0].Log();
    
    for (int i = 1; i < mats.num(); i++){
        out += ( weights[i] * mats[i].Log() );
    }
    return out.Exp();
}
//@}

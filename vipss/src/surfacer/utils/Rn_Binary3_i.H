// =========================================================
//
//	Binary operators on affine3 classes
//
// =========================================================

/** \ingroup SimpleGeometryAlgebraic */
//@{

/// Add two vectors
template<class Coord>
inline R3VectorTC<Coord>
operator + ( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return R3VectorTC<Coord>( v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] );
}

/// Subtract two vectors
template<class Coord>
inline R3VectorTC<Coord>
operator - ( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return R3VectorTC<Coord>( v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R3VectorTC<Coord>
operator * ( const Coord& s, const R3VectorTC<Coord>& v )
{
    return R3VectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R3VectorTC<Coord>
operator * ( const R3VectorTC<Coord>& v, const Coord& s )
{
    return R3VectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R3VectorTC<Coord>
operator / ( const R3VectorTC<Coord>& v, const Coord& s )
{
    ASSERT( !RNIsZero( s ) );

    Coord s1 = (Coord) 1.0 / s;
    return R3VectorTC<Coord>( v[0] * s1, v[1] * s1, v[2] * s1 );
}


/// Dot product of two vectors
template<class Coord>
inline Coord
Dot( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/// Dot product of a vector and a co vector
template<class Coord>
inline Coord
Dot( const R3VectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R3CoVectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/// Cross product
template<class Coord>
inline R3VectorTC<Coord>
Cross( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return R3VectorTC<Coord>(v1[1]*v2[2] - v1[2]*v2[1],
                            v1[2]*v2[0] - v1[0]*v2[2],
                            v1[0]*v2[1] - v1[1]*v2[0]);
}

/** /brief Projection: the projection of v2 onto v1
 Returns a vector that points in the direction of v1.
 The length of the vector is the amount of v2 that points in the direction of v1
 @param v1 vector to project v2 onto. 
 @param v2 vector to project onto v1
 @returns vector */
template<class Coord>
inline R3VectorTC<Coord>
Projection( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return (Dot(v1,v2)/Dot(v1,v1)) * v1;
}

/** /brief Rejection: the rejection of v2 on v1
 Returns a vector that does not point in the direction of v1.
 The length of the vector is the amount of v2 that does not point in the direction of v1
 @param v1 vector to project v2 onto. 
 @param v2 vector to project onto v1
 @returns vector pointing in the direction of v2 and not v1 */
template<class Coord>
inline R3VectorTC<Coord>
Rejection( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return v2 - Projection(v1,v2);
}

/// Outer product of two vectors
template<class Coord>
inline R3MatrixTC<Coord>
OuterProduct( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of a vector and a covector
template<class Coord>
inline R3MatrixTC<Coord>
OuterProduct( const R3CoVectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of a covector and a vector
template<class Coord>
inline R3MatrixTC<Coord>
OuterProduct( const R3VectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of two covectors
template<class Coord>
inline R3MatrixTC<Coord>
OuterProduct( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Alternating product of two vectors
template<class Coord>
inline R3MatrixTC<Coord>
AlternatingProduct( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   ret(1,2) = -ret(2,1);

    return ret;
}

/// Alternating product of a covector and a vector
template<class Coord>
inline R3MatrixTC<Coord>
AlternatingProduct( const R3VectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   ret(1,2) = -ret(2,1);

    return ret;
}

/// Alternating product of a covector and a vector
template<class Coord>
inline R3MatrixTC<Coord>
AlternatingProduct( const R3CoVectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   ret(1,2) = -ret(2,1);

    return ret;
}

/// Alternating product of two covectors
template<class Coord>
inline R3MatrixTC<Coord>
AlternatingProduct( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    R3MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = ret(2,2) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(2,0) = v1[2]*v2[0]-v1[0]*v2[2]; ret(2,1) = v1[2]*v2[1]-v1[1]*v2[2];
    ret(0,1) = -ret(1,0);
    ret(0,2) = -ret(2,0);   ret(1,2) = -ret(2,1);

    return ret;
}

// -------------------------------------
// binary operators
// -------------------------------------

/// Subtract two points
template<class Coord>
inline R3VectorTC<Coord>
operator - ( const R3PointTC<Coord>& p1, const R3PointTC<Coord>& p2 )
{
    return R3VectorTC<Coord>( p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2] );
}

/// Add two covectors
template<class Coord>
inline R3CoVectorTC<Coord>
operator + ( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return R3CoVectorTC<Coord>( v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] );
}

/// Subtract two covectors
template<class Coord>
inline R3CoVectorTC<Coord>
operator - ( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return R3CoVectorTC<Coord>( v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] );
}

/// Scale a covector
template<class Coord>
inline R3CoVectorTC<Coord>
operator * ( const Coord& s, const R3CoVectorTC<Coord>& v )
{
    return R3CoVectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s );
}

/// Scale a covector
template<class Coord>
inline R3CoVectorTC<Coord>
operator * ( const R3CoVectorTC<Coord>& v, const Coord& s )
{
    return R3CoVectorTC<Coord>( v[0] * s, v[1] * s, v[2] * s );
}

/// Scale a covector
template<class Coord>
inline R3CoVectorTC<Coord>
operator / ( const R3CoVectorTC<Coord>& v, const Coord& s )
{
    Coord s1 = Coord(1.0f) / s;
    return R3CoVectorTC<Coord>( v[0] * s1, v[1] * s1, v[2] * s1 );
}

/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/// Cross product of two CoVectors
template<class Coord>
inline R3CoVectorTC<Coord>
Cross( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return R3CoVectorTC<Coord>(v1[1]*v2[2] - v1[2]*v2[1],
                              v1[2]*v2[0] - v1[0]*v2[2],
                              v1[0]*v2[1] - v1[1]*v2[0]);
}


/** /brief Projection: the projection of v2 onto v1
 Returns a covector that points in the direction of v1.
 The length of the covector is the amount of v2 that points in the direction of v1
 @param v1 covector to project v2 onto. 
 @param v2 covector to project onto v1
 @returns covector */
template<class Coord>
inline R3CoVectorTC<Coord>
Projection( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return (Dot(v1,v2)/Dot(v1,v1)) * v1;
}


/// Add a point and a covector
template<class Coord>
inline R3PointTC<Coord>
operator + ( const R3PointTC<Coord>& p, const R3VectorTC<Coord>& vDelta )
{
    return R3PointTC<Coord>( p[0] + vDelta[0], p[1] + vDelta[1], p[2] + vDelta[2] );
}

/// Add two covectors
template<class Coord>
inline R3PointTC<Coord>
operator + ( const R3CoVectorTC<Coord>& vDelta, const R3PointTC<Coord>& p )
{
    return R3PointTC<Coord>( p[0] + vDelta[0], p[1] + vDelta[1], p[2] + vDelta[2] );
}

/// Subtrace two covectors
template<class Coord>
inline R3PointTC<Coord>
operator - ( const R3PointTC<Coord>& p, const R3VectorTC<Coord>& vDelta )
{
    return R3PointTC<Coord>( p[0] - vDelta[0], p[1] - vDelta[1], p[2] - vDelta[2] );
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
inline R3PointTC<Coord>
Lerp( const R3PointTC<Coord>&p1, const R3PointTC<Coord>&p2, const Coord s )
{
    R3VectorTC<Coord> p1p2 = p2 - p1;
    return R3PointTC<Coord>( p1 + s*p1p2 );
}

/** \brief Barycentric combination of three points
 * It's expected that 0 <= f,g <= 1 and f + g <= 1
 * For full barycentric coords, see UTILSBary in the utils/utils module
 * If f is one, get first point. If g is one, get second point. If both zero, get third point
 * @param p First point
 * @param q Second point
 * @param r Third point
 * @param dS Percentage along (q-p), zero gives p, one gives q
 * @param dT Percentage along (q-p), zero gives p, one gives r
 * @returns Point in p,q,r triangle
*/
template<class Coord>
inline R3PointTC<Coord>
Barycentric( const R3PointTC<Coord>& p, 
             const R3PointTC<Coord>& q, 
             const R3PointTC<Coord>& r,
             const Coord dS, const Coord dT)
{
    return p+dS*(q-p)+dT*(r-p);
}


/// Add two matrices
template<class Coord>
inline R3MatrixTC<Coord>
operator + ( const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2 )
{
    R3MatrixTC<Coord> ret;

    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = m1(i, j) + m2(i, j);
    return ret;
}

/// Subtract two matrics
template<class Coord>
inline R3MatrixTC<Coord>
operator - ( const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2 )
{
    R3MatrixTC<Coord> ret;

    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = m1(i, j) - m2(i, j);
    return ret;
}

/// Multiply a matrix by a scaler
template<class Coord>
inline R3MatrixTC<Coord>
operator * ( const R3MatrixTC<Coord>& m, const Coord s )
{
    R3MatrixTC<Coord> ret;

    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = m(i, j) * s;
    return ret;
}

/// Multiply a matrix by a scaler
template<class Coord>
inline R3MatrixTC<Coord>
operator * ( const Coord s, const R3MatrixTC<Coord>& m )
{
    R3MatrixTC<Coord> ret;

    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = m(i, j) * s;
    return ret;
}

/// Multiply a matrix by a matrix
template<class Coord>
inline R3MatrixTC<Coord> 
operator * ( const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2 )
{
    R3MatrixTC<Coord> ret;

    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = m1(i, 0) * m2(0,j) + m1(i,1) * m2(1, j) + m1(i,2) * m2(2,j);
    return ret;
}

/// Multiply a matrix by a scaler
template<class Coord>
inline R3MatrixTC<Coord>
operator / ( const R3MatrixTC<Coord>& m, const Coord s )
{
    R3MatrixTC<Coord> ret;

    ASSERT( s != 0.0 );

    Coord s1 = 1.0 / s;

    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = m(i, j) * s1;
    return ret;
}

/// Multiply a 3D vector by a 3D matrix (rhs multiplication)
template<class Coord>
inline R3VectorTC<Coord>
operator * ( const R3MatrixTC<Coord>& m, const R3VectorTC<Coord>& v)
{
    return R3VectorTC<Coord>(v[0]*m(0,0) + v[1]*m(0,1) + v[2]*m(0,2),
                             v[0]*m(1,0) + v[1]*m(1,1) + v[2]*m(1,2),
                             v[0]*m(2,0) + v[1]*m(2,1) + v[2]*m(2,2)
                             );
}

/// Multiply a 2D vector by a 3D matrix (rhs multiplication)
template<class Coord>
inline R2VectorTC<Coord>
operator * ( const R3MatrixTC<Coord>& m, const R2VectorTC<Coord>& v)
{
    R3VectorTC<Coord> vec (v[0]*m(0,0) + v[1]*m(0,1) + m(0,2),
                           v[0]*m(1,0) + v[1]*m(1,1) + m(1,2),
                           v[0]*m(2,0) + v[1]*m(2,1) + m(2,2)
                           );

    return R2VectorTC<Coord>( vec[0], vec[1] );
}

/// Multiply a 3D covector by a 3D matrix (lhs multiplication)
template<class Coord>
inline R3CoVectorTC<Coord>
operator * ( const R3CoVectorTC<Coord>& v, const R3MatrixTC<Coord>& m )
{
    return R3CoVectorTC<Coord>(v[0]*m(0,0) + v[1]*m(1,0) + v[2]*m(2,0),
                               v[0]*m(0,1) + v[1]*m(1,1) + v[2]*m(2,1),
                               v[0]*m(0,2) + v[1]*m(1,2) + v[2]*m(2,2)
                               );
}

/// Multiply a 2D vector by a 3D matrix (lhs multiplication)
template<class Coord>
inline R2CoVectorTC<Coord>
operator * ( const R2CoVectorTC<Coord>& v, const R3MatrixTC<Coord>& m )
{
    R3CoVectorTC<Coord> vec(v[0]*m(0,0) + v[1]*m(1,0) + m(2,0),
                            v[0]*m(0,1) + v[1]*m(1,1) + m(2,1),
                            v[0]*m(0,2) + v[1]*m(1,2) + m(2,2)
                            );

    return R2CoVectorTC<Coord>( vec[0], vec[1] );
}

/// Multiply a 3D point by a 3D matrix (rhs multiplication)
template<class Coord>
inline R3PointTC<Coord>
operator * ( const R3MatrixTC<Coord>& m, const R3PointTC<Coord>& p )
{
    return R3PointTC<Coord>(p[0]*m(0,0) + p[1]*m(0,1) + p[2]*m(0,2),
                            p[0]*m(1,0) + p[1]*m(1,1) + p[2]*m(1,2),
                            p[0]*m(2,0) + p[1]*m(2,1) + p[2]*m(2,2)
                            );
}


/// Multiply a 2D point by a 3D matrix (rhs multiplication)
template<class Coord>
inline R2PointTC<Coord>
operator * ( const R3MatrixTC<Coord>& m, const R2PointTC<Coord>& p )
{
    R3PointTC<Coord> pt( p[0]*m(0,0) + p[1]*m(0,1) + m(0,2),
                         p[0]*m(1,0) + p[1]*m(1,1) + m(1,2),
                         p[0]*m(2,0) + p[1]*m(2,1) + m(2,2)
                         );

    ASSERT( !RNIsZero( pt[2] ) );
    return R2PointTC<Coord>(pt[0] / pt[2], pt[1] / pt[2] );
}

//@}

/** \ingroup SimpleGeometryCompare */
//@{
/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R3VectorTC<Coord>& v1, const R3VectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) ? FALSE : TRUE;
}

/// Checks if elements within in_dEps of each other
template<class Coord>
inline WINbool
ApproxEqual(const R3VectorTC<Coord> & v1, const R3VectorTC<Coord> & v2, const Coord in_dEps )
{
    if ( RNApproxEqual( v1[0], v2[0], in_dEps ) && RNApproxEqual( v1[1], v2[1], in_dEps ) && RNApproxEqual( v1[2], v2[2], in_dEps ) )
        return TRUE;

    return FALSE;
}

/// Checks if elements within in_dEps of each other
inline WINbool
ApproxEqual(const R3VectorTC<double>& v1, const R3VectorTC<double>& v2 )
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_d  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_d) ? TRUE : FALSE; 
}


/// Checks if elements within in_fEps of each other
inline WINbool
ApproxEqual(const R3VectorTC<float>& v1, const R3VectorTC<float>& v2, const float in_fEps = RNEpsilon_f )
{
    return ( fabs( v1[0] - v2[0] ) < in_fEps  &&
             fabs( v1[1] - v2[1] ) < in_fEps  &&
             fabs( v1[2] - v2[2] ) < in_fEps) ? TRUE : FALSE; 
}
/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R3CoVectorTC<Coord>& v1, const R3CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]) ? FALSE : TRUE;
}

/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R3CoVectorTC<double>& v1, const R3CoVectorTC<double>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_d  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_d) ? TRUE : FALSE; 
}


/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R3CoVectorTC<float>& v1, const R3CoVectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_f  &&
             fabs( v1[2] - v2[2] ) < RNEpsilon_f) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R3PointTC<Coord>& p1, const R3PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R3PointTC<Coord>& p1, const R3PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]) ? FALSE : TRUE;
}


/// Checks if elements within in_eps of each other
template<class Coord>
inline WINbool
ApproxEqual(const R3PointTC<Coord>& p1, const R3PointTC<Coord>& p2, const Coord in_eps)
{
    return ( fabs( p1[0] - p2[0] ) < in_eps &&
             fabs( p1[1] - p2[1] ) < in_eps &&
             fabs( p1[2] - p2[2] ) < in_eps ) ? TRUE : FALSE; 
}

/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R3PointTC<double>& p1, const R3PointTC<double>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_d &&
             fabs( p1[1] - p2[1] ) < RNEpsilon_d &&
             fabs( p1[2] - p2[2] ) < RNEpsilon_d ) ? TRUE : FALSE; 
}

/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R3PointTC<float>& p1, const R3PointTC<float>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_f &&
             fabs( p1[1] - p2[1] ) < RNEpsilon_f &&
             fabs( p1[2] - p2[2] ) < RNEpsilon_f ) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2)
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            if ( m1(i,j) != m2(i,j) )
                return FALSE;
    return TRUE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2)
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            if ( m1(i,j) != m2(i,j) )
                return TRUE;
    return FALSE;
}


/// Checks if elements within RNEpsilon_d of each other
inline
WINbool 
ApproxEqual ( const R3MatrixTC<double> &m1, const R3MatrixTC<double> &m2 )
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            if ( fabs( m1(i,j) - m2(i,j) ) > RNEpsilon_d )
                return FALSE;
    return TRUE;
}

/// Checks if elements within in_dEps of each other
inline
WINbool 
ApproxEqual ( const R3MatrixTC<double> &m1, const R3MatrixTC<double> &m2, const double in_dEps )
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            if ( fabs( m1(i,j) - m2(i,j) ) > in_dEps )
                return FALSE;
    return TRUE;
}


/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual ( const R3MatrixTC<float>& m1, const R3MatrixTC<float>& m2)
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            if ( fabs( m1(i,j) - m2(i,j) ) > RNEpsilon_f )
                return FALSE;
    return TRUE;
}
/// Checks if elements within in_fEps of each other
inline WINbool
ApproxEqual ( const R3MatrixTC<float>& m1, const R3MatrixTC<float>& m2, const float in_fEps)
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
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
inline R3MatrixTC<Coord>
ScalarMult ( Coord s, const R3MatrixTC<Coord>& m1) 
{
    WINbool bSucLog = FALSE, bSucExp = FALSE;

    const R3MatrixTC<Coord> matRet = (s * m1.Log(bSucLog)).Exp(bSucExp);

    if ( bSucLog == FALSE || bSucExp == FALSE )
        ASSERT(FALSE);

    return matRet;
}

/** \brief Add together two matrix transformations
  * Used in Lerp, below
  * @param m1 Matrix 1
  * @param m2 Matrix 2
  * @returns m1 + m2
  */
template<class Coord>
inline R3MatrixTC<Coord>
LinearComb (const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2) {
   return (m1.Log() + m2.Log()).Exp();
}

/** \brief Take a linear combination of two matrix transformations

  * @param s Percentage (zero gives m1, one gives m2). Can use negative percentages.
  * @param m1 Matrix
  * @param m2 Matrix
  * @returns A matrix that will do the percentage transformation.
  */
template<class Coord>
inline R3MatrixTC<Coord>
Lerp( const R3MatrixTC<Coord>& m1, const R3MatrixTC<Coord>& m2, const Coord s )
{
   return LinearComb( ScalarMult(1.0-s, m1), ScalarMult(s, m2) );
}

/** \brief Take a weighted combination of n matrix transformations

  * @param weights Blend values. Should sum to one and be non-zero
  * @param mats Input matrices
  * @returns A matrix that will do the blended transformation.
  */
template<class Coord>
inline R3MatrixTC<Coord>
Blend( const Array< R3MatrixTC<Coord> >& mats, const Array<Coord>& weights )
{
    ASSERT(mats.num() == weights.num());
    
    R2MatrixTC<Coord> out;
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

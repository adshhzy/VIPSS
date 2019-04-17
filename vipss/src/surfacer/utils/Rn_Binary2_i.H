// =========================================================
//
//	Binary operators on affine2 classes
//
// =========================================================

/** \ingroup SimpleGeometryAlgebraic */
//@{

/// Add two vectors
template<class Coord>
inline R2VectorTC<Coord>
operator + ( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    return R2VectorTC<Coord>( v1[0] + v2[0], v1[1] + v2[1] );
}

/// Subtract two vectors
template<class Coord>
inline R2VectorTC<Coord>
operator - ( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    return R2VectorTC<Coord>( v1[0] - v2[0], v1[1] - v2[1] );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R2VectorTC<Coord>
operator * ( const Coord& s, const R2VectorTC<Coord>& v )
{
    return R2VectorTC<Coord>( v[0] * s, v[1] * s );
}

/// Multiply a vector by a scaler
template<class Coord>
inline R2VectorTC<Coord>
operator * ( const R2VectorTC<Coord>& v, const Coord& s )
{
    return R2VectorTC<Coord>( v[0] * s, v[1] * s );
}

/// Divide a vector by a scaler
template<class Coord>
inline R2VectorTC<Coord>
operator / ( const R2VectorTC<Coord>& v, const Coord& s )
{
    ASSERT( !RNIsZero( s ) );

    return R2VectorTC<Coord>( v[0] / s, v[1] / s );
}


/// Dot product of two vectors
template<class Coord>
inline Coord
Dot( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

/// Dot product of a vector and a co vector
template<class Coord>
inline Coord
Dot( const R2VectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R2CoVectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

/// Counter-clockwise multiplication of two vectors (signed z component of cross product)
template<class Coord>
inline Coord
CCW( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    return v1[0] * v2[1] - v1[1] * v2[0];
}

/** /brief Projection: the projection of v2 onto v1
 Returns a vector that points in the direction of v1.
 The length of the vector is the amount of v2 that points in the direction of v1
 @param v1 vector to project v2 onto. 
 @param v2 vector to project onto v1
 @returns vector */
template<class Coord>
inline R2VectorTC<Coord>
Projection( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
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
inline R2VectorTC<Coord>
Rejection( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    return v2 - Projection(v1,v2);
}

/// Outer product of two vectors
template<class Coord>
inline R2MatrixTC<Coord>
OuterProduct( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of a vector and a covector
template<class Coord>
inline R2MatrixTC<Coord>
OuterProduct( const R2VectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of a covector and a vector
template<class Coord>
inline R2MatrixTC<Coord>
OuterProduct( const R2CoVectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Outer product of two covectors
template<class Coord>
inline R2MatrixTC<Coord>
OuterProduct( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            ret(i,j)=v1[i]*v2[j];
    return ret;
}

/// Alternating product of two vectors
template<class Coord>
inline R2MatrixTC<Coord>
AlternatingProduct( const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(0,1) = -ret(1,0);

    return ret;
}

/// Alternating product of a covector and a vector
template<class Coord>
inline R2MatrixTC<Coord>
AlternatingProduct( const R2VectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(0,1) = -ret(1,0);

    return ret;
}

/// Alternating product of a covector and a vector
template<class Coord>
inline R2MatrixTC<Coord>
AlternatingProduct( const R2CoVectorTC<Coord>& v1, const R2VectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(0,1) = -ret(1,0);

    return ret;
}

/// Alternating product of two covectors
template<class Coord>
inline R2MatrixTC<Coord>
AlternatingProduct( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    R2MatrixTC<Coord> ret;

    ret(0,0) = ret(1,1) = 0.0f;
    ret(1,0) = v1[1]*v2[0]-v1[0]*v2[1];
    ret(0,1) = -ret(1,0);

    return ret;
}

/// Subtract two points to get a vector
template<class Coord>
inline R2VectorTC<Coord>
operator-( const R2PointTC<Coord>& p1, const R2PointTC<Coord>& p2)
{
    return R2VectorTC<Coord>( p1[0] - p2[0], p1[1] - p2[1] );
}

/// Add two covectors
template<class Coord>
inline R2CoVectorTC<Coord>
operator + ( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return R2CoVectorTC<Coord>( v1[0] + v2[0], v1[1] + v2[1] );
}

/// Subtract two covectors
template<class Coord>
inline R2CoVectorTC<Coord>
operator - ( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return R2CoVectorTC<Coord>( v1[0] - v2[0], v1[1] - v2[1] );
}

/// Scale a covector
template<class Coord>
inline R2CoVectorTC<Coord>
operator * ( const Coord& s, const R2CoVectorTC<Coord>& v )
{
    return R2CoVectorTC<Coord>( v[0] * s, v[1] * s );
}

/// Scale a covector
template<class Coord>
inline R2CoVectorTC<Coord>
operator * ( const R2CoVectorTC<Coord>& v, const Coord& s )
{
    return R2CoVectorTC<Coord>( v[0] * s, v[1] * s );
}

/// Scale a covector
template<class Coord>
inline R2CoVectorTC<Coord>
operator / ( const R2CoVectorTC<Coord>& v, const Coord& s )
{
    return R2CoVectorTC<Coord>(v[0] / s, v[1] / s);
}


/// Dot product of two covectors
template<class Coord>
inline Coord
Dot( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[0] + v1[1] * v2[1];
}

/// Counter-clockwise multiplication of two vectors (signed z component of cross product)
template<class Coord>
inline Coord
CCW( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return v1[0] * v2[1] - v1[1] * v2[0];
}

/** /brief Projection: the projection of v2 onto v1
 Returns a covector that points in the direction of v1.
 The length of the covector is the amount of v2 that points in the direction of v1
 @param v1 covector to project v2 onto. 
 @param v2 covector to project onto v1
 @returns covector */
template<class Coord>
inline R2CoVectorTC<Coord>
Projection( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return (Dot(v1,v2)/Dot(v1,v1)) * v1;
}

// -------------------------------------
// binary operators
// -------------------------------------

/// Point plus vector
template<class Coord>
inline R2PointTC<Coord>
operator + ( const R2PointTC<Coord>& p, const R2VectorTC<Coord>& vDelta )
{
    return R2PointTC<Coord>( p[0] + vDelta[0], p[1] + vDelta[1] );
}

/// Vector plus point
template<class Coord>
inline R2PointTC<Coord>
operator + ( const R2VectorTC<Coord>& vDelta, const R2PointTC<Coord>& p )
{
    return R2PointTC<Coord>( p[0] + vDelta[0], p[1] + vDelta[1] );
}

/// Point plus a (negative) vector
template<class Coord>
inline R2PointTC<Coord>
operator - ( const R2PointTC<Coord>& p, const R2VectorTC<Coord>& vDelta )
{
    return R2PointTC<Coord>( p[0] - vDelta[0], p[1] - vDelta[1] );
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
inline R2PointTC<Coord>
Lerp( const R2PointTC<Coord>&p1, const R2PointTC<Coord>&p2, const Coord& s )
{
    R2VectorTC<Coord> p1p2 = p2 - p1;
    return R2PointTC<Coord>( p1 + s*p1p2 );
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
inline R2PointTC<Coord>
Barycentric( const R2PointTC<Coord>& p, 
             const R2PointTC<Coord>& q, 
             const R2PointTC<Coord>& r,
             const Coord dS, const Coord dT)
{
    return p+dS*(q-p)+dT*(r-p);
}

/// Multiply a matrix by a matrix
template<class Coord>
inline R2MatrixTC<Coord>
operator * ( const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2 )
{
    R2MatrixTC<Coord> ret;

    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = m1(i, 0) * m2(0,j) + m1(i,1) * m2(1, j);
    return ret;
}


/// Add two matrices
template<class Coord>
inline R2MatrixTC<Coord>
operator + ( const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2 )
{
    R2MatrixTC<Coord> ret;

    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = m1(i, j) + m2(i, j);
    return ret;
}

/// Subtract two matrics
template<class Coord>
inline R2MatrixTC<Coord>
operator - ( const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2 )
{
    R2MatrixTC<Coord> ret;

    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = m1(i, j) - m2(i, j);
    return ret;
}

/// Scale a matrix
template<class Coord>
inline R2MatrixTC<Coord>
operator * ( const R2MatrixTC<Coord>& m, const Coord s )
{
    R2MatrixTC<Coord> ret;

    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = m(i, j) * s;
    return ret;
}

/// Scale a matrix
template<class Coord>
inline R2MatrixTC<Coord>
operator * ( const Coord s, const R2MatrixTC<Coord>& m )
{
    R2MatrixTC<Coord> ret;

    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = m(i, j) * s;
    return ret;
}


/// Scale a matrix
template<class Coord>
inline R2MatrixTC<Coord>
operator / ( const R2MatrixTC<Coord>& m, const Coord s )
{
    R2MatrixTC<Coord> ret;

    ASSERT( s != 0.0 );

    Coord s1 = 1.0 / s;

    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = m(i, j) * s1;
    return ret;
}

/// Multiply a 2D vector by a 2D matrix (rhs multiplication)
template<class Coord>
inline 
R2VectorTC<Coord>
operator * ( const R2MatrixTC<Coord> &m, const R2VectorTC<Coord> &vec )
{
    return R2VectorTC<Coord>( m(0, 0)*vec[0] + m(0, 1)*vec[1],
                              m(1, 0)*vec[0] + m(1, 1)*vec[1] );
}

/// Multiply a 2D covector by a 2D matrix (lhs multiplication)
template<class Coord>
inline 
R2CoVectorTC<Coord>
operator * ( const R2CoVectorTC<Coord> &cov, const R2MatrixTC<Coord>& m )
{
    return R2CoVectorTC<Coord>(cov[0]*m(0, 0) + cov[1]*m(1, 0),
                               cov[0]*m(0, 1) + cov[1]*m(1, 1));
}

/// Multiply a 2D point by a 2D matrix (rhs multiplication)
template<class Coord>
inline 
R2PointTC<Coord>
operator * ( const R2MatrixTC<Coord> &m, const R2PointTC<Coord> &p )
{
    return R2PointTC<Coord>( m(0, 0)*p[0] + m(0, 1)*p[1],
                             m(1, 0)*p[0] + m(1, 1)*p[1] );
}
//@}

/** \ingroup SimpleGeometryCompare */
/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator==(const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2)
{
    return (v1[0] == v2[0] && v1[1] == v2[1]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator!=(const R2VectorTC<Coord>& v1, const R2VectorTC<Coord>& v2)
{
    return (v1[0] == v2[0] && v1[1] == v2[1]) ? FALSE : TRUE;
}

/// Checks if elements within in_dEps of each other
template<class Coord>
inline WINbool
ApproxEqual(const R2VectorTC<Coord> & v1, const R2VectorTC<Coord> & v2, const Coord in_dEps )
{
    if ( ApproxEqual( v1[0], v2[0], in_dEps ) && ApproxEqual( v1[1], v2[1], in_dEps ) )
        return TRUE;

    return FALSE;
}
/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R2VectorTC<double>& v1, const R2VectorTC<double>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_d) ? TRUE : FALSE; 
}

/// Checks if elements within in_dEpsilon of each other
inline WINbool
ApproxEqual(const R2VectorTC<double>& v1, const R2VectorTC<double>& v2, const double in_dEpsilon)
{
    return ( fabs( v1[0] - v2[0] ) < in_dEpsilon  &&
             fabs( v1[1] - v2[1] ) < in_dEpsilon) ? TRUE : FALSE; 
}


/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R2VectorTC<float>& v1, const R2VectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_f) ? TRUE : FALSE; 
}
/// Checks if elements within in_fEpsilon of each other
inline WINbool
ApproxEqual(const R2VectorTC<float>& v1, const R2VectorTC<float>& v2, const float in_fEpsilon)
{
    return ( fabs( v1[0] - v2[0] ) < in_fEpsilon  &&
             fabs( v1[1] - v2[1] ) < in_fEpsilon) ? TRUE : FALSE; 
}
/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R2CoVectorTC<Coord>& v1, const R2CoVectorTC<Coord>& v2 )
{
    return (v1[0] == v2[0] && v1[1] == v2[1]) ? FALSE : TRUE;
}

/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R2CoVectorTC<double>& v1, const R2CoVectorTC<double>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_d  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_d) ? TRUE : FALSE; 
}

/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual(const R2CoVectorTC<float>& v1, const R2CoVectorTC<float>& v2)
{
    return ( fabs( v1[0] - v2[0] ) < RNEpsilon_f  &&
             fabs( v1[1] - v2[1] ) < RNEpsilon_f) ? TRUE : FALSE; 
}


/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R2PointTC<Coord>& p1, const R2PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0] && p1[1] == p2[1]) ? TRUE : FALSE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R2PointTC<Coord>& p1, const R2PointTC<Coord>& p2 )
{
    return (p1[0] == p2[0] && p1[1] == p2[1]) ? FALSE : TRUE;
}
/// Checks if elements within RNEpsilon_d of each other
inline WINbool
ApproxEqual(const R2PointTC<double>& p1, const R2PointTC<double>& p2)
{
    return ( fabs( p1[0] - p2[0] ) < RNEpsilon_d &&
             fabs( p1[1] - p2[1] ) < RNEpsilon_d ) ? TRUE : FALSE; 
}
/// Checks if elements within in_dEps of each other
inline WINbool
ApproxEqual(const R2PointTC<double>& p1, const R2PointTC<double>& p2, const double in_dEps)
{
    return ( fabs( p1[0] - p2[0] ) < in_dEps &&
             fabs( p1[1] - p2[1] ) < in_dEps ) ? TRUE : FALSE; 
}
/// Checks if elements within in_fEps of each other
inline WINbool
ApproxEqual(const R2PointTC<float>& p1, const R2PointTC<float>& p2, const float in_fEps)
{
    return ( fabs( p1[0] - p2[0] ) < in_fEps &&
             fabs( p1[1] - p2[1] ) < in_fEps ) ? TRUE : FALSE; 
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator == ( const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2)
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            if ( m1(i,j) != m2(i,j) )
                return FALSE;
    return TRUE;
}

/// Compares floating point numbers; use ApproxEqual for fuzzy comparison
template<class Coord>
inline WINbool
operator != ( const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2)
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            if ( m1(i,j) != m2(i,j) )
                return TRUE;
    return FALSE;
}

/// Checks if elements within RNEpsilon_d of each other
inline
WINbool 
ApproxEqual ( const R2MatrixTC<double> &m1, const R2MatrixTC<double> &m2 )
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            if ( fabs( m1(i, j) - m2(i, j) ) > RNEpsilon_d )
                return FALSE;
    return TRUE;
}

/// Checks if elements within in_dEps of each other
inline
WINbool 
ApproxEqual ( const R2MatrixTC<double> &m1, const R2MatrixTC<double> &m2, const double in_dEps )
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            if ( fabs( m1(i, j) - m2(i, j) ) > in_dEps )
                return FALSE;
    return TRUE;
}

/// Checks if elements within RNEpsilon_f of each other
inline WINbool
ApproxEqual ( const R2MatrixTC<float>& m1, const R2MatrixTC<float>& m2)
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            if ( fabs( m1(i,j) - m2(i,j) ) > RNEpsilon_f )
                return FALSE;
    return TRUE;
}

/// Checks if elements within in_fEps of each other
inline WINbool
ApproxEqual ( const R2MatrixTC<float>& m1, const R2MatrixTC<float>& m2, const float in_fEps)
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
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
inline R2MatrixTC<Coord>
ScalarMult ( const Coord s, const R2MatrixTC<Coord>& m1) {
   return (s * m1.Log()).Exp();
}

/** \brief Add together two matrix transformations
  * Used in Lerp, below
  * @param m1 Matrix 1
  * @param m2 Matrix 2
  * @returns m1 + m2
  */
template<class Coord>
inline R2MatrixTC<Coord>
LinearComb (const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2) {
   return (m1.Log() + m2.Log()).Exp();
}

/** \brief Take a linear combination of two matrix transformations

  * @param s Percentage (zero gives m1, one gives m2). Can use negative percentages.
  * @param m1 Matrix
  * @param m2 Matrix
  * @returns A matrix that will do the percentage transformation.
  */
template<class Coord>
inline R2MatrixTC<Coord>
Lerp( const R2MatrixTC<Coord>& m1, const R2MatrixTC<Coord>& m2, const Coord s )
{
   return LinearComb( ScalarMult(1.0-s, m1), ScalarMult(s, m2) );
}

/** \brief Take a weighted combination of n matrix transformations

  * @param weights Blend values. Should sum to one and be non-zero
  * @param mats Input matrices
  * @returns A matrix that will do the blended transformation.
  */
template<class Coord>
inline R2MatrixTC<Coord>
Blend( const Array< R2MatrixTC<Coord> >& mats, const Array<Coord>& weights )
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

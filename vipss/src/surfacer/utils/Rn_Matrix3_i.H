// =========================================================
//
//	Methods for R3MatrixTC<Coord>TC template class
//
// =========================================================


template<class Coord>
inline
R3MatrixTC<Coord>::R3MatrixTC( const R3VectorTC<Coord>& vCol1, 
                               const R3VectorTC<Coord>& vCol2, 
                               const R3VectorTC<Coord>& vCol3)
{
    m_mat[0][0] = vCol1[0]; m_mat[0][1] = vCol2[0]; m_mat[0][2] = vCol3[0];
    m_mat[1][0] = vCol1[1]; m_mat[1][1] = vCol2[1]; m_mat[1][2] = vCol3[1];
    m_mat[2][0] = vCol1[2]; m_mat[2][1] = vCol2[2]; m_mat[2][2] = vCol3[2];
}

template<class Coord>
inline
R3MatrixTC<Coord>::R3MatrixTC( const R3CoVectorTC<Coord>& cvRow1, 
                               const R3CoVectorTC<Coord>& cvRow2,
                               const R3CoVectorTC<Coord>& cvRow3 )
{
    m_mat[0][0] = cvRow1[0]; m_mat[0][1] = cvRow1[1]; m_mat[0][2] = cvRow1[2];
    m_mat[1][0] = cvRow2[0]; m_mat[1][1] = cvRow2[1]; m_mat[1][2] = cvRow2[2];
    m_mat[2][0] = cvRow3[0]; m_mat[2][1] = cvRow3[1]; m_mat[2][2] = cvRow3[2];
}

template<class Coord>
inline
R3MatrixTC<Coord>::R3MatrixTC( const R2MatrixTC<Coord>& mat )
{
    for (int i = 0; i < 2; i++) {
        m_mat[i][2] = 0;
        m_mat[2][i] = 0;
        for (int j = 0; j < 2; j++) {
            (*this)(i,j) = mat(i, j);
        }
    }
    m_mat[2][2] = 1.0;
}

template<class Coord>
inline R3MatrixTC<Coord>&
R3MatrixTC<Coord>::operator *= ( const R3MatrixTC<Coord>& m )
{
    (*this) = (*this) * m;
    return (*this);
}

template<class Coord>
inline R3MatrixTC<Coord>&
R3MatrixTC<Coord>::operator += ( const R3MatrixTC<Coord>& m )
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            m_mat[i][j] += m(i, j);
    return *this;
}

template<class Coord>
inline R3MatrixTC<Coord>&
R3MatrixTC<Coord>::operator -= ( const R3MatrixTC<Coord>& m )
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            m_mat[i][j] -= m(i, j);
    return *this;
}

template<class Coord>
inline R3MatrixTC<Coord>&
R3MatrixTC<Coord>::operator *= ( const Coord s )
{
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            m_mat[i][j] *= s;
    return *this;
}

template<class Coord>
inline R3MatrixTC<Coord>&
R3MatrixTC<Coord>::operator /= ( const Coord s )
{
    Coord s1 = 1.0f / s;
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            m_mat[i][j] *= s1;
    return *this;
}

template<class Coord>
inline R3MatrixTC<Coord>
R3MatrixTC<Coord>::operator + ( )
{
    return *this;
}

template<class Coord>
inline R3MatrixTC<Coord>
R3MatrixTC<Coord>::operator - ( )
{
    R3MatrixTC<Coord> ret;
    
    for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            ret(i, j) = -m_mat[i][j];
    return ret;
}

// -------------------------------------
// Read/write/print functions
// -------------------------------------
template<class Coord>
inline void R3MatrixTC<Coord>::Write(ofstream &out) const
{
    for (int i = 0; i < 3; i++) {
        for ( int j = 0; j < 3; j++) {
            out << (*this)(i,j) << " ";
        }
        out << "\n";
    }
    out << "\n";
}

template<class Coord>
inline void R3MatrixTC<Coord>::WriteBinary(ofstream &out) const
{
    out.write( (const char *) m_mat, 9 * sizeof(Coord) );
}

template<class Coord>
inline WINbool R3MatrixTC<Coord>::Read(ifstream &in) 
{
    for (int i = 0; i < 3; i++) {
        for ( int j = 0; j < 3; j++) {
            in >> (*this)(i,j);
        }
    }

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline WINbool R3MatrixTC<Coord>::ReadBinary(ifstream &in) 
{
    in.read( (char *) &m_mat, 9 * sizeof(Coord) );

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline void R3MatrixTC<Coord>::Print( const WINbool in_bDoReturn ) const
{
    for (int i = 0; i < 3; i++) {
        for ( int j = 0; j < 3; j++) {
            TRACE("%f ", (*this)(i,j) );
        }
        if ( in_bDoReturn == TRUE )
            TRACE("\n");
        else
            TRACE(" ");
    }
}

template<class Coord>
inline void R3MatrixTC<Coord>::MakeIdentity(  )
{
    for (int i = 0; i < 3; i++) {
        for ( int j = 0; j < 3; j++) {
            if ( i == j )
                m_mat[i][j] = 1;
            else
                m_mat[i][j] = 0;
        }
    }
}


template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Identity( )
{
    R3MatrixTC<Coord> out_mat;

    out_mat.MakeIdentity();
    return out_mat;
}


/// Matrix taking the from vectors to the to vectors
template<class Coord>
inline WINbool R3MatrixTC<Coord>::MatrixVecToVec( R3MatrixTC<Coord>        & out_mat, 
                                                  const R3VectorTC<Coord>  & in_vFromX, 
                                                  const R3VectorTC<Coord>  & in_vFromY,
                                                  const R3VectorTC<Coord>  & in_vToX,   
                                                  const R3VectorTC<Coord>  & in_vToY )
{
    if ( RNIsZero( Length( in_vFromX ) ) ) return FALSE;
    if ( RNIsZero( Length( in_vFromY ) ) ) return FALSE;
    if ( RNIsZero( Length( in_vToX ) ) ) return FALSE;
    if ( RNIsZero( Length( in_vToY ) ) ) return FALSE;
    if ( RNIsZero( 1.0 - fabs( Dot( Unit( in_vFromX ), Unit( in_vFromY ) ) ) ) ) return FALSE;
    if ( RNIsZero( 1.0 - fabs( Dot( Unit( in_vToX ), Unit( in_vToY ) ) ) ) ) return FALSE;

    R3VectorTC<Coord> vFromZ = Cross( in_vFromX, in_vFromY );
    R3VectorTC<Coord> vToZ = Cross( in_vToX, in_vToY );

    R3MatrixTC<Coord> f1;
    R3MatrixTC<Coord> f2;
    
    for (int i = 0; i < 3; i++) {
        f1.m_mat[i][0] = in_vFromX[i];
        f1.m_mat[i][1] = in_vFromY[i];
        f1.m_mat[i][2] = vFromZ[i];

        f2.m_mat[i][0] = in_vToX[i];
        f2.m_mat[i][1] = in_vToY[i];
        f2.m_mat[i][2] = vToZ[i];
    }
    out_mat = f2 * f1.Transpose( );

    return TRUE;
}

/// Matrix taking the from vector to the to vector, rotation around normal unspecified
template<class Coord>
inline WINbool R3MatrixTC<Coord>::MatrixVecToVec( R3MatrixTC<Coord>        & out_mat, 
                                                  const R3VectorTC<Coord>  & in_vFrom, 
                                                  const R3VectorTC<Coord>  & in_vTo )
{
    if ( RNIsZero( Length( in_vFrom ) ) ) return FALSE;
    if ( RNIsZero( Length( in_vTo ) ) ) return FALSE;

    R3VectorTC<Coord> vFrom( 0.0, in_vFrom[2], -in_vFrom[1] );
    R3VectorTC<Coord> vTo(   0.0, in_vTo[2],   -in_vTo[1] );

    if ( fabs( in_vFrom[1] ) < 1e-6f && fabs( in_vFrom[2] ) < 1e-6f ) {
        if ( fabs( in_vTo[1] ) < 1e-6f && fabs( in_vTo[0] ) < 1e-6f ) {
            vFrom = R3VectorTC<Coord>( in_vFrom[2], -in_vFrom[0], 0.0 );
            vTo = R3VectorTC<Coord>( in_vTo[2], -in_vTo[0], 0.0 );
        } else {
            vFrom = R3VectorTC<Coord>( in_vFrom[1], -in_vFrom[0], 0.0 );
            vTo = R3VectorTC<Coord>( in_vTo[1], -in_vTo[0], 0.0 );
        }
    }
	
	if ( fabs( in_vTo[1] ) < 1e-6f && fabs( in_vTo[2] ) < 1e-6f ) {
		if ( fabs( in_vFrom[1] ) < 1e-6f && fabs( in_vFrom[0] ) < 1e-6f ) {
			vFrom = R3VectorTC<Coord>( in_vFrom[2], -in_vFrom[0], 0.0 );
			vTo = R3VectorTC<Coord>( in_vTo[2], -in_vTo[0], 0.0 );
		} else {
			vFrom = R3VectorTC<Coord>( in_vFrom[1], -in_vFrom[0], 0.0 );
			vTo = R3VectorTC<Coord>( in_vTo[1], -in_vTo[0], 0.0 );
		}
	}
             
    return MatrixVecToVec( out_mat, 
                           in_vFrom, UnitSafe( vFrom ),
                           in_vTo, UnitSafe( vTo ) );
}

/// This is a rotation of the coordinate axis by the given angle (i.e,
///  the point rotates the opposite direction)
/// So, if you give it the y axis and a point, it will rotate the point
///   counter clockwise as viewed from (0,1,0) looking at the origin
template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Rotation(const int axis, const Coord angle)
{
    ASSERT( axis >= 0 && axis <= 2 );

    R3MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity();

    Coord c = cos(angle), s = sin(angle);
    switch (axis) {
    case 0:
        out_mat(0,0)=1; out_mat(1,1)=c; out_mat(1,2)=s; out_mat(2,1)=-s; out_mat(2,2)=c;
        break;
    case 1:
        out_mat(1,1)=1; out_mat(2,2)=c; out_mat(2,0)=s; out_mat(0,2)=-s; out_mat(0,0)=c;
        break;
    case 2:
        out_mat(2,2)=1; out_mat(0,0)=c; out_mat(0,1)=s; out_mat(1,0)=-s; out_mat(1,1)=c;
        break;
    }

    return out_mat;
}

///
template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Inverse( WINbool &out_bool ) const
{
    R3MatrixTC<Coord> m;
    const R3MatrixTC<Coord> &t = *this;

      // In general case, invert it directly

    Coord out_dDeterminant = t.Determinant();

    out_bool = TRUE;
    if (fabs(out_dDeterminant) < 1e-25) {
        out_bool = FALSE;
        m = R3MatrixTC<Coord>::Identity();
        return m;
    }
    m(0,0) = (-t(1,2) * t(2,1) + t(1,1) * t(2,2)) / out_dDeterminant;
    m(0,1) = ( t(0,2) * t(2,1) - t(0,1) * t(2,2)) / out_dDeterminant;
    m(0,2) = (-t(0,2) * t(1,1) + t(0,1) * t(1,2)) / out_dDeterminant; 
    m(1,0) = ( t(1,2) * t(2,0) - t(1,0) * t(2,2)) / out_dDeterminant;
    m(1,1) = (-t(0,2) * t(2,0) + t(0,0) * t(2,2)) / out_dDeterminant;
    m(1,2) = ( t(0,2) * t(1,0) - t(0,0) * t(1,2)) / out_dDeterminant;
    m(2,0) = (-t(1,1) * t(2,0) + t(1,0) * t(2,1)) / out_dDeterminant;
    m(2,1) = ( t(0,1) * t(2,0) - t(0,0) * t(2,1)) / out_dDeterminant;
    m(2,2) = (-t(0,1) * t(1,0) + t(0,0) * t(1,1)) / out_dDeterminant;

    return m;
}

template<class Coord>
inline Coord R3MatrixTC<Coord>::Determinant() const {
   return (-m_mat[0][2] * m_mat[1][1] * m_mat[2][0] + 
            m_mat[0][1] * m_mat[1][2] * m_mat[2][0] + 
            m_mat[0][2] * m_mat[1][0] * m_mat[2][1] - 
            m_mat[0][0] * m_mat[1][2] * m_mat[2][1] - 
            m_mat[0][1] * m_mat[1][0] * m_mat[2][2] + 
            m_mat[0][0] * m_mat[1][1] * m_mat[2][2]);
}


///
template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Rotation(const R3VectorTC<Coord> &in_vec, const Coord angle)
{
    const R3VectorTC<Coord> vecAng = UnitSafe( in_vec );
    const Coord yzAng = atan2( -in_vec[2], in_vec[0] );
    const R3MatrixTC<Coord> yRot = R3MatrixTC<Coord>::Rotation(1, yzAng);
    const R3MatrixTC<Coord> yRotInv = yRot.Transpose();
    const R3VectorTC<Coord> vecRot = yRot * vecAng;
    const Coord xzAng = atan2( vecRot[1], vecRot[0] );
    const R3MatrixTC<Coord> xRot = R3MatrixTC<Coord>::Rotation(2, xzAng );
    const R3MatrixTC<Coord> xRotInv = xRot.Transpose();
    const R3VectorTC<Coord> vecRot2 = xRot * vecRot;
    
    return yRotInv * xRotInv * R3MatrixTC<Coord>::Rotation(0, angle) * xRot * yRot;    
}

template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Rotation2D( const Coord angle)
{
    R3MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity();

    Coord c = cos(angle), s = sin(angle);
    out_mat(0,0) = c;     out_mat(0,1) = -s;
    out_mat(1,0) = s;     out_mat(1,1) = c;
    return out_mat;
}

///
template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Scaling( const Coord in_x, const Coord in_y, const Coord in_z )
{
    R3MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity(  );

    out_mat(0,0) = in_x;
    out_mat(1,1) = in_y;
    out_mat(2,2) = in_z;

    return out_mat;
}

template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Scaling( const Coord in_scl )
{
    return Scaling( in_scl, in_scl, in_scl );
}

template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Translation( const R2VectorTC<Coord> & in_vec )
{
    R3MatrixTC<Coord> mat;
    mat.MakeIdentity( );

    mat(0,2) = in_vec[0];
    mat(1,2) = in_vec[1];

    return mat;
}

template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Transpose(  ) const
{
	R3MatrixTC<Coord> out_mat;

    for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			out_mat(i,j) = (*this)(j,i);

	return out_mat;
}

/// a b c 
/// d e f
/// g h i
template<class T>
R3VectorTC< C2PointTC<T> > R3MatrixTC<T>::GetEigenVec( const C2PointTC<T> &in_dV ) const
{
    R3VectorTC< C2PointTC<T> > out_vec;
    const C2PointTC<T> a( (*this)(0,0) ), b( (*this)(0,1) ), c( (*this)(0,2) );
    const C2PointTC<T> d( (*this)(1,0) ), e( (*this)(1,1) ), f( (*this)(1,2) );
    const C2PointTC<T> g( (*this)(2,0) ), h( (*this)(2,1) ), i( (*this)(2,2) );
    
    //               -(b d) + (a - l) (e - l)        c d + f (-a + l)
    //Out[5]= {{z -> ------------------------, y -> ------------------}}
    //                  -(c e) + b f + c l          -(c e) + b f + c l

    out_vec[0] = - c * e + b * f + c * in_dV;
    out_vec[1] =   c * c + f * ( in_dV - a );
    out_vec[2] = - b * d + ( a - in_dV ) * ( e - in_dV );
    if ( !RNIsZero( out_vec[0] ) || !RNIsZero( out_vec[1] ) || !RNIsZero( out_vec[2] ) ) {
        const C2PointTC<T> dRes = g * out_vec[0] + h * out_vec[1] + ( i - in_dV ) * out_vec[2];
        if ( RNIsZero( dRes, 1e-6 ) ) {
            const C2PointTC<T> dLen = Length( out_vec );
            if ( !RNIsZero( dLen, 1e-6 ) )
                out_vec /= dLen;
            return out_vec;
        }
    }

    //                    f g + d (-i + l)              -(e g) + d h + g l
    //Out[22]= {{y -> ------------------------, z -> ------------------------}}
    //                -(f h) + (e - l) (i - l)       -(f h) + (e - l) (i - l)


    out_vec[0] = -f * h + (e - in_dV) * (i - in_dV);
    out_vec[1] =  f * g + d * ( in_dV - i );
    out_vec[2] = -e * g + d * h + g * in_dV;
    if ( !RNIsZero( out_vec[0] ) || !RNIsZero( out_vec[1] ) || !RNIsZero( out_vec[2] ) ) {
        const C2PointTC<T> dRes = ( a - in_dV ) * out_vec[0] + b * out_vec[1] + c * out_vec[2];
        if ( RNIsZero( dRes, 1e-6 ) ) {
            const C2PointTC<T> dLen = Length( out_vec );
            if ( !RNIsZero( dLen, 1e-6 ) )
                out_vec /= dLen;
            return out_vec;
        }
    }

    //                c g + (a - l) (-i + l)       b g + h (-a + l)
    //Out[26]= {{y -> ----------------------, z -> ----------------}}
    //                  -(c h) + b i - b l         c h - b i + b l

    out_vec[0] = -c * h + b * i - b * in_dV;
    out_vec[1] =  c * g + ( a - in_dV ) * ( in_dV - i );
    out_vec[2] = -b * g - h * ( in_dV - a );
    if ( !RNIsZero( out_vec[0] ) || !RNIsZero( out_vec[1] ) || !RNIsZero( out_vec[2] ) ) {
        const C2PointTC<T> dRes = d * out_vec[0] + ( e - in_dV ) * out_vec[1] + f * out_vec[2];
        if ( RNIsZero( dRes, 1e-6 ) ) {
            const C2PointTC<T> dLen = Length( out_vec );
            if ( !RNIsZero( dLen, 1e-6 ) )
                out_vec /= dLen;
            return out_vec;
        }
    }

    out_vec = R3VectorTC< C2PointTC<T> >(0,0,0);
    if ( RNIsZero( a - in_dV ) && RNIsZero(d) && RNIsZero(g) )      out_vec[0] = 1.0;
    else if ( RNIsZero( e - in_dV ) && RNIsZero(b) && RNIsZero(h) ) out_vec[1] = 1.0;
    else if ( RNIsZero( i - in_dV ) && RNIsZero(c) && RNIsZero(f) ) out_vec[2] = 1.0;
    else {
        //ASSERT(FALSE);
    }
    return out_vec;
}


// a b c
// d e f
// g h i
// Characteristic polynomial is x^3 - Trace x^2 - (a_ij * a_ji - a_ii a_jj)(1/2) - det(a)
// Use quadratic forumala to find answer
//    -b +- sqrt b^2 - 4ac
//              / 2a
template<class T>
int R3MatrixTC<T>::Eigen( C2PointTC< T > &out_dV1, C2PointTC< T > &out_dV2, C2PointTC< T > &out_dV3,
                          R3VectorTC< C2PointTC< T > > & out_vec1, 
                          R3VectorTC< C2PointTC< T > > & out_vec2, 
                          R3VectorTC< C2PointTC< T > > &out_vec3 ) const
{
    Array<T> adPoly(4);
    adPoly[3] = 1.0;
    adPoly[2] = -Trace();
    adPoly[1] = 0.0;
    adPoly[0] = -Determinant();

    for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) {
            if ( i != j )
                adPoly[1] += (*this)(i,j) * (*this)(j,i) - (*this)(i,i) * (*this)(j,j);
        }
    }
    adPoly[1] = -adPoly[1] / 2.0;

    Array< C2PointTC<T> > adRoots;
    C2FindRoots( adPoly, adRoots );

    out_dV3 = out_dV2 = out_dV1 = 0.0;
    out_vec3 = out_vec2 = out_vec1 = R3VectorTC< C2PointTC<T> >(0,0,0);

    if ( adRoots.num() == 1 ) {
        out_dV1 = adRoots[0];
        out_vec1 = GetEigenVec( out_dV1 );

        return 1;
    } else if ( adRoots.num() == 2 ) {
        if ( norm( adRoots[0] ) > norm( adRoots[1] ) ) {
            out_dV1 = adRoots[0];
            out_dV2 = adRoots[1];
        } else {
            out_dV2 = adRoots[0];
            out_dV1 = adRoots[1];
        }
        out_vec1 = GetEigenVec( out_dV1 );
        out_vec2 = GetEigenVec( out_dV2 );

        return 2;
    } else {
        if ( norm( adRoots[0] ) > norm( adRoots[1] ) && norm( adRoots[0] ) > norm( adRoots[2] ) ) {
            out_dV1 = adRoots[0];
            if ( norm( adRoots[1] ) > norm( adRoots[2] ) ) {
                out_dV2 = adRoots[1];
                out_dV3 = adRoots[2];
            } else {
                out_dV3 = adRoots[1];
                out_dV2 = adRoots[2];
            }
        } else if ( norm( adRoots[1] ) > norm( adRoots[2] ) ) {
            out_dV1 = adRoots[1];
            if ( norm( adRoots[0] ) > norm( adRoots[2] ) ) {
                out_dV2 = adRoots[0];
                out_dV3 = adRoots[2];
            } else {
                out_dV3 = adRoots[0];
                out_dV2 = adRoots[2];
            }
        } else {
            out_dV1 = adRoots[2];
            if ( norm( adRoots[0] ) > norm( adRoots[1] ) ) {
                out_dV2 = adRoots[0];
                out_dV3 = adRoots[1];
            } else {
                out_dV3 = adRoots[0];
                out_dV2 = adRoots[1];
            }
        }
        out_vec1 = GetEigenVec( out_dV1 );
        out_vec2 = GetEigenVec( out_dV2 );
        out_vec3 = GetEigenVec( out_dV3 );

        return 3;
    }
    return 0;
}


// Matrix Log and Exponent

template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Exp( WINbool &out_bSuc ) const
{
   R3MatrixTC<Coord> A = (*this);           // call it A to be like Alexa's pseudocode
   R3MatrixTC<Coord> X;  X.MakeIdentity();  // the current sum
   R3MatrixTC<Coord> D;  D.MakeIdentity();  // denominator
   R3MatrixTC<Coord> N;  N.MakeIdentity();  // numerator
   Coord c = 1.0;                           // coefficienty thing

   int j = (int) max(0.0, 1.0 + floor(log(A.NormF())/log(2.0)));  // gives log2(A.NormF())
   A = A * (Coord)pow(2.0,-j);

   int q = 6;      // supposedly 6 is a good number of iterations
   for (int k = 1; k <= q; k++) {
      c = c*(q - k + 1.0)/(Coord)(k*(2*q - k + 1.0));
      X = A*X;
      N = N + c*X;
      D = D + (Coord)pow(-1.0,k)*c*X;
   }
   
   X = D.Inverse( out_bSuc ) * N;
   const int p = (int)pow(2.0,j);
   X = X.Pow(p);
   return X;
}

template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Log( WINbool &out_bSuc ) const
{
   R3MatrixTC<Coord> A = (*this);           // call it A to be like Alexa's pseudocode
   R3MatrixTC<Coord> I; I.MakeIdentity();   // identity matrix

   int k = 0;
   const int arbitrary = 30;
   while ( (A-I).NormF() > 0.000001 ) {
      A = A.Sqrt( out_bSuc );
      k++;
      if (k > arbitrary || out_bSuc == FALSE) {
          out_bSuc = FALSE;
          printf("log: repeated square roots failed to converge after %d iterations\n", arbitrary);
          break;
      }
   }
   
   A = A - I;
   R3MatrixTC<Coord> Z = A;
   R3MatrixTC<Coord> X = A;
   Coord i = 1.0;
   double eps = 0.00000001;

   while (Z.NormF() > eps) {
      Z = Z*A;
      i++;
      X = X + Z/i;
      if (i > arbitrary) { 
          out_bSuc = FALSE;
          printf("log: failed to converge after %d iterations\n", arbitrary);
          break;
      }
   }

   X = (Coord)pow(2.0,k)*X;
   return X;
}


template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Sqrt( WINbool &out_bSuc ) const
{
   const R3MatrixTC<Coord> &A = (*this);          // call it A to be like Alexa's pseudocode
   R3MatrixTC<Coord> X = (*this);
   R3MatrixTC<Coord> Y; Y.MakeIdentity();

   out_bSuc = TRUE;

   const double eps = 0.000000001;
   const int arbitrary = 20;
   WINbool bSuc = FALSE;
   int i = 0;
   while ((X*X - A).NormF() > eps) {
      R3MatrixTC<Coord> iX = X.Inverse(bSuc);
      if ( bSuc == FALSE ) out_bSuc = FALSE;
      R3MatrixTC<Coord> iY = Y.Inverse(bSuc);
      if ( bSuc == FALSE ) out_bSuc = FALSE;

      X = (X + iY)/(Coord)2.0;
      Y = (Y + iX)/(Coord)2.0;
      i++;
      if (i > arbitrary) { 
          out_bSuc = FALSE;
          printf("sqrt: failed to converge after %d iterations\n", arbitrary);
          break;
      }
   }
   //printf("sqrt: converged with i=%d\n", i);
   return X;
}

// returns the Frobenius norm
template<class Coord>
inline double R3MatrixTC<Coord>::NormF( ) const
{
   double out = 0;

   for (int i = 0; i < 3; i += 1)
        for (int j = 0; j < 3; j += 1)
            out += pow(m_mat[i][j], 2);

   return sqrt(out);
}

// returns the 1-Norm
template<class Coord>
inline double R3MatrixTC<Coord>::Norm1( ) const
{
   double sum = 0;
   double bestsum = 0;
   for (int j = 0; j < 3; j += 1) {
      for (int i = 0; i < 3; i += 1) {
         sum += abs(m_mat[i][j]);
      }
      if (sum > bestsum) bestsum = sum;
      sum = 0;
   }
   return bestsum;
}

// returns the infinity-Norm
template<class Coord>
inline double R3MatrixTC<Coord>::NormInf( ) const
{
   double sum = 0;
   double bestsum = 0;
   for (int i = 0; i < 3; i += 1) {
      for (int j = 0; j < 3; j += 1) {
         sum += abs(m_mat[i][j]);
      }
      if (sum > bestsum) bestsum = sum;
      sum = 0;
   }
   return bestsum;
}

// makes use of the inequalities given here: 
// http://planetmath.org/encyclopedia/MatrixPNorm.html
// to squeeze down the bounds on the 2-Norm of a matrix
// returns the average of the best bounds it can find
template<class Coord>
inline double R3MatrixTC<Coord>::Norm( ) const
{
   double norm1 = (*this).Norm1();
   double lowbound1 = 1.0/sqrt(2.0) * norm1;
   double highbound1 = sqrt(2.0) * norm1;

   double norminf = (*this).NormInf();
   double lowboundinf = 1.0/sqrt(2.0) * norminf; 
   double highboundinf = sqrt(2.0) * norminf;

   double normf = (*this).NormF();
   double lowboundf = normf/sqrt(2.0);
   double highboundf = normf;

   double lowbound = max(lowbound1, max(lowboundinf, lowboundf));
   double highbound = min(highbound1, min(highboundinf, highboundf));

   double err = abs(lowbound - highbound)/2.0;
   if (err > 1.0) printf("warning: norm error = %f\n", err);
   return (lowbound + highbound)/2.0;
}

// raises the matrix to the p power
template<class Coord>
inline R3MatrixTC<Coord> R3MatrixTC<Coord>::Pow(int p) const
{
   R3MatrixTC<Coord> out_mat = (*this);
   R3MatrixTC<Coord> temp = (*this);

   // don't deal with negative exponents, just return identity
   if (p <= 0) {
      out_mat.MakeIdentity();
      return out_mat;
   }

   for (int i = 0; i < p - 1; i++) {
      out_mat *= temp;
   }
   return out_mat;
}


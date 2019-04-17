// =========================================================
//
//	Methods for Point2TC template class
//
// =========================================================


// -------------------------------------
// constructors
// -------------------------------------

template<class Coord>
inline
R2MatrixTC<Coord>::R2MatrixTC( const R2VectorTC<Coord>& vCol1, const R2VectorTC<Coord>& vCol2 )
{
    mat[0][0] = vCol1[0]; mat[0][1] = vCol2[0];
    mat[1][0] = vCol1[1]; mat[1][1] = vCol2[1];
}

template<class Coord>
inline
R2MatrixTC<Coord>::R2MatrixTC( const R2CoVectorTC<Coord>& covRow1, 
                               const R2CoVectorTC<Coord>& covRow2 )
{
    mat[0][0] = covRow1[0]; mat[0][1] = covRow1[1];
    mat[1][0] = covRow2[0]; mat[1][1] = covRow2[1];
}

template<class Coord>
inline R2MatrixTC<Coord>&
R2MatrixTC<Coord>::operator *= ( const R2MatrixTC<Coord>& m )
{
    (*this) = (*this) * m;
    return *this;
}

template<class Coord>
inline R2MatrixTC<Coord>&
R2MatrixTC<Coord>::operator += ( const R2MatrixTC<Coord>& m )
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            mat[i][j] += m(i, j);
    return *this;
}

template<class Coord>
inline R2MatrixTC<Coord>&
R2MatrixTC<Coord>::operator -= ( const R2MatrixTC<Coord>& m )
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            mat[i][j] -= m(i, j);
    return *this;
}

template<class Coord>
inline R2MatrixTC<Coord>&
R2MatrixTC<Coord>::operator *= ( Coord s )
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            mat[i][j] *= s;
    return *this;
}

template<class Coord>
inline R2MatrixTC<Coord>&
R2MatrixTC<Coord>::operator /= ( Coord s )
{
    ASSERT( s != 0.0 );

    Coord s1 = 1.0 / s;
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            mat[i][j] *= s1;
    return *this;
}

template<class Coord>
inline R2MatrixTC<Coord>
R2MatrixTC<Coord>::operator + ( )
{
    return *this;
}

template<class Coord>
inline R2MatrixTC<Coord>
R2MatrixTC<Coord>::operator - ( )
{
    R2MatrixTC<Coord> ret;
    
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            ret(i, j) = -mat[i][j];
    return ret;
}


template<class Coord>
inline void
R2MatrixTC<Coord>::SetIdentity( )
{
    for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            mat[i][j] = 0;
    mat[0][0] = 1.0;
    mat[1][1] = 1.0;
}




template<class Coord>
inline R2MatrixTC<Coord>
R2MatrixTC<Coord>::Inverse( Coord &out_dDeterminant ) const
{
    R2MatrixTC<Coord> out_m;

    out_dDeterminant = this->Determinant();
    if (out_dDeterminant == 0.0) {
       out_m.SetIdentity();
       return out_m;
    }

    out_m.mat[0][0] = mat[1][1] / out_dDeterminant;
    out_m.mat[1][1] = mat[0][0] / out_dDeterminant;
    out_m.mat[0][1] = -mat[0][1] / out_dDeterminant;
    out_m.mat[1][0] = -mat[1][0] / out_dDeterminant;

    return out_m;
}

template<class Coord>
inline Coord
R2MatrixTC<Coord>::Determinant() const {
   return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

/// Angle for a Givens rotation of the matrix (the rotation component)
template<class Coord>
inline Coord R2MatrixTC<Coord>::Givens() const
{
   if( RNApproxEqual( (*this)(0,0), (*this)(1,1) ) ) {
	   if( (*this)(1,0) < 0.0 ) {
		   return -M_PI / 2.0;
	   } else {
		   return  M_PI / 2.0;
	   }
   }
   return 0.5 * atan(2.0 * (*this)(1,0) / ( (*this)(1,1)-(*this)(0,0)) );
}

template<class Coord>
inline WINbool
R2MatrixTC<Coord>::SetRotation( const Coord in_fAngle )
{
    mat[0][0] = mat[1][1] = cos( in_fAngle );
    mat[1][0] = -sin( in_fAngle );
    mat[0][1] = sin( in_fAngle );

    return TRUE;
}


template<class Coord>
inline R2MatrixTC<Coord>
R2MatrixTC<Coord>::Rotation( const Coord in_fAngle )
{
    R2MatrixTC<Coord> mat;
    mat.SetRotation( in_fAngle );

    return mat;
}

template<class Coord>
inline R2MatrixTC<Coord>
R2MatrixTC<Coord>::Scaling( const Coord in_fScaleU, const Coord in_fScaleV )
{
    R2MatrixTC<Coord> mat;
    mat.SetScale( in_fScaleU, in_fScaleV );

    return mat;
}



template<class Coord>
inline WINbool
R2MatrixTC<Coord>::SetScale( const Coord in_fScaleU, const Coord in_fScaleV )
{
    mat[0][0] = in_fScaleU;
    mat[1][1] = in_fScaleV;
    mat[1][0] = mat[0][1] = 0;

    return TRUE;
}

// -------------------------------------
// Read/write/print functions
// -------------------------------------
template<class Coord>
inline void R2MatrixTC<Coord>::Write(ofstream &out) const
{
    out << mat[0][0] << " " << mat[0][1] << "\n";
    out << mat[1][0] << " " << mat[1][1] << "\n";
}

template<class Coord>
inline WINbool R2MatrixTC<Coord>::Read(ifstream &in) 
{
    in >> mat[0][0] >> mat[0][1] >> mat[1][0] >> mat[1][1];

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline void R2MatrixTC<Coord>::Print( const WINbool in_bDoReturn ) const
{
    cout << mat[0][0] << " " << mat[0][1];
    if ( in_bDoReturn == TRUE )
        cout << "\n";
    else
        cout << " ";

    cout << mat[1][0] << " " << mat[1][1];
    if ( in_bDoReturn == TRUE )
        cout << "\n";
    else
        cout << " ";
}

template<class Coord>
inline R2MatrixTC<Coord> R2MatrixTC<Coord>::Transpose( ) const
{
	R2MatrixTC<Coord> out_mat;

    for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			out_mat(i,j) = (*this)(j,i);

	return out_mat;
}

// a b
// c d
// Characteristic polynomial is x^2 - (a + d) x + (ad - bc)
// Use quadratic forumala to find answer
//    -b +- sqrt b^2 - 4ac
//              / 2a
template<class T>
R2VectorTC< C2PointTC<T> > R2MatrixTC< T >::GetEigenVec( const C2PointTC<T> &in_dV ) const
{
    R2VectorTC< C2PointTC<T> > out_vec;

    out_vec[0] = (*this)(0,1); // a x + b y = l x, x = 1, y = (l-a) / b
    out_vec[1] = ( in_dV - (*this)(0,0) );

    if ( !RNIsZero( out_vec[0] ) || !RNIsZero( out_vec[1] ) ) {
        const C2PointTC<T> dRes = out_vec[0] * (*this)(1,0) + out_vec[1] * ( (C2PointTC<T>) (*this)(1,1) - in_dV );
        if ( RNIsZero( dRes, 1e-6 ) ) {
            const C2PointTC<T> dLen = Length( out_vec );
            if ( !RNIsZero( dLen, 1e-6 ) )
                out_vec /= dLen;
            return out_vec;
        }
    }

    out_vec[1] = (*this)(1,0); // c x + d y = l y, y = 1, x = (l-d) / c
    out_vec[0] = ( in_dV - (*this)(1,1) );

    if ( !RNIsZero( out_vec[0] ) || !RNIsZero( out_vec[1] ) ) {
        const C2PointTC<T> dRes = ( (C2PointTC<T>) (*this)(0,0) - in_dV ) * out_vec[0] + out_vec[1] * (*this)(0,1);
        if ( RNIsZero( dRes, 1e-6 ) ) {
            const C2PointTC<T> dLen = Length( out_vec );
            if ( !RNIsZero( dLen, 1e-12 ) )
                out_vec /= dLen;
            return out_vec;
        }
    }


    out_vec[0] = 1.0;
    out_vec[1] = 0.0;
    return out_vec;
}

template<class T>
int R2MatrixTC<T>::Eigen( C2PointTC<T> &out_dV1, 
                              C2PointTC<T> &out_dV2, 
                              R2VectorTC< C2PointTC<T> > & out_vec1, 
                              R2VectorTC< C2PointTC<T> > & out_vec2 ) const
{
    Array<T> adPoly(3);
    adPoly[2] = 1.0;
    adPoly[1] = -Trace();
    adPoly[0] = Determinant();

    Array< C2PointTC<T> > adRoots;
    C2FindRoots( adPoly, adRoots );

    out_dV2 = out_dV1 = 0.0;
    out_vec2 = out_vec1 = R2VectorTC< C2PointTC<T> >(0,0);

    if ( adRoots.num() == 1 ) {
        out_dV1 = adRoots[0];
        out_vec1 = GetEigenVec( out_dV1 );

        return 1;
    } else if ( adRoots.num() == 2 ) {
        if ( norm( adRoots[0] ) > norm( adRoots[1] ) ) {
            out_dV1 = adRoots[0];
            out_dV2 = adRoots[1];
        } else if ( norm( adRoots[0] ) < norm( adRoots[1] ) ) {
            out_dV2 = adRoots[0];
            out_dV1 = adRoots[1];
        } else {
            out_dV1 = adRoots[0];
            out_dV2 = adRoots[1];

            if ( ApproxEqual( adRoots[0], adRoots[1] ) ) {
                out_vec1 = GetEigenVec( out_dV1 );
                out_vec2[0] = out_vec1[1];
                out_vec2[1] = out_vec1[0];

                return 2;
            }
        }

        out_vec1 = GetEigenVec( out_dV1 );
        out_vec2 = GetEigenVec( out_dV2 );

        return 2;
    }
    return 0;
}

template<class Coord>
int R2MatrixTC<Coord>::EigenReal( Coord &out_dE1, Coord &out_dE2, R2VectorTC< Coord > & out_vec1, R2VectorTC< Coord > & out_vec2, const Coord in_dEps ) const
{
    out_dE1 = 1.0;
    out_dE2 = 0.0;
    out_vec1 = R2VectorTC<Coord>(1,0);        
    out_vec2 = R2VectorTC<Coord>(0,1);

    // polynomial to find roots for
    const Coord a = 1.0;
    const Coord b = -Trace();
    const Coord c = Determinant();

    double dRoot = b*b - 4.0*a*c;
	if( dRoot < -in_dEps ) {
		TRACE( "Complex eigenvalues: A=%lf, B=%lf, C=%lf, R=%lf\n",
				a, b, c, dRoot );
		return 0;
	}

    // Force the square root to be zero
	if( dRoot < in_dEps ) dRoot = 0.0;

	const Coord v1 = ( -b + sqrt( dRoot ) ) / ( 2.0 * a );
	const Coord v2 = ( -b - sqrt( dRoot ) ) / ( 2.0 * a );

	if( v1*v1 > v2*v2 ) {
		out_dE1 = v1;
		out_dE2 = v2;
	} else {
		out_dE1 = v2;
		out_dE2 = v1;
	}

    const Coord x1 = (*this)(0,0) - out_dE1;
    const Coord x2 = (*this)(1,0);
    const Coord x3 = (*this)(0,1);
    const Coord x4 = (*this)(1,1) - out_dE1;
    const Coord x5 = (*this)(0,0) - out_dE2;
    const Coord x6 = (*this)(1,1) - out_dE2;

	// Special cases
	if ( RNIsZero( x1 ) && RNIsZero( x2 ) && RNIsZero( x3 ) && RNIsZero( x4 ) && RNIsZero( x5 ) && RNIsZero( x6 ) ) {
		return 1;
	}

	if( RNApproxEqual( out_dE1, out_dE2 ) ) {
		return 2;
	}
	
	// First eigenvector
	if( !RNIsZero( x1 ) || !RNIsZero( x2 ) ) {
		if( x1*x1 > x2*x2 ) {
			out_vec1[0] = -x2/x1;
			out_vec1[1] = 1.0;
		} else {
			out_vec1[0] = 1.0;
			out_vec1[1] = -x1/x2;
		}
	} else {  // x3 or x4 must be non-zero
		if( x3*x3 > x4*x4 ) {
			out_vec1[0] = -x4/x3;
			out_vec1[1] = 1.0;
		} else {
			out_vec1[0] = 1.0;
			out_vec1[1] = -x3/x4;
		}
	}

	// Second eigenvector
	if(  !RNIsZero( x5 ) || !RNIsZero( x2 )  ) {
		if( x5*x5 > x2*x2 ) {
			out_vec2[0] = -x2/x5;
			out_vec2[1] = 1.0;
		} else {
			out_vec2[0] = 1.0;
			out_vec2[1] = -x5/x2;
		}
	} else {  // x3 or x6 must be non-zero
		if( x3*x3 > x6*x6 ) {
			out_vec2[0] = -x6/x3;
			out_vec2[1] = 1.0;
		} else {
			out_vec2[0] = 1.0;
			out_vec2[1] = -x3/x6;
		}
	}
	return 2;
}


template<class Coord>
WINbool R2MatrixTC<Coord>::SVD( R2MatrixTC<Coord> &out_matU, R2MatrixTC<Coord> &out_matDiag, R2MatrixTC<Coord> &out_matV ) const
{
    R2VectorTC< C2PointTC<Coord> > vecU1, vecU2, vecV1, vecV2;
    C2PointTC<Coord> d1, d2;

    const R2MatrixTC<Coord> matAAT = (*this) * Transpose();
    const R2MatrixTC<Coord> matATA = Transpose() * (*this);

    const int iNU = matAAT.Eigen( d1, d2, vecU1, vecU2 );
    if ( iNU != 2 ) return FALSE;
    const int iNV = matATA.Eigen( d1, d2, vecV1, vecV2 );
    if ( iNV != 2 ) return FALSE;

    if ( !RNIsZero( d1.imag(), 1e-6 ) ) return FALSE;
    if ( !RNIsZero( d2.imag(), 1e-6 ) ) return FALSE;
    if (  d1.real() < 0.0 ) return FALSE;
    if (  d2.real() < 0.0 ) return FALSE;
    out_matDiag = R2MatrixTC<Coord>( R2VectorTC<Coord>( sqrt( d1.real() ), 0.0 ),
                                     R2VectorTC<Coord>( 0.0, sqrt( d2.real()) ) );

    /*
    const Coord dDot = Dot( R2Vec( vecU1[0].real(), vecU1[1].real() ), R2Vec( vecV1[0].real(), vecV1[1].real() ) );
    if ( dDot < 0.0 )
        vecU1 = vecU1 * C2PointTC<Coord>( -1.0 );

    const Coord dCrossU = vecU1[0].real() * vecU2[1].real() - vecU1[1].real() * vecU2[0].real();
    const Coord dCrossV = vecV1[0].real() * vecV2[1].real() - vecV1[1].real() * vecV2[0].real();
    if ( dCrossU < 0.0 )
        vecU2 = vecU2 * C2PointTC<Coord>( -1.0 );
    if ( dCrossV < 0.0 )
        vecV2 = vecV2 * C2PointTC<Coord>( -1.0 );
    */

    for ( int i = 0; i < 4; i++ ) {
        const int i1 = (i / 2 == 0) ? 1 : -1;
        const int i2 = (i % 2 == 0) ? 1 : -1;

        for ( int j = 0; j < 4; j++ ) {
            const int j1 = (j / 2 == 0) ? 1 : -1;
            const int j2 = (j % 2 == 0) ? 1 : -1;
            out_matU = R2MatrixTC<Coord>( R2VectorTC<Coord>( i1 * vecU1[0].real(), i1 * vecU1[1].real() ),
                                          R2VectorTC<Coord>( i2 * vecU2[0].real(), i2 * vecU2[1].real() ) );
            out_matV = R2MatrixTC<Coord>( R2VectorTC<Coord>( j1 * vecV1[0].real(), j1 * vecV1[1].real() ),
                                          R2VectorTC<Coord>( j2 * vecV2[0].real(), j2 * vecV2[1].real() ) );

            if ( ApproxEqual( out_matU * out_matDiag * out_matV.Transpose(), (*this), 1e-6 ) )
                //TRACE("foo %d %d  %d %d\n", i1, i2, j1, j2);
                return TRUE;
        }
    }
    return TRUE;
}

template<class Coord>
void R2MatrixTC<Coord>::QRDecomposition( R2MatrixTC<Coord> &out_matQ, R2MatrixTC<Coord> &out_matR ) const
{
    const R2VectorTC<Coord> vecU1( (*this)(0,0), (*this)(1,0) );
    const R2VectorTC<Coord> vecA2( (*this)(0,1), (*this)(1,1) );
    const R2VectorTC<Coord> vecE1 = UnitSafe( vecU1 );
    const R2VectorTC<Coord> vecU2 = vecA2 - Projection( vecE1, vecA2 );
    const R2VectorTC<Coord> vecE2 = UnitSafe( vecU2 );
    
    out_matQ(0,0) = vecE1[0];
    out_matQ(1,0) = vecE1[1];
    out_matQ(0,1) = vecE2[0];
    out_matQ(1,1) = vecE2[1];
    
    out_matR = out_matQ.Transpose() * (*this);
}


///
template<class Coord>
WINbool R2MatrixTC<Coord>::Takagi( R2MatrixTC<Coord> &out_matU, R2MatrixTC<Coord> &out_matDiag ) const
{
    const Coord dA = (*this)(0,0);
    const Coord dB = (*this)(0,1);
    const Coord dC = (*this)(1,1);

    const Coord dSqrt = sqrt( pow(dA, 2) + 4.0 * pow( dB, 2 ) - 2.0 * dA * dC + pow( dC, 2.0 ) );
    const Coord d1 = ( dA + dC - dSqrt ) / 2.0;
    const Coord d2 = ( dA + dC + dSqrt ) / 2.0;
    const Coord dSqrt2 = sqrt( -4.0 * pow(dB,2) + pow( (d1 - d2), 2 ) );
    const Coord dTerm = sqrt( -( dSqrt2 + d2 - d1 ) / (d1 - d2) );

    const Coord dT1 = ( (d1 + dSqrt2 - d2) * dTerm ) / ( 2.0 * sqrt(2.0) * dB );
    const Coord dT2 = dTerm / sqrt(2.0);

    out_matDiag = R2MatrixTC<Coord>( R2VectorTC<Coord>( d1, 0.0 ),
                                     R2VectorTC<Coord>( 0.0, d2 ) );

    for ( int i = 0; i < 4; i++ ) {
        const int i1 = (i / 2 == 0) ? 1 : -1;
        const int i2 = (i % 2 == 0) ? 1 : -1;
        for ( int j = 0; j < 4; j++ ) {
            const int j1 = (j / 2 == 0) ? 1 : -1;
            const int j2 = (j % 2 == 0) ? 1 : -1;

            out_matU = R2MatrixTC<Coord>( R2VectorTC<Coord>( i1 * dT1, i2 * dT2 ), R2VectorTC<Coord>( j1 * dT2, j2 * dT1 ) );

            if ( ApproxEqual( out_matU * out_matDiag * out_matU.Transpose(), *this, 1e-6 ) )
                return TRUE;
        }
    }

    for ( FORINT i = 0; i < 4; i++ ) {
        const int i1 = (i / 2 == 0) ? 1 : -1;
        const int i2 = (i % 2 == 0) ? 1 : -1;
        for ( int j = 0; j < 4; j++ ) {
            const int j1 = (j / 2 == 0) ? 1 : -1;
            const int j2 = (j % 2 == 0) ? 1 : -1;

            out_matU = R2MatrixTC<Coord>( R2VectorTC<Coord>( i1 * dT2, i2 * dT1 ), R2VectorTC<Coord>( j1 * dT1, j2 * dT2 ) );

            if ( ApproxEqual( out_matU * out_matDiag * out_matU.Transpose(), *this, 1e-6 ) )
                return TRUE;
        }
    }
    return TRUE;
}



// Matrix Log and Exponent

template<class Coord>
inline R2MatrixTC<Coord> R2MatrixTC<Coord>::Exp( ) const
{
   R2MatrixTC<Coord> A = (*this);          // call it A to be like Alexa's pseudocode
   R2MatrixTC<Coord> X;  X.SetIdentity();  // the current sum
   R2MatrixTC<Coord> D;  D.SetIdentity();  // denominator
   R2MatrixTC<Coord> N;  N.SetIdentity();  // numerator
   Coord c = 1.0;                          // coefficienty thing

   int j = (int) max(0.0, 1.0 + floor(log(A.Norm())/log(2.0)));  // gives log2(A.Norm())
   A = A * (Coord)pow(2.0,-j);

   int q = 6;      // supposedly 6 is a good number of iterations
   for (int k = 1; k <= q; k++) {
      c = c*(q - k + 1.0)/(Coord)(k*(2*q - k + 1.0));
      X = A*X;
      N = N + c*X;
      D = D + (Coord)pow(-1.0,k)*c*X;
   }
   
   Coord det;
   X = D.Inverse(det) * N;
   int p = (int)pow(2.0,j);
   X = X.Pow(p);
   return X;
}

template<class Coord>
inline R2MatrixTC<Coord> R2MatrixTC<Coord>::Log( ) const
{
   R2MatrixTC<Coord> A = (*this);          // call it A to be like Alexa's pseudocode
   R2MatrixTC<Coord> I; I.SetIdentity();   // identity matrix

   int k = 0;
   int arbitrary = 30;
   while ( (A-I).NormF() > 0.000001 ) {
      A = A.Sqrt();
      k++;
      if (k > arbitrary) {
         printf("log: repeated square roots failed to converge after %d iterations\n", arbitrary);
         break;
      }
   }
   
   A = A - I;
   R2MatrixTC<Coord> Z = A;
   R2MatrixTC<Coord> X = A;
   Coord i = 1.0;
   double eps = 0.00000001;

   while (Z.NormF() > eps) {
      Z = Z*A;
      i++;
      X = X + Z/i;
      if (i > arbitrary) { 
         printf("log: failed to converge after %d iterations\n", arbitrary);
         break;
      }
   }

   X = (Coord)pow(2.0,k)*X;
   return X;
}

template<class Coord>
inline R2MatrixTC<Coord> R2MatrixTC<Coord>::Sqrt( ) const
{
   R2MatrixTC<Coord> A = (*this);          // call it A to be like Alexa's pseudocode
   R2MatrixTC<Coord> X = (*this);
   R2MatrixTC<Coord> Y; Y.SetIdentity();

   double eps = 0.000000001;
   Coord det;
   int i = 0;
   int arbitrary = 20;
   while ((X*X - A).Norm() > eps) {
      R2MatrixTC<Coord> iX = X.Inverse(det);
      R2MatrixTC<Coord> iY = Y.Inverse(det);
      X = (X + iY)/(Coord)2.0;
      Y = (Y + iX)/(Coord)2.0;
      i++;
      if (i > arbitrary) { 
         printf("sqrt: failed to converge after %d iterations\n", arbitrary);
         break;
      }
   }
   //printf("sqrt: converged with i=%d\n", i);
   return X;
}

// returns the Frobenius norm
template<class Coord>
inline Coord R2MatrixTC<Coord>::NormF( ) const
{
   Coord out = 0;

   for (int i = 0; i < 2; i += 1)
        for (int j = 0; j < 2; j += 1)
            out += pow(mat[i][j], 2);

   return sqrt(out);
}

// returns the 1-Norm
template<class Coord>
inline Coord R2MatrixTC<Coord>::Norm1( ) const
{
   Coord sum = 0;
   Coord bestsum = 0;
   for (int j = 0; j < 2; j += 1) {
      for (int i = 0; i < 2; i += 1) {
         sum += abs(mat[i][j]);
      }
      if (sum > bestsum) bestsum = sum;
      sum = 0;
   }
   return bestsum;
}

// returns the infinity-Norm
template<class Coord>
inline Coord R2MatrixTC<Coord>::NormInf( ) const
{
   Coord sum = 0;
   Coord bestsum = 0;
   for (int i = 0; i < 2; i += 1) {
      for (int j = 0; j < 2; j += 1) {
         sum += abs(mat[i][j]);
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
inline Coord R2MatrixTC<Coord>::Norm( ) const
{
   const Coord norm1 = (*this).Norm1();
   const Coord lowbound1 = 1.0/sqrt(2.0) * norm1;
   const Coord highbound1 = sqrt(2.0) * norm1;

   const Coord norminf = (*this).NormInf();
   const Coord lowboundinf = 1.0/sqrt(2.0) * norminf; 
   const Coord highboundinf = sqrt(2.0) * norminf;

   const Coord normf = (*this).NormF();
   const Coord lowboundf = normf/sqrt(2.0);
   const Coord highboundf = normf;

   //printf("Matrix 2-norm is between %f and %f.\n", lowbound1, highbound1);
   //printf("Matrix 2-norm is between %f and %f.\n", lowboundinf, highboundinf);
   //printf("Matrix 2-norm is between %f and %f.\n", lowboundf, highboundf);

   const double lowbound = max(lowbound1, max(lowboundinf, lowboundf));
   const double highbound = min(highbound1, min(highboundinf, highboundf));

   //printf("Matrix 2-norm returning %f, error is +/-%f\n", (lowbound + highbound)/2.0, abs(lowbound - highbound)/2.0);
   const double err = abs(lowbound - highbound)/2.0;
   if (err > 1.0) printf("warning: norm error = %f\n", err);
   return (lowbound + highbound)/2.0;
}

// raises the matrix to the p power
template<class Coord>
inline R2MatrixTC<Coord> R2MatrixTC<Coord>::Pow(int p) const
{
   R2MatrixTC<Coord> out_mat = (*this);
   R2MatrixTC<Coord> temp = (*this);

   // don't deal with negative exponents, just return identity
   if (p <= 0) {
      out_mat.SetIdentity();
      return out_mat;
   }

   for (int i = 0; i < p - 1; i++) {
      out_mat *= temp;
   }
   return out_mat;
}


// =========================================================
//
//	Methods for R4MatrixTC<Coord>TC template class
//
// =========================================================


template<class Coord>
inline
R4MatrixTC<Coord>::R4MatrixTC( const R4VectorTC<Coord>& vRow1, 
                               const R4VectorTC<Coord>& vRow2, 
                               const R4VectorTC<Coord>& vRow3,
                               const R4VectorTC<Coord>& vRow4)
{
    // it says "Row" but uses the input vectors as columns! BEWARE
    m_mat[0][0] = vRow1[0]; m_mat[1][0] = vRow2[0]; m_mat[2][0] = vRow3[0]; m_mat[3][0] = vRow4[0];
    m_mat[0][1] = vRow1[1]; m_mat[1][1] = vRow2[1]; m_mat[2][1] = vRow3[1]; m_mat[3][1] = vRow4[1];
    m_mat[0][2] = vRow1[2]; m_mat[1][2] = vRow2[2]; m_mat[2][2] = vRow3[2]; m_mat[3][2] = vRow4[2];
    m_mat[0][3] = vRow1[3]; m_mat[1][3] = vRow2[3]; m_mat[2][3] = vRow3[3]; m_mat[3][3] = vRow4[3];
}

template<class Coord>
inline
R4MatrixTC<Coord>::R4MatrixTC( const R4CoVectorTC<Coord>& cvCol1, 
                               const R4CoVectorTC<Coord>& cvCol2,
                               const R4CoVectorTC<Coord>& cvCol3,
                               const R4CoVectorTC<Coord>& cvCol4 )
{
    // yep, this one is backwards also. BEWARE BEWARE
    m_mat[0][0] = cvCol1[0]; m_mat[1][0] = cvCol1[1]; m_mat[2][0] = cvCol1[2]; m_mat[3][0] = cvCol1[3];
    m_mat[0][1] = cvCol2[0]; m_mat[1][1] = cvCol2[1]; m_mat[2][1] = cvCol2[2]; m_mat[3][1] = cvCol2[3];
    m_mat[0][2] = cvCol3[0]; m_mat[1][2] = cvCol3[1]; m_mat[2][2] = cvCol3[2]; m_mat[3][2] = cvCol3[3];
    m_mat[0][3] = cvCol4[0]; m_mat[1][3] = cvCol4[1]; m_mat[2][3] = cvCol4[2]; m_mat[3][3] = cvCol4[3];
}

template<class Coord>
inline
R4MatrixTC<Coord>::R4MatrixTC( const R3MatrixTC<Coord>& mat )
{
    for (int i = 0; i < 3; i++) {
        m_mat[i][3] = 0;
        m_mat[3][i] = 0;
        for (int j = 0; j < 3; j++) {
            (*this)(i,j) = mat(i, j);
        }
    }
    m_mat[3][3] = 1.0;
}

template<class Coord>
inline R4MatrixTC<Coord>&
R4MatrixTC<Coord>::operator *= ( const R4MatrixTC<Coord>& m )
{
    (*this) = (*this) * m;
    return (*this);
}

template<class Coord>
inline R4MatrixTC<Coord>&
R4MatrixTC<Coord>::operator += ( const R4MatrixTC<Coord>& m )
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            (*this)(i,j) += m(i, j);
    return *this;
}

template<class Coord>
inline R4MatrixTC<Coord>&
R4MatrixTC<Coord>::operator -= ( const R4MatrixTC<Coord>& m )
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            (*this)(i,j) -= m(i, j);
    return *this;
}

template<class Coord>
inline R4MatrixTC<Coord>&
R4MatrixTC<Coord>::operator *= ( const Coord s )
{
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            (*this)(i,j) *= s;
    return *this;
}

template<class Coord>
inline R4MatrixTC<Coord>&
R4MatrixTC<Coord>::operator /= ( const Coord s )
{
    Coord s1 = 1.0f / s;
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            (*this)(i,j) *= s1;
    return *this;
}

template<class Coord>
inline R4MatrixTC<Coord>
R4MatrixTC<Coord>::operator + ( )
{
    return *this;
}

template<class Coord>
inline R4MatrixTC<Coord>
R4MatrixTC<Coord>::operator - ( )
{
    R4MatrixTC<Coord> ret;
    
    for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            ret(i, j) = -(*this)(i,j);
    return ret;
}

///
template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Inverse( WINbool &out_bSuccess ) const
{
    R4MatrixTC<Coord> mo;
    const R4MatrixTC<Coord> &mi = *this;
    
    Coord t[4][8], save;
    int i,j;
    // copy mi into t
    // set t's last four rows to the identity matrix
    for (i=0;i<4;i++)
        for (j=0;j<4;j++) {
            t[i][j]=mi(i, j);
            t[i][j+4]=0;
        }
    for (i=0;i<4;i++)
        t[i][i+4]=1;
    
    mo.MakeIdentity();
    for (i=0;i<4;i++) {
        if (i<3) {	// swap row with largest front coefficient
            Coord a=fabs(t[i][i]),ab;
            int m=i;
            for (int l=i+1;l<4;l++) {
                if ((ab=fabs(t[l][i]))>a) {
                    a=ab,m=l;
                }
            }
            if (m!=i) {
                for (j=0;j<8;j++) {
                    save = t[i][j];
                    t[i][j] = t[m][j];
                    t[m][j] = save;
                }
            }
        }
        if (!t[i][i]) {
            out_bSuccess = FALSE ;
            return mo;
        }
        for (j=0;j<4;j++) {
			if (i==j) {
                   Coord a=1/t[i][i];
                   for (int k=0;k<8;k++)
                       t[j][k]*=a;
			} else {
                   Coord a=-t[j][i]/t[i][i];
                   for (int k=0;k<8;k++)
                       t[j][k]+=a*t[i][k];
			}
        }
    }
    for (i=0;i<4;i++)
        for (j=0;j<4;j++)
            mo(i,j)=t[i][j+4];
    
    out_bSuccess = TRUE ;
    return mo;
}

// -------------------------------------
// Read/write/print functions
// -------------------------------------
template<class Coord>
inline void R4MatrixTC<Coord>::Write(ofstream &out) const
{
    for (int i = 0; i < 4; i++) {
        for ( int j = 0; j < 4; j++) {
            out << (*this)(i,j) << " ";
        }
        out << "\n";
    }
    out << "\n";
}

template<class Coord>
inline WINbool R4MatrixTC<Coord>::Read(ifstream &in) 
{
    for (int i = 0; i < 4; i++) {
        for ( int j = 0; j < 4; j++) {
            in >> (*this)(i,j);
        }
    }

    return in.good() ? TRUE : FALSE;
}

template<class Coord>
inline void R4MatrixTC<Coord>::Print( const WINbool in_bDoReturn ) const
{
    for (int i = 0; i < 4; i++) {
        for ( int j = 0; j < 4; j++) {
            cout << (*this)(i,j) << " ";
        }
        if ( in_bDoReturn == TRUE )
            cout << "\n";
        else
            cout << " ";
    }
}

template<class Coord>
inline void R4MatrixTC<Coord>::PrintMatlab( ) const
{
    cout << "[" ;
    for (int i = 0; i < 4; i++) {
        for ( int j = 0; j < 4; j++) {
            if (j < 3) cout << (*this)(i,j) << ", ";
            else cout << (*this)(i,j) << "\n";
        }
    }
    cout << "];";
}


template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Identity( )
{
    R4MatrixTC<Coord> out_mat;

    out_mat.MakeIdentity();
    return out_mat;
}

template<class Coord>
inline void R4MatrixTC<Coord>::MakeIdentity( )
{
    for ( int i = 0; i < 4; i++ )
		for ( int j = 0; j < 4; j++ )
			m_mat[i][j] = (i == j) ? 1.0 : 0.0;
}

/// 
template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Rotation(const int axis, const Coord angle)
{
    ASSERT( axis >= 0 && axis <= 2 );

    R4MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity();

    double c = cos(angle), s = sin(angle);
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
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Rotation(const R3VectorTC<Coord> &in_vec, const Coord angle)
{
    const R3VectorTC<Coord> vecAng = UnitSafe( in_vec );
    const Coord yzAng = atan2( -in_vec[2], in_vec[0] );
    const R4MatrixTC<Coord> yRot = R4MatrixTC<Coord>::Rotation(1, yzAng);
    const R4MatrixTC<Coord> yRotInv = yRot.Transpose();
    const R3VectorTC<Coord> vecRot = yRot * vecAng;
    const Coord xzAng = atan2( vecRot[1], vecRot[0] );
    const R4MatrixTC<Coord> xRot = R4MatrixTC<Coord>::Rotation(2, xzAng );
    const R4MatrixTC<Coord> xRotInv = xRot.Transpose();
    const R3VectorTC<Coord> vecRot2 = xRot * vecRot;
    
    return yRotInv * xRotInv * R4MatrixTC<Coord>::Rotation(0, angle) * xRot * yRot;    
}

///
template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Scaling( const Coord in_x, const Coord in_y, const Coord in_z, const Coord in_w )
{
    R4MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity(  );

    out_mat(0,0) = in_x;
    out_mat(1,1) = in_y;
    out_mat(2,2) = in_z;
    out_mat(3,3) = in_w;

    return out_mat;
}

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Scaling( const Coord in_scl )
{
    return Scaling( in_scl, in_scl, in_scl, 1.0 );
}

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Translation( const R3VectorTC<Coord> &in_vec )
{
    R4MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity( );

    out_mat(0,3) = in_vec[0];
    out_mat(1,3) = in_vec[1];
    out_mat(2,3) = in_vec[2];

    return out_mat;
}


template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Translation( const R3PointTC<Coord> &in_pt )
{
    R4MatrixTC<Coord> out_mat;
    out_mat.MakeIdentity( );

    out_mat(0,3) = in_pt[0];
    out_mat(1,3) = in_pt[1];
    out_mat(2,3) = in_pt[2];

    return out_mat;
}

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Transpose(  ) const
{
	R4MatrixTC<Coord> out_mat;

    for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			out_mat(i,j) = (*this)(j,i);

	return out_mat;
}

template<class Coord>
inline Coord R4MatrixTC<Coord>::Determinant() const {
  R3MatrixTC<Coord> mat00, mat01, mat02, mat03;

  {
    R3VectorTC<Coord> u0(m_mat[1][1],m_mat[1][2],m_mat[1][3]);
    R3VectorTC<Coord> u1(m_mat[2][1],m_mat[2][2],m_mat[2][3]);
    R3VectorTC<Coord> u2(m_mat[3][1],m_mat[3][2],m_mat[3][3]);
    mat00 = R3Matrix(u0, u1, u2);
  }

  {
    R3VectorTC<Coord> u0(m_mat[1][0],m_mat[1][2],m_mat[1][3]);
    R3VectorTC<Coord> u1(m_mat[2][0],m_mat[2][2],m_mat[2][3]);
    R3VectorTC<Coord> u2(m_mat[3][0],m_mat[3][2],m_mat[3][3]);
    mat01 = R3Matrix(u0, u1, u2);
  }

//   R3Matrix mat01 = R3Matrix( R3Vec(m_mat[1][0],m_mat[1][2],m_mat[1][3]),
// 			     R3Vec(m_mat[2][0],m_mat[2][2],m_mat[2][3]),
// 			     R3Vec(m_mat[3][0],m_mat[3][2],m_mat[3][3]) );

  {
    R3VectorTC<Coord> u0(m_mat[1][0],m_mat[1][1],m_mat[1][3]);
    R3VectorTC<Coord> u1(m_mat[2][0],m_mat[2][1],m_mat[2][3]);
    R3VectorTC<Coord> u2(m_mat[3][0],m_mat[3][1],m_mat[3][3]);
    mat02 = R3Matrix(u0, u1, u2);
  }

//   R3Matrix mat02 = R3Matrix( R3Vec(m_mat[1][0],m_mat[1][1],m_mat[1][3]),
// 			     R3Vec(m_mat[2][0],m_mat[2][1],m_mat[2][3]),
// 			     R3Vec(m_mat[3][0],m_mat[3][1],m_mat[3][3]) );

  {
    R3VectorTC<Coord> u0(m_mat[1][0],m_mat[1][1],m_mat[1][2]);
    R3VectorTC<Coord> u1(m_mat[2][0],m_mat[2][1],m_mat[2][2]);
    R3VectorTC<Coord> u2(m_mat[3][0],m_mat[3][1],m_mat[3][2]);
    mat03 = R3Matrix(u0, u1, u2);
  }

//   R3Matrix mat03 = R3Matrix( R3Vec(m_mat[1][0],m_mat[1][1],m_mat[1][2]),
// 			     R3Vec(m_mat[2][0],m_mat[2][1],m_mat[2][2]),
// 			     R3Vec(m_mat[3][0],m_mat[3][1],m_mat[3][2]) );

   return ( + m_mat[0][0] * mat00.Determinant() 
            - m_mat[0][1] * mat01.Determinant() 
            + m_mat[0][2] * mat02.Determinant() 
            - m_mat[0][3] * mat03.Determinant()) ;
}

// Matrix Log and Exponent

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Exp( ) const
{
   R4MatrixTC<Coord> A = (*this);           // call it A to be like Alexa's pseudocode
   R4MatrixTC<Coord> X;  X.MakeIdentity();  // the current sum
   R4MatrixTC<Coord> D;  D.MakeIdentity();  // denominator
   R4MatrixTC<Coord> N;  N.MakeIdentity();  // numerator
   Coord c = 1.0;                           // coefficienty thing

   int j = (int) max(0.0, 1.0 + floor(log(A.NormF())/log(2.0)));  // gives logbase2(A.Norm())
   A = A * (Coord)pow(2.0,-j);

   int q = 6;      // supposedly 6 is a good number of iterations
   for (int k = 1; k <= q; k++) {
      c = c*(q - k + 1.0)/(Coord)(k*(2*q - k + 1.0));
      X = A*X;
      N = N + c*X;
      D = D + (Coord)pow(-1.0,k)*c*X;
   }
   
   WINbool bSuc = FALSE;
   X = D.Inverse(bSuc) * N;
   int p = (int)pow(2.0,j);
   X = X.Pow(p);
   return X;
}

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Log( /*int in_n , float &out_id*/ ) const
{
   R4MatrixTC<Coord> A = (*this);           // call it A to be like Alexa's pseudocode
   R4MatrixTC<Coord> I; I.MakeIdentity();   // identity matrix
   //A.PrintMatlab();

   int k = 0;
   int n1 = 30;
   double eps1 = 0.0001;
   while ((A-I).NormF() > eps1 /*&& k < in_n*/ ) {
      //double error = (A-I).NormF();
      A = A.Sqrt();
      k++;

      if (k > n1) {
         printf("log: repeated square roots failed to converge after %d iterations\n", n1);
         break;
      }
   }

   A = A - I;
   R4MatrixTC<Coord> Z = A;
   R4MatrixTC<Coord> X = A;
   Coord i = 1.0;
   double eps2 = 0.000000001;
   int n2 = 7;

   while ( Z.NormF() > eps2 ) {
      Z = Z*A;
      i++;
      X = X + Z/i;
      if (i > n2) { 
         printf("log: failed to converge after %d iterations\n", n2);
         break;
      }
   }

   X = (Coord)pow(2.0,k)*X;
   return X;
}

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Sqrt( ) const
{
   R4MatrixTC<Coord> A = (*this);          // call it A to be like Alexa's pseudocode
   R4MatrixTC<Coord> X = (*this);
   R4MatrixTC<Coord> Y; Y.MakeIdentity();

   WINbool bSuc = FALSE;
   int i = 0;
   double eps3 = 0.000001;
   int n3 = 10;

   while ((X*X - A).NormF() > eps3) {
      double error = (X*X - A).NormF();
      R4MatrixTC<Coord> iX = X.Inverse(bSuc);
      R4MatrixTC<Coord> iY = Y.Inverse(bSuc);
      X = (X + iY)/(Coord)2.0;
      Y = (Y + iX)/(Coord)2.0;
      i++;
      if (i > n3) { 
         if (error > 0.01)
             printf("sqrt: failed to converge, error = %f\n", error);
         break;
      }
   }
   return X;
}

// returns the Frobenius norm
template<class Coord>
inline double R4MatrixTC<Coord>::NormF( ) const
{
   double out = 0;

   for (int i = 0; i < 4; i += 1)
        for (int j = 0; j < 4; j += 1)
            out += pow(m_mat[i][j], 2);

   out = sqrt(out);
   return out;
}

// returns the 1-Norm
template<class Coord>
inline double R4MatrixTC<Coord>::Norm1( ) const
{
   double sum = 0;
   double bestsum = 0;
   for (int j = 0; j < 4; j += 1) {
      for (int i = 0; i < 4; i += 1) {
         sum += abs(m_mat[i][j]);
      }
      if (sum > bestsum) bestsum = sum;
      sum = 0;
   }
   return bestsum;
}

// returns the infinity-Norm
template<class Coord>
inline double R4MatrixTC<Coord>::NormInf( ) const
{
   double sum = 0;
   double bestsum = 0;
   for (int i = 0; i < 4; i += 1) {
      for (int j = 0; j < 4; j += 1) {
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
inline double R4MatrixTC<Coord>::Norm( ) const
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
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::Pow(int p) const
{
   R4MatrixTC<Coord> out_mat = (*this);
   R4MatrixTC<Coord> temp = (*this);

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

template<class Coord>
inline R4MatrixTC<Coord> R4MatrixTC<Coord>::MakeRandomMatrix(/*int axis, double angle*/) const 
{
    int n;
    R4MatrixTC<Coord> R; 
    R4MatrixTC<Coord> T; 
    R4MatrixTC<Coord> S; 

    // make the rotation matrix
    n = rand();
    if (n > RAND_MAX/3) {
        //int axis = rand() % 3;
        Coord v1 = (Coord)rand() / 1000.0;
        Coord v2 = (Coord)rand() / 1000.0;
        Coord v3 = (Coord)rand() / 1000.0;
        Coord angle = 1.57;
        while (angle >= 1.57) angle = (Coord)rand() / 10000.0;
        if (rand() > RAND_MAX/2) angle = -angle;
        //TRACE("angle: %f\n", angle);
        R = Rotation(R3Vec(v1, v2, v3), angle);
    }
    else R.MakeIdentity();

    // make the scale matrix
    n = rand();
    if (n > RAND_MAX/3) {
        Coord s1 = (Coord)rand() / 1000.0 + 0.01;
        Coord s2 = (Coord)rand() / 1000.0 + 0.01;
        Coord s3 = (Coord)rand() / 1000.0 + 0.01;
        //TRACE("scale: %f %f %f\n", s1, s2, s3);
        S = Scaling(s1, s2, s3, 1.0);
    }
    else S.MakeIdentity();

    // make the translation matrix
    n = rand();
    if (n > RAND_MAX/3) {
        Coord t1 = (Coord)rand() / 1000.0;
        Coord t2 = (Coord)rand() / 1000.0;
        Coord t3 = (Coord)rand() / 1000.0;
        T = Translation(R3VectorTC<Coord>(t1, t2, t3));
    }
    else T.MakeIdentity();

    //TRACE("\n");
    return S*R*T;
}

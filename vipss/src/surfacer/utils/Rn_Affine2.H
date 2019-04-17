#ifndef __BASETYPE_AFFINE2_H
#define __BASETYPE_AFFINE2_H


#include <WINSystemDefines.H>
#include <utils/Mesh_Array.H>

template <class Coord> class R2CoVectorTC;
template <class Coord> class R2PointTC;
template <class Coord> class R2MatrixTC;

extern double RNEpsilon_d;
extern float RNEpsilon_f;


/** 
\class R2VectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 2D Column vector.

  Two numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R2VectorTC
{
public:
    /**@name Constructors */
    //@{
    ///
    inline R2VectorTC() {}
    /// Initialize with given values
    inline R2VectorTC( const Coord& du, const Coord& dv );
    //@}

    /**@name Access grants */
    //@{
    /// Return dimension (always 2)
    inline int Dim() const { return 2; }

    /// Does a check that int is in bounds, i.e., 0,1
    inline Coord& operator[] ( const int i ) { ASSERT(i == 0 || i == 1); return (&du)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1
    inline const Coord& operator[] (  const int i ) const { ASSERT(i == 0 || i == 1); return (&du)[i]; }
    //@}

    /**@name Self-editing */
    //@{
    /// Divide by Length; sets length to 1.0 if not zero. Fails if length is zero.
    inline Coord Normalize();
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R2VectorTC<Coord>& operator += ( const R2VectorTC<Coord>& );
    ///
    inline R2VectorTC<Coord>& operator -= ( const R2VectorTC<Coord>& );
    ///
    inline R2VectorTC<Coord>& operator *= ( const Coord& );
    ///
    inline R2VectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    inline R2VectorTC<Coord>  operator + () const;
    ///
    inline R2VectorTC<Coord>  operator - () const;
    //@}


    /**@name Read/write/print */
    //@{
    /// Writes number [space] number [space]
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] number [space]; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) ;
    
    //@}

private:
    ///
    Coord du, dv;
};

/** 
\class R2CoVectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 2D Row vector.

  Two numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R2CoVectorTC
{
public:
    /**@name Constructors */
    //@{
    ///
    R2CoVectorTC( ) {}
    /// Initialize with given values
    R2CoVectorTC( const Coord &du, const Coord& dv );
    //@}

    /**@name Access grants */
    //@{
    /// Return dimension (always 2)
    int Dim() const { return 2; }

    /// Does a check that int is in bounds, i.e., 0,1
    Coord& operator[] (  const int i ) { ASSERT(i == 0 || i == 1); return (&du)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1
    const Coord& operator[] (  const int i ) const { ASSERT(i == 0 || i == 1); return (&du)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    R2CoVectorTC<Coord>& operator += ( const R2CoVectorTC<Coord>& );
    ///
    R2CoVectorTC<Coord>& operator -= ( const R2CoVectorTC<Coord>& );
    ///
    R2CoVectorTC<Coord>& operator *= ( const Coord& );
    ///
    R2CoVectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    R2CoVectorTC<Coord>  operator + () const;
    ///
    R2CoVectorTC<Coord>  operator - () const;
    //@}

    /**@name Read/write/print */
    //@{
    /// Writes number [space] number [space]
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] number [space]; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}

private:
    ///
    Coord du, dv;
};


/** 
\class R2PointTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 2D Point

  Two numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R2PointTC
{
public:
    /**@name Constructors */
    //@{
    ///
    R2PointTC() { }
    /// Initialize with given values
    R2PointTC(const Coord u, const Coord v);
    //@}

    /**@name Access grants */
    //@{
    /// Return dimension (always 2)
    int Dim() const { return 2; }

    /// Does a check that int is in bounds, i.e., 0,1
    Coord& operator[] (  const int i ) { ASSERT(i == 0 || i == 1); return (&u)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1
    const Coord& operator[] (  const int i ) const { ASSERT(i == 0 || i == 1); return (&u)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R2PointTC<Coord>& operator += ( const R2VectorTC<Coord>& );
    ///
    inline R2PointTC<Coord>& operator -= ( const R2VectorTC<Coord>& );
    //@}

    /**@name Comparison operators */
    //@{
    /// 
    bool operator < ( const R2PointTC<Coord>& p ) const;
    ///
    bool operator <= ( const R2PointTC<Coord>& p ) const;
    //@}

    /**@name Read/write/print */
    //@{
    /// Writes number [space] number [space]
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] number [space]; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}

private:
    ///
    Coord u, v;
};

/** 
\class R2RayTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 2D Column vector.

  Two numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template< class Coord >
class R2RayTC
{
private:
    ///
    R2PointTC<Coord> pOrg;
    ///
    R2VectorTC<Coord> vDir;

public:
    /**@name Access grants */
    //@{
    ///
    const R2PointTC<Coord> &Pt() const { return pOrg; }
    ///
    const R2VectorTC<Coord> &Dir() const { return vDir; }
    //@}



    /**@name Assignment operators */
    //@{
    ///
    R2PointTC<Coord> &Pt() { return pOrg; }
    ///
    R2VectorTC<Coord> &Dir() { return vDir; }
    //@}

    /**@name Constructor/destructor */
    //@{
    ///
    R2RayTC() { }
    ///
    /// Initialize with given elements
    R2RayTC( const R2PointTC<Coord> & _pOrg, const R2VectorTC<Coord>& _vDir )
        { pOrg = _pOrg; vDir = _vDir; }
    ///
    ~R2RayTC() {}
    //@}

    /**@name Read/write/print */
    //@{
    /// Writes point vector
    inline void    Write(ofstream &out) const { pOrg.Write(out); vDir.Write(out); }
    ///
    inline WINbool Read(ifstream &in) { pOrg.Read(in); return vDir.Read(in); }
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is point [space] vector [space]; if true, replaces spaces with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const { pOrg.Print( in_bDoReturn ); vDir.Print( in_bDoReturn ); }
    //@}

};

/** 
\class R2MatrixTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 2X2 matrix

  Four numbers.  Type-safe addition/subtraction/multiplication for matrices, points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template <class T> class C2PointTC;
template<class Coord>
class R2MatrixTC
{
public:
    /**@name Constructors */
    //@{
    ///
    inline R2MatrixTC() { mat[0][0] = mat[1][1] = 1.0; mat[0][1] = mat[1][0] = 0.0; }
    ///
    /// Initialize columns with given vectors
    inline R2MatrixTC( const R2VectorTC<Coord>& vCol1, const R2VectorTC<Coord>& vCol2 );
    /// Initialize rows with given covectors
    inline R2MatrixTC( const R2CoVectorTC<Coord>& covRow1, const R2CoVectorTC<Coord>& covRow2 );
    /// Row-column order
    inline R2MatrixTC( const Coord d00, const Coord d01, const Coord d10, const Coord d11 )
    {
        mat[0][0] = d00; mat[0][1] = d01; mat[1][0] = d10; mat[1][1] = d11;
    }
    //@}

    /**@name Access grants */
    //@{
    ///
    inline Coord & operator () ( const int i, const int j ) { return mat[i][j]; }
    ///
    inline const Coord & operator () ( const int i, const int j ) const { return mat[i][j]; }
    //@}

    /**@name Assignment operators 
    */
    //@{
    /// 
    inline R2MatrixTC<Coord> & operator *= ( const R2MatrixTC<Coord>& );
    /// Componant-wise addition
    inline R2MatrixTC<Coord> & operator += ( const R2MatrixTC<Coord>& );
    /// Componant-wise subtraction
    inline R2MatrixTC<Coord> & operator -= ( const R2MatrixTC<Coord>& );
    /// Componant-wise multiplication
    inline R2MatrixTC<Coord> & operator *= ( const Coord );
    /// Componant-wise division (will not divide if input is zero)
    inline R2MatrixTC<Coord> & operator /= ( const Coord );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    inline R2MatrixTC<Coord> operator + ( );
    ///
    inline R2MatrixTC<Coord> operator - ( );
    //@}

    /**@name Matrix operation entry points.
    */
    //@{
    /// Sets self to identity matrix
    inline void SetIdentity(  );
    /// Returns the inverse of this matrix, if exists; @param out_dfDeterminant is the determinant
    inline R2MatrixTC<Coord> Inverse( Coord &out_dfDeterminant ) const;
    ///
    inline Coord Determinant () const;
    /// Angle for a Givens rotation of the matrix (the rotation component)
    inline Coord Givens() const;
    ///
    inline Coord Trace () const { return mat[0][0] + mat[1][1]; }
    /** \brief The eigen vectors and eigen values
      * This only seems to work some of the time. If you want a more stable version, see fitting/FITTools_SVD
      @param out_cptL1 The first Eigen (largest) value, may be complex
      @param out_cptL2 The second Eigen (smallest) value, may be complex
      @param out_vec1 The first Eigen vector, may be complex
      @param out_vec2 The second Eigen vector, may be complex
      @returns The number of valid Eigen Values */
    inline int Eigen( C2PointTC<Coord> &out_cptL1, C2PointTC<Coord> &out_cptL2, R2VectorTC< C2PointTC<Coord> > & out_vec1, R2VectorTC< C2PointTC< Coord> > & out_vec2 ) const;
    /** \brief The eigen vectors and eigen values
      * Assumes matrix has positive, real eigen values and vectors
      @param out_dE1 The first Eigen (largest) value,
      @param out_dE2 The second Eigen (smallest) value
      @param out_vec1 The first Eigen vector
      @param out_vec2 The second Eigen vector
      @returns The number of valid Eigen Values */
    inline int EigenReal( Coord &out_dE1, Coord &out_dE2, R2VectorTC< Coord > & out_vec1, R2VectorTC< Coord > & out_vec2, const Coord in_dEps = 1e-6 ) const;
    /** \brief Singular value decomposition
      * This only seems to work some of the time. If you want a more stable version, see fitting/FITTools_SVD
      @param out_matU A 2X2 orthonormal matrix with the two eigen vectors
      @param out_matDiag A 2X2 diagonal matrix with the eigen values
      @param out_matVTransp A 2X2 orthonormal matrix which is the transpose of the v vectors
      @returns True if two valid eigen values/vectors */
    inline WINbool SVD( R2MatrixTC<Coord> &out_matU, R2MatrixTC<Coord> &out_matDiag, R2MatrixTC<Coord> &out_matVTransp ) const;
    ///
    inline void QRDecomposition( R2MatrixTC<Coord> &out_matQ, R2MatrixTC<Coord> &out_matR ) const;
    /// Only works for symmetric matrices
    inline WINbool Takagi( R2MatrixTC<Coord> &out_matU, R2MatrixTC<Coord> &out_matDiag ) const;
    ///
    inline R2MatrixTC<Coord> Transpose(  ) const;
    /// Set self to be a 2D rotation matrix
    inline WINbool SetRotation( const Coord fAngle );
    /// Set self to be a 2D scale matrix
    inline WINbool SetScale( const Coord fScaleU, const Coord fScaleV );
    /// Static member function returning a 2D matrix that performs the given rotation
    static R2MatrixTC<Coord> Rotation(const Coord angle);
    /// Static member function returning a 2D matrix that performs the given scaling
    static R2MatrixTC<Coord> Scaling(const Coord  fScaleU, const Coord fScaleV);
    //@}

    /**@name Matrix log and exponential */
    //@{
    /// Exponential. Found by optimization.
    inline R2MatrixTC<Coord> Exp( ) const;
    /// Log. Found by optimization
    inline R2MatrixTC<Coord> Log( ) const;
    /// Squareroot. Found by optimization
    inline R2MatrixTC<Coord> Sqrt( ) const;
    /// the Frobenius norm
    inline Coord NormF( ) const;
    /// the 1-norm
    inline Coord Norm1( ) const;
    /// the infinity-norm
    inline Coord NormInf( ) const;
    /// an approximated 2-norm
    inline Coord Norm( ) const;
    /// Raise matrix to given power
    inline R2MatrixTC<Coord> Pow(int) const;
    //@}

    /**@name Read/write/print */
    //@{
    /// Writes 2 X 2 matrix, formatted
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is formatted; if true, replaces carriage returns with spaces */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;
    //@}

private:
    ///
    R2VectorTC< C2PointTC< Coord> > GetEigenVec( const C2PointTC<Coord> &in_dV ) const;
    ///
    Coord mat[2][2];
};

/** \ingroup SimpleGeometryElems Typedefs */
//@{
/// From R2PointTC
typedef R2PointTC<double> R2Pt;
/// From R2PointTC
typedef R2PointTC<float>  R2Pt_f;
/// From R2PointTC
typedef R2PointTC<int>    R2Pt_i;

/// From R2VectorTC
typedef R2VectorTC<double> R2Vec;
/// From R2VectorTC
typedef R2VectorTC<float>  R2Vec_f;
/// From R2VectorTC
typedef R2VectorTC<int>    R2Vec_i;

/// From R2CoVectorTC
typedef R2CoVectorTC<double> R2CoVec;

/// From R2RayTC
typedef R2RayTC<double> R2Ray;

/// From R2MatrixTC
typedef R2MatrixTC<double> R2Matrix;
//@}

#include "Rn_Vector2_i.H"
#include "Rn_CoVector2_i.H"
#include "Rn_Point2_i.H"
#include "Rn_Matrix2_i.H"
#include "Rn_Unary2_i.H"
#include "Rn_Binary2_i.H"
#include "Rn_Io2_i.H"

#endif
// #ifndef __BASETYPE_AFFINE2_H


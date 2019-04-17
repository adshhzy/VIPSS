#if _MSC_VER > 1000
#pragma once
#endif
// _MSC_VER > 1000

#ifndef __BASETYPE_AFFINE3_H
#define __BASETYPE_AFFINE3_H

#include <WINSystemDefines.H>

#include <utils/Rn_Affine2.H>

template <class Coord> class R3CoVectorTC;
template <class Coord> class R3PointTC;
template <class Coord> class R3MatrixTC;

/** 
\class R3VectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 3D Column vector.

  Three numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R3VectorTC
{
public:
    ///
    R3VectorTC() {}
    /// Initialize with given values
    inline R3VectorTC( const Coord& dx, const Coord& dy, const Coord& dz );

    /// Return dimension (always 3)
    int Dim() const { return 3; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0,1,2
    Coord& operator[] ( const int i ) { ASSERT(i >= 0 && i <= 2 ); return (&dx)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1,2
    const Coord& operator[] ( const int i ) const { ASSERT(i >= 0 && i <= 2 ); return (&dx)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R3VectorTC<Coord>& operator += ( const R3VectorTC<Coord>& );
    ///
    inline R3VectorTC<Coord>& operator -= ( const R3VectorTC<Coord>& );
    ///
    inline R3VectorTC<Coord>& operator *= ( const Coord& );
    ///
    inline R3VectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Self-editing */
    //@{
    /// Divide by Length; sets length to 1.0 if not zero. Fails if length is zero.
    inline Coord Normalize();
    //@}

    /**@name Unary operators */
    //@{
    /// 
    inline R3VectorTC<Coord>  operator + () const;
    ///
    inline R3VectorTC<Coord>  operator - () const;
    /// Return a vector that is perpendicular to this one (not unique)
    inline R3VectorTC<Coord>  Perpendicular() const;
    //@}

    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] ...; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}
private:
    ///
    Coord dx, dy, dz;
};

/** 
\class R3CoVectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 3D Row vector.

  Three numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R3CoVectorTC
{
public:
    ///
    R3CoVectorTC( ) {}
    ///
    /// Initialize with given values
    R3CoVectorTC( const Coord& dx, const Coord& dy, const Coord& dz );

    /// Return dimension (always 3)
    int Dim() const { return 3; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0,1,2
    Coord& operator[] ( const int i ) { ASSERT(i >= 0 && i <= 2 ); return (&dx)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1,2
    const Coord& operator[] ( const int i ) const { ASSERT(i >= 0 && i <= 2 ); return (&dx)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    R3CoVectorTC<Coord>& operator += ( const R3CoVectorTC<Coord>& );
    ///
    R3CoVectorTC<Coord>& operator -= ( const R3CoVectorTC<Coord>& );
    ///
    R3CoVectorTC<Coord>& operator *= ( const Coord& );
    ///
    R3CoVectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    R3CoVectorTC<Coord>  operator + () const;
    ///
    R3CoVectorTC<Coord>  operator - () const;
    //@}


    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] ...; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}
private:
    ///
    Coord dx, dy, dz;
};


/** 
\class R3PointTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 3D point.

  Three numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R3PointTC
{
public:
    /// Default constructor
    R3PointTC() { }
    /// Initialize with given values
    inline R3PointTC( const Coord& x, const Coord& y, const Coord& z );

    /// Return dimension (always 3)
    int Dim() const { return 3; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0,1,2
    Coord& operator[] ( const int i ) { ASSERT(i >= 0 && i <= 2 ); return (&x)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1,2
    const Coord& operator[] ( const int i ) const { ASSERT(i >= 0 && i <= 2 ); return (&x)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R3PointTC<Coord>& operator += ( const R3VectorTC<Coord>& );
    ///
    inline R3PointTC<Coord>& operator -= ( const R3VectorTC<Coord>& );
    //@}


    /**@name Point dominance */
    //@{
    /// 
    inline WINbool operator < ( const R3PointTC<Coord>& ) const;
    ///
    inline WINbool operator <= ( const R3PointTC<Coord>& ) const;
    //@}

    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] ...; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}
private:
    ///
    Coord x, y, z;
};


/** 
\class R3RayTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 3D ray (point plus vector)

  Three numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template< class Coord >
class R3RayTC
{
private:
    ///
    R3PointTC<Coord> pOrg;
    ///
    R3VectorTC<Coord> vDir;

public:
    /**@name Access grants */
    //@{
    ///
    const R3PointTC<Coord> &Pt() const { return pOrg; }
    ///
    const R3VectorTC<Coord> &Dir() const { return vDir; }
    //@}

    /**@name Assignment operators */
    //@{
    ///
    R3PointTC<Coord> &Pt() { return pOrg; }
    ///
    R3VectorTC<Coord> &Dir() { return vDir; }
    ///
    R3PointTC<Coord> operator()( const double in_dT ) const { return pOrg + vDir * in_dT; }
    //@}

    /**@name Constructor/destructor */
    //@{
    ///
    R3RayTC() { }
    ///
    /// Initialize with given elements
    R3RayTC( const R3PointTC<Coord> & _pOrg, const R3VectorTC<Coord>& _vDir )
        { pOrg = _pOrg; vDir = _vDir; }
    ///
    ~R3RayTC() {}
    //@}

    /**@name Read/write/print */
    //@{
    ///Write out in format: Pt Dir
    inline void    Write(ofstream &out) const { pOrg.Write(out); vDir.Write(out); }
    ///Read in format: Pt Dir
    inline WINbool Read(ifstream &in) { pOrg.Read(in); return vDir.Read(in); }
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) { pOrg.ReadBinary(in); return vDir.ReadBinary(in); }
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const { pOrg.WriteBinary(out); vDir.WriteBinary(out); }
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is spaces between pt and vec ...; if true, replaces spaces with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const { pOrg.Print( in_bDoReturn ); vDir.Print( in_bDoReturn ); }
    //@}

};

/** 
\class R3MatrixTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 3X3 matrix.

  Type-safe addition/subtraction/multiplication for matrices,,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template< class Coord >
class R3MatrixTC
{
public:
    ///
    R3MatrixTC() { m_mat[0][0] = m_mat[1][1] = m_mat[2][2] = 1.0; m_mat[0][1] = m_mat[0][2] = m_mat[1][0] = m_mat[1][2] = m_mat[2][0] = m_mat[2][1] = 0.0; }
    ///
    /// Initialize columns with given vectors
    inline R3MatrixTC( const R3VectorTC<Coord>& vCol1, const R3VectorTC<Coord>& vCol2, const R3VectorTC<Coord>& vCol3 );
    ///
    /// Initialize rows with given covectors
    inline R3MatrixTC( const R3CoVectorTC<Coord>& cvRow1, const R3CoVectorTC<Coord>& cvRow2, const R3CoVectorTC<Coord>& cvRow3 );
    /// Convert an R2Matrix to an R3 matrix by adding an identity row/column
    inline R3MatrixTC( const R2MatrixTC<Coord> &in_oMat );

    /**@name Access grants */
    //@{
    ///
    Coord & operator () ( const int i, const int j ) { return m_mat[i][j]; }
    ///
    const Coord & operator () ( const int i, const int j ) const { return m_mat[i][j]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R3MatrixTC<Coord>& operator *= ( const R3MatrixTC<Coord>& );
    ///
    inline R3MatrixTC<Coord>& operator += ( const R3MatrixTC<Coord>& );
    ///
    inline R3MatrixTC<Coord>& operator -= ( const R3MatrixTC<Coord>& );
    ///
    inline R3MatrixTC<Coord>& operator *= ( const Coord );
    ///
    inline R3MatrixTC<Coord>& operator /= ( const Coord );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    inline R3MatrixTC<Coord> operator + ( );
    ///
    inline R3MatrixTC<Coord> operator - ( );
    //@}

    /**@name Matrix construction operators */
    //@{
    /// Set self to identity
    void MakeIdentity( );
    /// Static method that returns an identity matrix
    static R3MatrixTC<Coord> Identity( );
    /** \brief Static method that returns a rotation matrix 
      @param axis is one of 0,1,2 (x,y,z)
      @param angle Angle in Radians
      @returns matrix */
    static R3MatrixTC<Coord> Rotation(const int axis, const Coord angle);
    /** \brief Static method that returns a rotation matrix around the vector by the given angle

     Will return a bad result if vector is zero length
     @param in_vec Vector (vector will be normalized)
     @param angle Angle in Radians
     */
    static R3MatrixTC<Coord> Rotation(const R3VectorTC<Coord> &in_vec, const Coord angle);
    /** \brief Static method with 2x2 rotation matrix in upper left
     @param angle 2D rotation angle (i.e., rotation around z axis)
     */
    static R3MatrixTC<Coord> Rotation2D( const Coord angle);
    /// Static method that returns a scale matrix
    static R3MatrixTC<Coord> Scaling( const Coord in_x, const Coord in_y, const Coord in_z );
    /// Static method that returns a uniform scaling matrix
    static R3MatrixTC<Coord> Scaling( const Coord in_scl );
    /// Static method that returns a translation matrix for a 2D point/vector
    static R3MatrixTC<Coord> Translation( const R2VectorTC<Coord> & );

    /** \brief Matrix taking the from vectors to the to vectors

     Produces a rotation matrix that takes the first two vectors to the second two. Assumes from (and to) vectors are orthogonal;
     constructs the third vector as the cross product of the first two.
     Fails if the input vectors are not unit length or not perpendicular
     @param out_mat The rotation matrix
     @param in_vFromX First input vector of pair one. Should be unit length.
     @param in_vFromY Second input vector of pair one. Should be unit length and perpendicular to in_vFromX
     @param in_vToX First input vector of pair two. Should be unit length.
     @param in_vToY Second input vector of pair two. Should be unit length and perpendicular to in_vToX
     @returns FALSE if no valid rotation matrix */
    static WINbool MatrixVecToVec( R3MatrixTC<Coord>     & out_mat, 
                                   const R3VectorTC<Coord>  & in_vFromX, 
                                   const R3VectorTC<Coord>  & in_vFromY,
                                   const R3VectorTC<Coord>  & in_vToX,   
                                   const R3VectorTC<Coord>  & in_vToY );

    /** \brief Matrix taking the from vector to the to vector, 
    
     Rotation around normal unspecified
     Produces a rotation matrix that takes the first vector to the second. Assumes both vectors are unit length. 
     Builds two perpendicular vectors from the input ones.
     @param out_mat Rotation matrix taking in_vFromX to in_vToX
     @param in_vFromX Unit length vector
     @param in_vToX Unit length vector
     @returns FALSE if vectors not unit length*/
    static WINbool MatrixVecToVec( R3MatrixTC<Coord>     & out_mat, 
                                   const R3VectorTC<Coord>  & in_vFromX, 
                                   const R3VectorTC<Coord>  & in_vToX );
    //@}

    /**@name Matrix operations */
    //@{
    /** \brief Take the inverse of the matrix
      @param out_bool If the determinant was non-zero, i.e., inverse exists
      @returns matrix inverse */
    inline R3MatrixTC<Coord> Inverse( WINbool &out_bool ) const;
    ///
    inline Coord Determinant () const;
    ///
    inline Coord Trace () const { return m_mat[0][0] + m_mat[1][1] + m_mat[2][2]; }
    /// 
    inline R3MatrixTC<Coord> Transpose(  ) const;
    /** \brief The eigen vectors and eigen values

      * This only seems to work some of the time. If you want a more stable version, see fitting/FITTools_SVD
      @param out_cptL1 The first (largest) Eigen value, may be complex
      @param out_cptL2 The second Eigen value, may be complex
      @param out_cptL3 The third (smallest) Eigen value, may be complex
      @param out_vec1 The first Eigen vector, may be complex
      @param out_vec2 The second Eigen vector, may be complex
      @param out_vec3 The third Eigen vector, may be complex
      @returns The number of valid Eigen Values */
    int Eigen( C2PointTC< Coord > &out_cptL1, C2PointTC< Coord > &out_cptL2, C2PointTC< Coord > &out_cptL3, 
               R3VectorTC< C2PointTC< Coord > > & out_vec1, 
               R3VectorTC< C2PointTC< Coord > > & out_vec2, 
               R3VectorTC< C2PointTC< Coord > > & out_vec3 ) const;
    /// Only works for symmetric matrices
    WINbool Takagi( R3MatrixTC<Coord> &out_matU, R3MatrixTC<Coord> &out_matDiag ) const;
    //@}

    /**@name Matrix log and exponential */
    //@{
    /// Exponential. Found by optimization.
    inline R3MatrixTC<Coord> Exp( WINbool &out_bSuc ) const;
    /// Log. Found by optimization
    inline R3MatrixTC<Coord> Log( WINbool &out_bSuc ) const;
    /// Squareroot. Found by optimization
    inline R3MatrixTC<Coord> Sqrt( WINbool &out_bSuc ) const;
    /// the Frobenius norm
    inline double NormF( ) const;
    /// the 1-norm
    inline double Norm1( ) const;
    /// the infinity-norm
    inline double NormInf( ) const;
    /// an approximated 2-norm
    inline double Norm( ) const;
    /// Raise matrix to given power
    inline R3MatrixTC<Coord> Pow(int) const;
    //@}

    /**@name Read/write/print */
    //@{
    /// Writes 3 X 3 matrix, formatted
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline WINbool ReadBinary(ifstream &in) ;
    /// Will not work across platforms that have different floating point conventions (big vs little indian)
    inline void    WriteBinary(ofstream &out) const;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is formatted; if true, replaces carriage returns with spaces< */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;
    //@}

private:
    ///
    R3VectorTC< C2PointTC<Coord> > GetEigenVec( const C2PointTC<Coord> &in_dV ) const;
    ///
    Coord m_mat[3][3];
};

/** \ingroup SimpleGeometryElems Typedefs */
//@{
/// From R3PointTC
typedef R3PointTC<double> R3Pt;
/// From R3PointTC
typedef R3PointTC<float>  R3Pt_f;
/// From R3PointTC
typedef R3PointTC<int>    R3Pt_i;

/// From R3VectorTC
typedef R3VectorTC<double> R3Vec;
/// From R3VectorTC
typedef R3VectorTC<float>  R3Vec_f;
/// From R3VectorTC
typedef R3VectorTC<int>    R3Vec_i;

/// From R3CoVectorTC
typedef R3CoVectorTC<double> R3CoVec;

/// From R3RayTC
typedef R3RayTC<double> R3Ray;

/// From R3MatrixTC
typedef R3MatrixTC<double> R3Matrix;
//@}

#include "Rn_Vector3_i.H"
#include "Rn_CoVector3_i.H"
#include "Rn_Point3_i.H"
#include "Rn_Matrix3_i.H"
#include "Rn_Unary3_i.H"
#include "Rn_Binary3_i.H"
#include "Rn_Io3_i.H"

#endif
// #ifndef __BASETYPE_AFFINE3_H

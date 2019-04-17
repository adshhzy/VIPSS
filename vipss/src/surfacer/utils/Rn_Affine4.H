#if _MSC_VER > 1000
#pragma once
#endif
// _MSC_VER > 1000

#ifndef __BASETYPE_AFFINE4_H
#define __BASETYPE_AFFINE4_H

#include <WINSystemDefines.H>

template <class Coord> class R4CoVectorTC;
template <class Coord> class R4PointTC;
template <class Coord> class R4MatrixTC;

/** 
\class R4VectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 4D Column vector.

  Four numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R4VectorTC
{
public:
    ///
    R4VectorTC() {}
    /// Initialize with given values
    inline R4VectorTC( const Coord& , const Coord& , const Coord& , const Coord & );

    /// Return dimension (always 4)
    int Dim() const { return 4; }


    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0,1,2,3
    Coord& operator[] ( const int i ) { ASSERT(i >= 0 && i <= 3 ); return (&dx)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1,2,3
    const Coord& operator[] ( const int i ) const { ASSERT(i >= 0 && i <= 3 ); return (&dx)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R4VectorTC<Coord>& operator += ( const R4VectorTC<Coord>& );
    ///
    inline R4VectorTC<Coord>& operator -= ( const R4VectorTC<Coord>& );
    ///
    inline R4VectorTC<Coord>& operator *= ( const Coord& );
    ///
    inline R4VectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Self-editing */
    //@{
    /// Divide by Length; sets length to 1.0 if not zero. Fails if length is zero.
    inline Coord Normalize();
    //@}

    /**@name Unary operators */
    //@{
    /// 
    inline R4VectorTC<Coord>  operator + () const;
    ///
    inline R4VectorTC<Coord>  operator - () const;
    //@}


    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] ...; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;
    //@}

protected:
    /// Usually used as a homogeneous R3 point
    Coord dx, dy, dz, dw;
};

/** 
\class R4CoVectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 4D Row vector.

  Four numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R4CoVectorTC
{
public:
    ///
    R4CoVectorTC( ) {}

    /// Initialize with given values
    R4CoVectorTC( const Coord& dx, const Coord& dy, const Coord& dz, const Coord& dw );

    /// Return dimension (always 4)
    int Dim() const { return 4; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0,1,2,3
    Coord& operator[] ( const int i ) { ASSERT(i >= 0 && i <= 3 ); return (&dx)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1,2,3
    const Coord& operator[] ( const int i ) const { ASSERT(i >= 0 && i <= 3 ); return (&dx)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    R4CoVectorTC<Coord>& operator += ( const R4CoVectorTC<Coord>& );
    ///
    R4CoVectorTC<Coord>& operator -= ( const R4CoVectorTC<Coord>& );
    ///
    R4CoVectorTC<Coord>& operator *= ( const Coord& );
    ///
    R4CoVectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    R4CoVectorTC<Coord>  operator + () const;
    ///
    R4CoVectorTC<Coord>  operator - () const;
    //@}


    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] ...; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}
private:
    ///
    Coord dx, dy, dz, dw;
};


/** 
\class R4PointTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 4D point

  Four numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R4PointTC
{
public:
    ///
    R4PointTC() { }
    ///
    /// Initialize with given values
    inline R4PointTC( const Coord& x, const Coord& y, const Coord& z, const Coord & w );

    /// Return dimension (always 4)
    int Dim() const { return 4; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0,1,2,3
    Coord& operator[] ( const int i ) { ASSERT(i >= 0 && i <= 3 ); return (&x)[i]; }
    /// Does a check that int is in bounds, i.e., 0,1,2,3
    const Coord& operator[] ( const int i ) const { ASSERT(i >= 0 && i <= 3 ); return (&x)[i]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R4PointTC<Coord>& operator += ( const R4VectorTC<Coord>& );
    ///
    inline R4PointTC<Coord>& operator -= ( const R4VectorTC<Coord>& );
    //@}


    /**@name Point dominance */
    //@{
    /// 
    inline WINbool operator < ( const R4PointTC<Coord>& ) const;
    ///
    inline WINbool operator <= ( const R4PointTC<Coord>& ) const;
    //@}

    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space] ...; if true, replaces last space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}
private:
    ///
    Coord x, y, z, w;
};


/** 
\class R4RayTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 4D Point and Ray.

  Four numbers.  Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template< class Coord >
class R4RayTC
{
private:
    ///
    R4PointTC<Coord> pOrg;
    ///
    R4VectorTC<Coord> vDir;

public:
    /**@name Access grants */
    //@{
    ///
    const R4PointTC<Coord> &Pt() const { return pOrg; }
    ///
    const R4VectorTC<Coord> &Dir() const { return vDir; }
    //@}

    /**@name Assignment operators */
    //@{
    ///
    R4PointTC<Coord> &Pt() { return pOrg; }
    ///
    R4VectorTC<Coord> &Dir() { return vDir; }
    //@}

    /**@name Constructor/destructor */
    //@{
    ///
    R4RayTC() { }
    ///
    /// Initialize with given elements
    R4RayTC( const R4PointTC<Coord> & _pOrg, const R4VectorTC<Coord>& _vDir )
        { pOrg = _pOrg; vDir = _vDir; }
    ///
    ~R4RayTC() {}
    //@}

    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const { pOrg.Write(out); vDir.Write(out); }
    ///
    inline WINbool Read(ifstream &in) { pOrg.Read(in); return vDir.Read(in); }
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is spaces between pt and vec ...; if true, replaces spaces with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const { pOrg.Print( in_bDoReturn ); vDir.Print( in_bDoReturn ); }
    //@}

};

/** 
\class R4MatrixTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 4x4 matrix

  Four numbers. Data must be in column major order to match OpenGL Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template< class Coord >
class R4MatrixTC
{
public:
    ///
    R4MatrixTC() {}
    /// Initialize columns with given vectors
    inline R4MatrixTC( const R4VectorTC<Coord>& vCol1, 
                       const R4VectorTC<Coord>& vCol2, 
                       const R4VectorTC<Coord>& vCol3,
                       const R4VectorTC<Coord>& vCol4 );
    /// Initialize rows with given covectors
    inline R4MatrixTC( const R4CoVectorTC<Coord>& cvRow1, 
                       const R4CoVectorTC<Coord>& cvRow2, 
                       const R4CoVectorTC<Coord>& cvRow3,
                       const R4CoVectorTC<Coord>& cvRow4 );

    /// Convert an R3Matrix to an R4 matrix by adding an identity row/column
    inline R4MatrixTC( const R3MatrixTC<Coord> &in_oMat );

    /**@name Access grants */
    //@{
    ///
    inline Coord & operator () ( const int iRow, const int iCol ) { return m_mat[iCol][iRow]; }
    ///
    inline const Coord & operator () ( const int iRow, const int iCol ) const { return m_mat[iCol][iRow]; }
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    inline R4MatrixTC<Coord>& operator *= ( const R4MatrixTC<Coord>& );
    ///
    inline R4MatrixTC<Coord>& operator += ( const R4MatrixTC<Coord>& );
    ///
    inline R4MatrixTC<Coord>& operator -= ( const R4MatrixTC<Coord>& );
    ///
    inline R4MatrixTC<Coord>& operator *= ( const Coord );
    ///
    inline R4MatrixTC<Coord>& operator /= ( const Coord );
    //@}

    /**@name Unary operators */
    //@{
    /// 
    inline R4MatrixTC<Coord> operator + ( );
    ///
    inline R4MatrixTC<Coord> operator - ( );
    //@}

    /**@name Matrix construction operators */
    //@{
    /// Set self to identity
    void MakeIdentity( );
    /// Static method that returns an identity matrix
    static R4MatrixTC<Coord> Identity( );
    /** \brief Static method that returns a rotation matrix 
      @param axis is one of 0,1,2 (x,y,z)
      @param angle Angle in Radians
      @returns matrix */
    static R4MatrixTC<Coord> Rotation(const int axis, const Coord angle);
    /** \brief Static method that returns a rotation matrix around the vector by the given angle

     Will return a bad result if vector is zero length
     @param in_vec Vector (vector will be normalized)
     @param angle Angle in Radians
     */
    static R4MatrixTC<Coord>Rotation(const R3VectorTC<Coord> &in_vec, const Coord angle);
    /// Static method that returns a scale matrix
    static R4MatrixTC<Coord> Scaling( const Coord in_x, const Coord in_y, const Coord in_z, const Coord in_w );
    /// Static method that returns a uniform scaling matrix
    static R4MatrixTC<Coord> Scaling( const Coord in_scl );
    /// Static method that returns a translation matrix for a 3D point/vector
    static R4MatrixTC<Coord> Translation( const R3VectorTC<Coord> &in_vec );
    /// Static method that returns a translation matrix for a 3D point/vector
    static R4MatrixTC<Coord> Translation( const R3PointTC<Coord> &in_pt );
    //@}

    /**@name Matrix operations */
    //@{
    /** \brief Take the inverse of the matrix
      @param out_bSuccess If the determinant was non-zero, i.e., inverse exists
      @returns matrix inverse */
    inline R4MatrixTC<Coord> Inverse( WINbool &out_bSuccess ) const;
    /// 
    inline R4MatrixTC<Coord> Transpose(  ) const;
	///
    inline Coord Determinant () const;
    //@}

    /**@name Matrix log and exponential */
    //@{
    /// Exponential. Found by optimization.
    inline R4MatrixTC<Coord> Exp( ) const;
    /// Log. Found by optimization
    inline R4MatrixTC<Coord> Log(  ) const;
    /// Squareroot. Found by optimization
    inline R4MatrixTC<Coord> Sqrt( ) const;
    /// the Frobenius norm
    inline double NormF( ) const;
    /// the 1-norm
    inline double Norm1( ) const;
    /// the infinity-norm
    inline double NormInf( ) const;
    /// an approximated 2-norm
    inline double Norm( ) const;
    /// Raise matrix to given power
    inline R4MatrixTC<Coord> Pow(int) const;
    ///
    inline R4MatrixTC<Coord> MakeRandomMatrix() const;
    //@}


    /**@name Read/write/print */
    //@{
    /// Writes 3 X 3 matrix, formatted
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is formatted; if true, replaces spaces with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;
    ///
    inline void    PrintMatlab( ) const;
    //@}

private:
    /// column row
    Coord m_mat[4][4];
};

/* -------------------  Typedefs  ------------------------------- */

/** \ingroup SimpleGeometryElems Typedefs */
//@{
/// From R4PointTC
typedef R4PointTC<double> R4Pt;
/// From R4PointTC
typedef R4PointTC<float>  R4Pt_f;
/// From R4PointTC
typedef R4PointTC<int>    R4Pt_i;

/// From R4VectorTC
typedef R4VectorTC<double> R4Vec;
/// From R4VectorTC
typedef R4VectorTC<float>  R4Vec_f;
/// From R4VectorTC
typedef R4VectorTC<int>    R4Vec_i;

/// From R4CoVectorTC
typedef R4CoVectorTC<double> R4CoVec;

/// From R4RayTC
typedef R4RayTC<double> R4Ray;

/// From R4MatrixTC
typedef R4MatrixTC<double> R4Matrix;

//@}

#include "Rn_Vector4_i.H"
#include "Rn_CoVector4_i.H"
#include "Rn_Point4_i.H"
#include "Rn_Matrix4_i.H"
#include "Rn_Unary4_i.H"
#include "Rn_Binary4_i.H"
#include "Rn_Io4_i.H"

#endif
// #ifndef __BASETYPE_AFFINE4_H

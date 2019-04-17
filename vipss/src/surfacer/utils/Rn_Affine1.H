#ifndef __BASETYPE_AFFINE1_H
#define __BASETYPE_AFFINE1_H

#include <WINSystemDefines.H>


template <class Coord> class R1CoVectorTC;
template <class Coord> class R1PointTC;

extern double RNEpsilon_d;
extern float RNEpsilon_f;


/** 
\class R1VectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems
\brief A 1D Column vector.

  Essentially a single number. Included for completeness'
  sake. Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R1VectorTC
{
public:
    /**@name Constructors */
    //@{
    ///
    inline R1VectorTC() {}
    /// Initialize with value
    inline R1VectorTC( const Coord d );
    //@}
    
    /**@name Access */
    //@{

    /// Return dimension (always returns 1)
    inline int Dim() const { return 1; }

    /// Does a check that int is in bounds, i.e., 0
    inline Coord& operator[] ( const int DECLAREVAR(i) ) { ASSERT(i == 0); return d; }
    /// Does a check that int is in bounds, i.e., 0 
    inline const Coord& operator[] ( const int DECLAREVAR(i) ) const { ASSERT(i == 0); return d; }
    //@}

    /**@name Assignment operators */
    //@{
    ///
    inline R1VectorTC<Coord>& operator += ( const R1VectorTC<Coord>& );
    ///
    inline R1VectorTC<Coord>& operator -= ( const R1VectorTC<Coord>& );
    ///
    inline R1VectorTC<Coord>& operator *= ( const Coord& );
    ///
    inline R1VectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Self-editing */
    //@{
    /// Divide by self; sets to 1.0 if not zero. Fails if value is zero.
    inline double Normalize();
    //@}

    /**@name Unary operators */
    //@{
    ///
    inline R1VectorTC<Coord>  operator + () const;
    ///
    inline R1VectorTC<Coord>  operator - () const;
    //@}

    /**@name Read/write/print */
    //@{
    /// Writes number [space]
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space]; if true, replaces space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}

private:
    ///
    Coord d;
};


/** 
\brief A 1D Row vector.
\class R1CoVectorTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems


  Essentially a single number. Included for completeness'
  sake. Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R1CoVectorTC
{
public:
    ///
    inline R1CoVectorTC( ) {}
    /// Initialize with value
    inline R1CoVectorTC( const Coord d );

    /// Return dimension of the type
    inline int Dim() const { return 1; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0 
    inline Coord& operator[] ( const int DECLAREVAR(i) ) { ASSERT(i == 0); return d; }
    /// Does a check that int is in bounds, i.e., 0 
    inline const Coord& operator[] ( const int DECLAREVAR(i) ) const { ASSERT(i == 0); return d; }
    //@}

    /**@name Assignment operators */
    //@{
    ///
    inline R1CoVectorTC<Coord>& operator += ( const R1CoVectorTC<Coord>& );
    ///
    inline R1CoVectorTC<Coord>& operator -= ( const R1CoVectorTC<Coord>& );
    ///
    inline R1CoVectorTC<Coord>& operator *= ( const Coord& );
    ///
    inline R1CoVectorTC<Coord>& operator /= ( const Coord& );
    //@}

    /**@name Unary operators */
    //@{
    ///
    inline R1CoVectorTC<Coord>  operator + () const;
    ///
    inline R1CoVectorTC<Coord>  operator - () const;
    //@}


    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space]; if true, replaces space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;

    //@}

private:
    ///
    Coord d;
};


/** 
\class R1PointTC Rn_Defs.H utils/Rn_Defs.H
\brief A 1D Point
\ingroup SimpleGeometryElems
  

  Essentially a single number. Included for completeness'
  sake. Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template<class Coord>
class R1PointTC
{
public:
    ///
    R1PointTC() { }
    /// Initialize with value
    R1PointTC(Coord d);
    /// Return dimension of the type
    int Dim() const { return 1; }

    /**@name Access grants */
    //@{
    /// Does a check that int is in bounds, i.e., 0 
    inline Coord& operator[] ( const int DECLAREVAR(i) ) { ASSERT( i == 0 ); return u; }
    /// Does a check that int is in bounds, i.e., 0 
    inline const Coord& operator[] ( const int DECLAREVAR(i) ) const { ASSERT(i == 0); return u; }
    //@}

    /**@name Assignment operators */
    //@{
    ///
    inline R1PointTC<Coord>& operator += ( const R1VectorTC<Coord>& );
    ///
    inline R1PointTC<Coord>& operator -= ( const R1VectorTC<Coord>& );
    //@}

    /**@name Point dominance */
    //@{
    /// Returns true if value is less than value in in_pt
    inline WINbool operator < ( const R1PointTC<Coord>&in_pt ) const;
    /// Returns true if value is less than or equal to value in in_pt
    inline WINbool operator <= ( const R1PointTC<Coord>&in_pt ) const;
    //@}


    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const;
    ///
    inline WINbool Read(ifstream &in) ;
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is number [space]; if true, replaces space with a carriage return */
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const;
    //@}

private:
    ///
    Coord u;
};

/** 
\brief A 1D Ray (point plus vector)
\class R1RayTC Rn_Defs.H utils/Rn_Defs.H
\ingroup SimpleGeometryElems


  Essentially a single number. Included for completeness'
  sake. Type-safe addition/subtraction/multiplication for points,
  vectors, and co-vectors. See operators in documentation for
  non-member functions.

 */
template< class Coord >
class R1RayTC
{
private:
    ///
    R1PointTC<Coord> pOrg;
    ///
    R1VectorTC<Coord> vDir;

public:
    /**@name Access grants */
    //@{
    ///
    const R1PointTC<Coord> &Pt() const { return pOrg; }
    ///
    const R1VectorTC<Coord> &Dir() const { return vDir; }
    //@}

    /**@name Assignment operators */
    //@{
    ///
    R1PointTC<Coord> &Pt() { return pOrg; }
    ///
    R1VectorTC<Coord> &Dir() { return vDir; }
    //@}

    /**@name Constructor/destructor */
    //@{
    ///
    R1RayTC() { }
    ///
    R1RayTC( const R1PointTC<Coord> & _pOrg, const R1VectorTC<Coord>& _vDir )
        { pOrg = _pOrg; vDir = _vDir; }
    ///
    ~R1RayTC() {}
    //@}

    /**@name Read/write/print */
    //@{
    ///
    inline void    Write(ofstream &out) const { pOrg.Write(out); vDir.Write(out); }
    ///
    inline WINbool Read(ifstream &in) { pOrg.Read(in); return vDir.Read(in); }
    /** \brief Prints to standard out.
     * @param in_bDoReturn Default is point [space] ray [space]; if true, replaces spaces with carriage returns */
    ///
    inline void    Print( const WINbool in_bDoReturn = FALSE ) const { pOrg.Print( in_bDoReturn ); vDir.Print( in_bDoReturn ); }
    //@}

};



#include "Rn_Vector1_i.H"
#include "Rn_Point1_i.H"
#include "Rn_CoVector1_i.H"
#include "Rn_Unary1_i.H"
#include "Rn_Binary1_i.H"
#include "Rn_Io1_i.H"

/** \ingroup SimpleGeometryElems Typedefs */
//@{
/// From R1PointTC
typedef R1PointTC<double> R1Pt;
/// From R1PointTC
typedef R1PointTC<float>  R1Pt_f;
/// From R1PointTC
typedef R1PointTC<int>    R1Pt_i;

/// From R1VectorTC
typedef R1VectorTC<double> R1Vec;
/// From R1VectorTC
typedef R1VectorTC<float>  R1Vec_f;
/// From R1VectorTC
typedef R1VectorTC<int>    R1Vec_i;

/// From R1CoVectorTC
typedef R1CoVectorTC<double> R1CoVec;
/// From R1RayTC
typedef R1RayTC<double> R1Ray;

//@}

#endif
// #ifndef __BASETYPE_AFFINE1_H


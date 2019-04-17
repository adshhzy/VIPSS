#ifndef _C2_LFTTC_H_DEFS_
#define _C2_LFTTC_H_DEFS_

/** Linear fractional transform
  * Essentially a 2x2 matrix made of complex numbers
  */
template <class T>
class C2LFTTC {
protected:
    C2PointTC<T> m_mat[2][2];

public:
    /**@name Accessors */
    //@{
    /// Matrix entry
    const C2PointTC<T> &operator()( const int i, const int j ) const { return m_mat[i][j]; }
    //@}

    /**@name Queries */
    //@{
    /// 
    WINbool IsIdentity( const T dEps = RNEpsilon_d ) const;
    //@}

    //@name Assignment */
    //@{
    /// Matrix entry
    C2PointTC<T> &operator()( const int i, const int j ) { return m_mat[i][j]; }
                   
    /// Identity
    ///  1 0
    ///  0 1
    void SetIdentity(  );
    /// Scale the point up (s is the scale factor)
    ///  s 0
    ///  0 1
    void SetScale( const C2PointTC<T> &in_pt );
    /// translate 
    ///  1 t
    ///  0 1
    void SetTrans( const C2PointTC<T> &in_pt );
    /// Invert matrix
    ///   a b ^-1 = d -b
    //    c d       -c a
    void SetInvert( const C2LFTTC<T> &in_mat );
    /// Rotation
    ///   1  r
    ///   r* 1
    void SetFlip( const C2PointTC<T> &in_pt );
    //@}

    //@name operators */
    //@{
    /// Multiplication by a point (like multiplication of a homogeneous point by a matrix)
    ///   a b   z  = (a z + b) / (c z + d)
    ///   c d   1
    C2PointTC<T> operator *( const C2PointTC<T> &in_pt ) const;

    /// Multiplication by a matrix (like a regular matrix)
    C2LFTTC<T> operator *( const C2LFTTC<T> &in_mat ) const;
    //@}

    /**@name Comparisonn operators */
    //@{
    /// Comparison
    WINbool operator==( const C2LFTTC<T> &in_mat ) const;
    ///
    WINbool operator!=( const C2LFTTC<T> &in_mat ) const;
    //@}


    /**@name Constructors */
    //@{
    C2LFTTC( const C2PointTC<T> &in_a00,
             const C2PointTC<T> &in_a01,
             const C2PointTC<T> &in_a10,
             const C2PointTC<T> &in_a11 );
    ///
    C2LFTTC(  ) { SetIdentity(); }
    ///
    C2LFTTC( const C2LFTTC<T> &in_mat ) { *this = in_mat; }
    ~C2LFTTC() { }
    //@}

    ///
    void Write(ofstream &out) const ;
    ///
    void Read(ifstream &in) ;
    ///
    void WriteBinary(ofstream &out) const ;
    ///
    void ReadBinary(ifstream &in) ;

    static void TestLFT();
};

#include "C2_LFTTC_i.H"

typedef C2LFTTC<double>   C2LFT;

#endif

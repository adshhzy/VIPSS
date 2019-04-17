#ifndef _C2_POINT_TC_H_DEFS_
#define _C2_POINT_TC_H_DEFS_

#include <utils/Rn_Defs.H>

/**name
 * A complex number.
 * Mirrors the complex number class that windows has yet to deign to add.
 */
template <class T>
class C2PointTC {
protected :
    T m_dRe, m_dIm;

public:
    /**@name Access */
    //@{
    ///
    T real() const { return m_dRe; }
    ///
    T imag() const { return m_dIm; }
    ///
    T Length() const { return sqrt( m_dRe * m_dRe + m_dIm * m_dIm ); }
    ///
    WINbool IsReal( const T in_dEps = RNEpsilon_d) const { return fabs( m_dIm ) < in_dEps ? TRUE : FALSE; }
    //@}

    /**@name Setting */
    //@{
    ///
    T &real() { return m_dRe; }
    ///
    T &imag() { return m_dIm; }
    ///
    void Set( const T in_dRe, const T in_dIm );
    ///
    void SetRadiusTheta( const T in_dRadius, const T in_dTheta );
    //@}

    /**@name Assignment operators */
    //@{
    /// 
    C2PointTC<T>& operator += ( const C2PointTC<T> & );
    ///
    C2PointTC<T>& operator -= ( const C2PointTC<T> & );
    ///
    C2PointTC<T>& operator *= ( const C2PointTC<T> & );
    ///
    C2PointTC<T>& operator /= ( const C2PointTC<T> & );
    //@}

    /**@name operators */
    //@{
    /// 
    C2PointTC<T> operator+( const C2PointTC<T> &in_pt ) const;
    /// 
    C2PointTC<T> operator+( const T &in_pt ) const;
    ///
    C2PointTC<T> operator-( const C2PointTC<T> &in_pt ) const;
    /// 
    C2PointTC<T> operator-( const T &in_pt ) const;
    ///
    C2PointTC<T> operator-( ) const;
    ///
    C2PointTC<T> operator*( const C2PointTC<T> &in_pt ) const;
    ///
    C2PointTC<T> operator/( const C2PointTC<T> &in_pt ) const;
    //@}

    /**@name Comparison operators */
    //@{
    ///
    WINbool operator==( const C2PointTC<T> &in_pt ) const;
    ///
    WINbool operator!=( const C2PointTC<T> &in_pt ) const;
    //@}

    ///
    C2PointTC( const T in_dRe, const T in_dIm = 0.0 ) { Set(in_dRe, in_dIm); }
    ///
    C2PointTC( const C2PointTC<T> &in_pt ) { Set(in_pt.real(), in_pt.imag()); }
    ///
    C2PointTC() : m_dRe(0), m_dIm(0) {}
    ///
    ~C2PointTC() {}

    ///
    void Write(ostream &out) const { out << m_dRe << " " << m_dIm << "\n"; }
    ///
    void Read(ifstream &in) { in >> m_dRe >> m_dIm; }

    static WINbool Test();
};

#include "C2_PointTC_i.H"

typedef C2PointTC<double> C2Pt;

#endif

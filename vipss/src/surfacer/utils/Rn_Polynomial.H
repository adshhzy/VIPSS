#ifndef _RN_DEFS_POLYNOMIAL_H
#define _RN_DEFS_POLYNOMIAL_H

#include <utils/Rn_Defs.H>
#include <utils/Mesh_Array.H>

/** \class RNPolynomial Rn_Polynomial.H utils/Rn_Polynomial.H
 \ingroup Polynomials
 \brief  A polynomial in one or two variables.

 Defined by the dimension (number of variables) and the coefficients (degree + 1). <br>
 Files: <br>
  - include/utils/Rn_Polynomial.H <br>
  - include/utils/R2_Polynomial_i.H <br>
  - src/utils/geometry/RNPolynomial.cpp <br>
  */
class RNPolynomial : private Array<double> {
private:
    /// number or variables (x,y, etc)
    int m_iDim;
    /// Degree of polynomial
    int m_iDegree;
    
    void Compress(); // remove excess degrees

public:
    /**@name Access to data */
    //@{
    /// Number of variables.
    inline int     Dim()            const { return m_iDim; }
    /// Degree of polynomial +1 raised to dimension.
    inline int     Num_coefs()      const { return num(); }
    /// Highest power of polynomial
    inline int     Degree()         const { return m_iDegree; }
    /// Assumes polynomial in one variable
    inline double  operator[](int i) const ;
    /// Assumes polynomial in one variable
    inline double &operator[](int i) ;

    /// Assumes polynomial in one variable
    inline double  Coef(int which) const;
    /// Assumes polynomial in one variable
    inline double &Coef(int which) ;
    /// Assumes polynomial in two variables
    inline double  Coef(int which_s, int which_t) const;
    /// Assumes polynomial in two variables
    inline double &Coef(int which_s, int which_t) ;
    //@}

    /**@name Evaluations on polynomial */
    //@{
    /// Evaluate one variable polynomial
    double  operator()(double in_t) const;
    /// Evaluate two variable polynomial
    double  operator()(double in_s, double in_t) const;
    //@}

    /**@name Operations on polynomials
       Polynomials must have same dimension */
    //@{
    /// 
    RNPolynomial operator+(const RNPolynomial &) const ;
    ///
    RNPolynomial operator-(const RNPolynomial &) const ;
    ///
    RNPolynomial operator-() const ;
    ///
    RNPolynomial operator*(const RNPolynomial &) const ;
    ///
    RNPolynomial Divide(const RNPolynomial &, RNPolynomial &out_polyRemainder ) const ;
    /// 
    RNPolynomial Integrate( double in_dConstant ) const ;
    /// Integrate from t0 to t1
    double       Integrate( double t0, double t1 ) const ;
    ///
    RNPolynomial Differentiate() const ;
    /// Composition
    RNPolynomial operator() ( const RNPolynomial & ) const ;
    ///
    RNPolynomial operator*  ( double in_dScale ) const ;
    //@}

    /**@name Comparison routine */
    //@{
    ///
    WINbool       operator==(const RNPolynomial & ) const ;
    //@}

    /**@name Constructors and assignments */
    //@{
    ///
    RNPolynomial &operator=( const RNPolynomial & ) ;
    ///
    RNPolynomial( const RNPolynomial & );
    ///
    RNPolynomial( const Array<double> &in_adCoefs );
    /// Polynomial with in_iDim variables, highest power in_iDegree
    RNPolynomial(int in_iDim, int in_iDegree);
    /// Polynomial of dimension 1, degree in_iNumCoefs - 1
    RNPolynomial(int in_iNumCoefs);
    /// Default polynomial, evaluates to 0
    RNPolynomial();
    //@}

    ///
    virtual ~RNPolynomial();

    /**@name Read, write, print */
    //@{
    /// Debugging print.
    void    Print() const ;
    /// Writes dimension, degree, and coefficients
    void    Write ( ofstream &out ) const ;
    /// Reads dimension, degree, and coefficients
    WINbool Read  ( ifstream &in );
    //@}

    ///
    static WINbool Test();
};

#include "Rn_Polynomial_i.H"
#endif

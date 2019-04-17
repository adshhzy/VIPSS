/* Copyright 1995, Brown Computer Graphics Group.  All Rights Reserved. */


#ifndef CRV_BASIS_DEFS
#define CRV_BASIS_DEFS

#include <utils/Rn_Polynomial.H>

/* -------------------  Constants  ------------------------------- */

/* ---------------------  Classes  ------------------------------- */


/** \class CRVSpace Crv_Basis.H utils/Crv_Basis.H
   \ingroup SplinesBSplines

   \brief The order (C^k) of a curve 

  This is just used to distinguish between the number of control points (int) and the degree of the curve (k)

<br> Files:
- include/utils/Crv_Basis.H
- src/utils/splines/CRVSpace.cpp
- src/utils/splines/CRVBasis.cpp
*/
class CRVSpace {
protected:
    int m_iOrder;
    
public:
    ///
    inline int Order() const { return m_iOrder; }
    
    ///
    inline CRVSpace(const CRVSpace &in_crvSpace) 
        : m_iOrder(in_crvSpace.Order()) { ASSERT(m_iOrder > -2); }
    ///
    inline CRVSpace(int in_iOrder) : m_iOrder(in_iOrder) { ASSERT(m_iOrder > -2); }
    ///
    virtual ~CRVSpace() { }

    ///
    CRVSpace(ifstream &in, WINbool bIs_binary = FALSE);
};


/** \class CRVBasis Crv_Basis.H utils/Crv_Basis.H
   \ingroup SplinesBSplinesKnot

   \brief A single basis funtion for a spline curve.

    A single non-uniform basis function. <br>
    Has an order+1 knot vector, going from Start_t to End_t. <br>
    Normally, these are created by the curve knot classes and don't
    need to be created by hand. <br>

    Coef of polynomials. coefs[j][i] is for x^i in range 
    knot[j] to knot[j+1]
    Recalculated when knot vector changes
*/
class CRVBasis : public CRVSpace {
private:
    /// The part of the knot vector that this is non-zero over (order k + 1 knots)
    Array<double>              m_adKnots;
    /// The actual polynomials of the individual segments (k segments)
    Array<RNPolynomial>        m_aPolys;

    inline RNPolynomial       &Poly(const int wh) { return m_aPolys[wh]; }
    inline const RNPolynomial &Poly(const int wh) const { return m_aPolys[wh]; }
    
    /// Fill in the coefs vec (recursive call)
    void Calc_coefs(const CRVBasis &ith, const CRVBasis &ithPlus1);
    
    static int m_iDefaultConstructOrder;
    
public:
    /// Use to set the default order when constructing, if not specified. Useful for Array constructions.
    static void SetDefaultOrder(int in_iOrder) {m_iDefaultConstructOrder = in_iOrder;}
    
    /**@name Access */
    //@{
    /// Order of curve + 1
    inline int Num_knots() const { return m_adKnots.num(); }
    ///
    inline const Array<double> &Knots() const { return m_adKnots; }
    ///
    inline double Knot(const int in_i) const { return m_adKnots[in_i]; }
    /// Beginning of knot vector
    inline double Start_t() const { return m_adKnots[0]; }
    /// End of knot vector
    inline double End_t() const { return m_adKnots.last(); }
    //@}

    /**@name Evaluation */
    //@{
    /// Evaluate the polynomial at t (goes to zero outside of Start_t() End_t()).
    double operator()(const double in_dT) const;
    /// Return the polynomial for the given segment
    inline const RNPolynomial &Segment(const int in_i) const { return Poly(in_i); }
    ///
    double Derivative(const double in_dT) const;
    /// If order < 2, returns 0
    double SecondDerivative(const double in_dT) const;
    /// Eval, and calculate the first and second derivative
    void AllEval( const double in_dT, double &out_dEval, double &out_dDeriv, double &out_dSecondDeriv ) const;
    /// Eval, and calculate the first derivative
    void AllEval( const double in_dT, double &out_dEval, double &out_dDeriv ) const;
    ///
    double Integrate(const double in_dT0, const double in_dT1) const;
    /// Make a basis function which is the derivative of this one.
    CRVBasis Derivative() const;
    //@}

    /**@name Never, ever call this by hand */
    //@{
    /** Change the knot vector. 
        The knot vector is taken from the vector in_adKnots by taking the
         ith - ith+order+1 elements.  (This makes it easy for a curve
         to create all of it's basis functions.)  
         @param in_adKnots The entire knot vector.
         @param in_iStart The knot where this basis function starts.
    */
    void Change_knot(const Array<double> &in_adKnots, const int in_iStart = 0);
    //@}

    /**@name Constructors and assignments */
    //@{
    /// Copy basis
    CRVBasis &operator=(const CRVBasis &);
    /// 
    WINbool operator==(const CRVBasis &) const ;
    ///
    CRVBasis(const CRVBasis &);
    
    /** Create a basis with the given knot vector.
        (which must be of dimension order+1) */
    CRVBasis(const int in_iOrder, 
             const Array<double> &in_adKnot, int in_iStartArrayIndex);
    
    /** Evenly spaced basis.
        Create a basis whose t values fall between left_end and right_end and are evenly spaced */
    CRVBasis( const int in_iOrder, const double in_dLeftEnd,  const double in_dRightEnd);
    
    /// Create a basis with the order in defaultConstructOrder
    CRVBasis();
    
    ///
    virtual ~CRVBasis() {}
    //@}
   
    ///
    virtual void Print() const;
    ///
    void Write( ofstream & ) const;
    ///
    void Read( ifstream & ) ;
    ///
    void WriteBinary( ofstream & ) const;
    ///
    void ReadBinary( ifstream & ) ;
};

#endif


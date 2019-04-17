/* -------------------------------------------------------------------------
 * DESCR   :	Compare two points
 * ------------------------------------------------------------------------- */
inline
WINbool R2Ellipse::operator==(const R2Ellipse &in_Ellipse) const 
{
    return ApproxEqual( *this, in_Ellipse );
}

inline
WINbool R2Ellipse::ApproxEqual(const R2Ellipse &in_ell1, const R2Ellipse &in_ell2, const double in_dEps)  
{
    if (::ApproxEqual(in_ell1.Center(), in_ell2.Center(), in_dEps) == FALSE) return FALSE;
    
    if ( RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation(), in_dEps ) || 
         RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation() - M_PI, in_dEps ) ||
         RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation() + M_PI, in_dEps ) ) {

        if (!RNIsZero(in_ell2.XRadius() - in_ell1.XRadius(), in_dEps)) return FALSE;    
        if (!RNIsZero(in_ell2.YRadius() - in_ell1.YRadius(), in_dEps)) return FALSE;

        return TRUE;
    }

    if ( RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation() - M_PI / 2.0, in_dEps ) || 
         RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation() + M_PI / 2.0, in_dEps ) ||
         RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation() + 3.0 * M_PI / 2.0, in_dEps ) ||
         RNApproxEqual( in_ell1.Rotation(), in_ell2.Rotation() - 3.0 * M_PI / 2.0, in_dEps ) ) {

        if (!RNIsZero(in_ell2.XRadius() - in_ell1.YRadius(), in_dEps)) return FALSE;    
        if (!RNIsZero(in_ell2.YRadius() - in_ell1.XRadius(), in_dEps)) return FALSE;

        return TRUE;
    }

    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Evaluate implicit equation
 * ------------------------------------------------------------------------- */
inline
double  R2Ellipse::operator()( const R2Pt &in_pt ) const
{
    return m_adCoef[0] * in_pt[0] * in_pt[0] +
           m_adCoef[1] * in_pt[0] * in_pt[1] + 
           m_adCoef[2] * in_pt[1] * in_pt[1] +
           m_adCoef[3] * in_pt[0] +
           m_adCoef[4] * in_pt[1] +
           m_adCoef[5];

}


/* -------------------------------------------------------------------------
 * DESCR   :	Are we inside the circle (open disk)?
 * ------------------------------------------------------------------------- */
inline
WINbool R2Ellipse::Inside(const R2Pt &in_pt) const
{
    if ( (*this)(in_pt) <= 0 ) return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we on the circle ?
 * ------------------------------------------------------------------------- */
inline
WINbool R2Ellipse::On(const R2Pt &in_pt) const
{
    return RNIsZero( (*this)(in_pt) );
}





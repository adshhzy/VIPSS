/* -------------------------------------------------------------------------
 * DESCR   :	Compare two points
 * ------------------------------------------------------------------------- */
inline
WINbool R3Ellipse::operator==(const R3Ellipse &in_Ellipse) const 
{
    if (ApproxEqual(Center(), in_Ellipse.Center()) == FALSE) return FALSE;
    
    if (!RNIsZero(in_Ellipse.m_dA - m_dA )) return FALSE;
    
    if (!RNIsZero(in_Ellipse.m_dB - m_dB )) return FALSE;
    
    if (!RNIsZero(in_Ellipse.m_dC - m_dC )) return FALSE;
    
    if (!RNIsZero(in_Ellipse.m_dZRot - m_dZRot )) return FALSE;
    if (!RNIsZero(in_Ellipse.m_dYRot - m_dYRot )) return FALSE;
    if (!RNIsZero(in_Ellipse.m_dXRot - m_dXRot )) return FALSE;
    
    return TRUE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Evaluate implicit equation
 * ------------------------------------------------------------------------- */
inline
double  R3Ellipse::operator()( const R3Pt &in_pt ) const
{
    const R4Matrix mat = ToNormal();

    const R3Pt ptNormal = mat * in_pt;

    return (ptNormal[0] * ptNormal[0]) / (m_dA * m_dA) + 
           (ptNormal[1] * ptNormal[1]) / (m_dB * m_dB) +
           (ptNormal[2] * ptNormal[2]) / (m_dC * m_dC) - 1.0;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we inside the circle (open disk)?
 * ------------------------------------------------------------------------- */
inline
WINbool R3Ellipse::Inside(const R3Pt &in_pt) const
{
    if ( (*this)( in_pt ) <= 0 ) return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we on the circle ?
 * ------------------------------------------------------------------------- */
inline
WINbool R3Ellipse::On(const R3Pt &in_pt) const
{
    if (RNIsZero( (*this)(in_pt) ) ) return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Set two points to be equal
 * ------------------------------------------------------------------------- */
inline
R3Ellipse &
R3Ellipse::operator=(const R3Ellipse &in_s)
{
    m_ptCenter = in_s.Center();
    m_dA = in_s.m_dA;
    m_dB = in_s.m_dB;
    m_dC = in_s.m_dC;
    m_dXRot = in_s.m_dXRot;
    m_dYRot = in_s.m_dYRot;
    m_dZRot = in_s.m_dZRot;
    
    return *this;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Make a Ellipse from a Ellipse
 * ------------------------------------------------------------------------- */
inline
R3Ellipse::R3Ellipse(const R3Pt &in_ptCenter)
    : m_ptCenter(in_ptCenter)
{
    m_dA = m_dB = m_dC = 1.0;
    m_dXRot = m_dYRot = m_dZRot = 0.0;
}

inline
R3Ellipse::R3Ellipse() 
{
    RNZero(m_ptCenter);
    m_dA = m_dB = m_dC = 1.0;
    m_dXRot = m_dYRot = m_dZRot = 0.0;
}



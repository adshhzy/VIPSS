/* -------------------------------------------------------------------------
 * DESCR   :	Compare two points
 * ------------------------------------------------------------------------- */
inline
WINbool R2Sphere::operator==(const R2Sphere &in_sphere) const 
{
    if (ApproxEqual(Center(), in_sphere.Center()) == FALSE) return FALSE;
    
    if (!RNIsZero(in_sphere.Radius() - Radius())) return FALSE;
    
    return TRUE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we inside the circle (open disk)?
 * ------------------------------------------------------------------------- */
inline
WINbool R2Sphere::Inside(const R2Pt &in_pt) const
{
    if (Length(in_pt - Center()) < Radius()) return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we inside the circle (open disk)?
 * ------------------------------------------------------------------------- */
inline
WINbool R2Sphere::OnOrInside(const R2Pt &in_pt) const
{
    const double dLen = Length(in_pt - Center());
    if ( dLen <= Radius() ) return TRUE;
    if ( RNApproxEqual( dLen, Radius() ) )
        return TRUE;
    
    return FALSE;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Are we on the circle ?
 * ------------------------------------------------------------------------- */
inline
WINbool R2Sphere::On(const R2Pt &in_pt) const
{
    if (RNIsZero(Radius() - Length(in_pt - Center()))) return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Distance between centers
 * ------------------------------------------------------------------------- */
inline
WINbool R2Sphere::Overlaps( const R2Sphere & in_circ ) const
{
    const double dLen = Length( in_circ.Center() - Center() );
    if ( dLen <= in_circ.Radius() + Radius() )
        return TRUE;

    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Set two points to be equal
 * ------------------------------------------------------------------------- */
inline
R2Sphere &
R2Sphere::operator=(const R2Sphere &in_s)
{
    m_ptCenter = in_s.Center();
    m_dRadius = in_s.Radius();
    
    return *this;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Make a sphere from a sphere
 * ------------------------------------------------------------------------- */
inline
R2Sphere::R2Sphere(const R2Pt &in_ptCenter, double in_dR)
    : m_dRadius(in_dR), m_ptCenter(in_ptCenter)
{
}

inline
R2Sphere::R2Sphere(double in_dR) : m_dRadius(in_dR)
{
    RNZero(m_ptCenter);
}



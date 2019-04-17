/* -------------------------------------------------------------------------
 * DESCR   :	Interior area
 * ------------------------------------------------------------------------- */
inline double R3Sphere::Area() const  
{ 
    return 4.0 / 3.0 * M_PI * Radius() * Radius() * Radius(); 
}

/* -------------------------------------------------------------------------
 * DESCR   :	Compare two points
 * ------------------------------------------------------------------------- */
inline
WINbool R3Sphere::operator==(const R3Sphere &in_sphere) const 
{
    if (ApproxEqual(Center(), in_sphere.Center()) == FALSE) return FALSE;
    
    if (!RNIsZero(in_sphere.Radius() - Radius())) return FALSE;
    
    return TRUE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we inside the circle (open disk)?
 * ------------------------------------------------------------------------- */
inline
WINbool R3Sphere::Inside(const R3Pt &in_pt) const
{
    if (Length(in_pt - Center()) < Radius()) return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Are we on the circle ?
 * ------------------------------------------------------------------------- */
inline
WINbool R3Sphere::On(const R3Pt &in_pt) const
{
    const double dLen = Length(in_pt - Center());
    if (RNIsZero(Radius() - dLen)) 
        return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Set two points to be equal
 * ------------------------------------------------------------------------- */
inline
R3Sphere &
R3Sphere::operator=(const R3Sphere &in_s)
{
    m_ptCenter = in_s.Center();
    m_dRadius = in_s.Radius();
    
    return *this;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Make a sphere from a sphere
 * ------------------------------------------------------------------------- */
inline
R3Sphere::R3Sphere(const R3Pt &in_ptCenter, double in_dR)
    : m_dRadius(in_dR), m_ptCenter(in_ptCenter)
{
}

inline
R3Sphere::R3Sphere(double in_dR) : m_dRadius(in_dR)
{
    RNZero(m_ptCenter);
}



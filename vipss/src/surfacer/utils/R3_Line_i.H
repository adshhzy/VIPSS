/* -------------------------------------------------------------------------
 * DESCR   :	Variations on FindPtOnLine
 * ------------------------------------------------------------------------- */
inline double R3Line::Dist_to_line(const R3Pt &in_pt) const
{
   R3Pt out_pt;
   double out_t, out_d;
   FindPtOnLine(in_pt, out_pt, out_t, out_d);

   return out_d;
}

inline R3Pt R3Line::Project(const R3Pt &in_pt) const
{
    R3Pt out_pt;
    double out_t, out_d;
    FindPtOnLine(in_pt, out_pt, out_t, out_d);

    return out_pt;
}

inline
double R3Line::Projected_dist_on_line(const R3Pt &in_pt) const
{
    R3Pt out_pt;
    double out_t, out_d;
    FindPtOnLine(in_pt, out_pt, out_t, out_d);

    return out_t;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Determine if the two lines are parallel
 * ------------------------------------------------------------------------- */
inline
WINbool IsParallel(const R3Line &in_l1, const R3Line &in_l3) 
{
    double res = 1.0f - (double) fabs( Dot(in_l1.Vec(), in_l3.Vec()) );
    if ( RNIsZero( res ) )
        return TRUE;

    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Determine if the two lines are perpendicular
 * ------------------------------------------------------------------------- */
inline
WINbool IsPerpendicular(const R3Line &in_l1, const R3Line &in_l3)
{
    if ( RNIsZero( Dot( in_l1.Vec(), in_l3.Vec() ) ) )
        return TRUE;
    
    return FALSE;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Comparison
 * ------------------------------------------------------------------------- */
inline WINbool R3Line::operator==(const R3Line &in_l) const 
{
    if (ApproxEqual(Vec(), in_l.Vec()) == FALSE) return FALSE;
    
    if (ApproxEqual(Pt(), in_l.Pt()) == TRUE) return TRUE;
    
    return FALSE;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Constructors
 * ------------------------------------------------------------------------- */
inline R3Line &R3Line::operator=(const R3Line &in_l)
{
    m_pt = in_l.m_pt;
    m_vec = in_l.m_vec;

    return *this;
}

inline void R3Line::SetVec(const R3Vec &in_vec) 
{ 
    //ASSERT( ! RNIsZero( Length( m_vec ) ) );
    m_vec = UnitSafe(in_vec); 
}

inline R3Line::R3Line(const R3Pt &in_pt, const R3Vec &in_vec) 
    : m_pt(in_pt), m_vec( UnitSafe( in_vec ) )
{
    //ASSERT( ! RNIsZero( Length( m_vec ) ) );
}

inline R3Line::R3Line(const R3Pt &in_pt1, const R3Pt  &in_pt2) 
    : m_pt( in_pt1 ), m_vec( UnitSafe( in_pt2 - in_pt1 ) )
{
    //ASSERT( ! RNIsZero( Length( m_vec ) ) );
}

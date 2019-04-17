inline
WINbool R3Line_seg::IsPtOnSeg(const R3Pt &in_pt) const
{
    R3Pt out_pt;
    double out_t, out_d;
    return FindPtOnSeg(in_pt, out_pt, out_t, out_d);
}

inline
double R3Line_seg::Dist_to_seg(const R3Pt &in_pt) const
{
    R3Pt out_pt;
    double out_t, out_d;
    FindPtOnSeg(in_pt, out_pt, out_t, out_d);

    return out_d;
}

inline
R3Pt R3Line_seg::Project(const R3Pt &in_pt) const
{
    R3Pt out_pt;
    double out_t, out_d;
    FindPtOnSeg(in_pt, out_pt, out_t, out_d);

    return out_pt;
}

inline
double 
R3Line_seg::Projected_dist_on_seg(const R3Pt &in_pt) const
{
    R3Pt out_pt;
    double out_t, out_d;
    FindPtOnSeg(in_pt, out_pt, out_t, out_d);

    return out_t;
}


/* -------------------------------------------------------------------------
 * DESCR   :	Parallel & perpendicular - use RNline<double> routines
 * ------------------------------------------------------------------------- */
inline
WINbool
IsParallelSeg(const R3Line_seg &in_l1, const R3Line_seg &in_l2)
{
    R3Line l1(in_l1.P1(), in_l1.P2() - in_l1.P1());
    R3Line l2(in_l2.P1(), in_l2.P2() - in_l2.P1());
    
    return IsParallel(l1, l2);
}

inline
WINbool
IsPerpendicularSeg(const R3Line_seg &in_l1, const R3Line_seg &in_l2) 
{
    R3Line l1(in_l1.P1(), in_l1.P2() - in_l1.P1());
    R3Line l2(in_l2.P1(), in_l2.P2() - in_l2.P1());
    
    return IsPerpendicular(l1, l2);
}



/* -------------------------------------------------------------------------
 * DESCR   :	Compare two lines
 * ------------------------------------------------------------------------- */
inline
WINbool R3Line_seg::operator==(const R3Line_seg &in_l) const 
{
    if (P1() == in_l.P1() && P2() == in_l.P2())
        return TRUE;
    
    if (P1() == in_l.P2() && P2() == in_l.P1())
        return TRUE;
    
    return FALSE;
}

/* -------------------------------------------------------------------------
 * DESCR   :	Calculate the length of the line segment
 * ------------------------------------------------------------------------- */
inline
double R3Line_seg::Length() const
{
    return ::Length(P2() - P1() );
}


inline
R3Pt R3Line_seg::operator()( double s) const 
{ 
    return m_pt1 + (m_pt2 - m_pt1) * s; 
}

inline
R3Line_seg  &R3Line_seg::operator=(const R3Line_seg &in_l)
{
    m_pt1 = in_l.m_pt1;
    m_pt2 = in_l.m_pt2;

    return *this;
}

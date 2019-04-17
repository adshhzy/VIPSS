/*
 *  Utils_EdgeFilter.H
 *  utils
 *
 *  Created by Cindy Grimm on 12/16/11.
 *  Copyright 2011 Washington University in St. Louis. All rights reserved.
 *
 */

#ifndef _UTILS_EDGE_FILTER_DEFS_H_
#define _UTILS_EDGE_FILTER_DEFS_H_

#include <WINSystemDefines.H>

/* Edge filter function
 * Difference of Gaussians in the x direction (e1(x) - e2(x)) multiplied
 * by a Gaussian in the y direction that is designed to fall off at roughly 3 times 
 * the width of the x direction filter
 * Ie, if you clip the filter at 2 times sigma, the filter will be mostly zero by then, with 1/2 of the
 * radius a "big" peak
 *
 * sigma 1 is the 
 */
class Utils_EdgeFilter {
private:
    double m_dTheta, m_dSigma1, m_dSigma2, m_dSigmaY;
    double m_dRotCos, m_dRotSin;
    double m_dScaleFactorFirst, m_dScaleFactorSecond, m_dScaleFactorY; 
    double m_dGaussScaleFactorFirst, m_dGaussScaleFactorSecond, m_dGaussScaleFactorY; 
    
public:
    /** Set the filter up
     * @param in_dTheta Angle of edge
     * @param in_dSigma Width of wider Gaussian
     * @param in_dSigmaRatio (should be < 1) width of skinnier Gaussian. Default is 0.25
     */
    void SetFilterParams( const double in_dTheta, const double in_dSigma, const double in_dSigmaRatio )
    {
        m_dTheta = in_dTheta;
        m_dRotCos = cos(m_dTheta);
        m_dRotSin = sin(m_dTheta);
        m_dSigma1 = in_dSigma < 0.0 ? 1.0 : in_dSigma;
        m_dSigma2 = m_dSigma1 * ( in_dSigmaRatio < 0.0 ? 0.25 : in_dSigmaRatio);
        m_dSigmaY = 3.0 * m_dSigma1;
        m_dScaleFactorFirst = 1.0 / ( m_dSigma1 * std::sqrt(2.0*M_PI) );
        m_dScaleFactorSecond = 1.0 / ( m_dSigma2 * std::sqrt(2.0*M_PI) );
        m_dScaleFactorY = 1.0 / ( m_dSigmaY * std::sqrt(2.0*M_PI) );
        m_dGaussScaleFactorFirst = (2.0*m_dSigma1*m_dSigma1);
        m_dGaussScaleFactorSecond = (2.0*m_dSigma2*m_dSigma2);
        m_dGaussScaleFactorY = (2.0*m_dSigmaY*m_dSigmaY);
    }
    
    double operator()( const double in_dX, const double in_dY ) const
    {
        const double dX = m_dRotCos * in_dX + m_dRotSin * in_dY;
        const double dY = -m_dRotSin * in_dX + m_dRotCos * in_dY;
        const double dYFalloff = m_dScaleFactorY * std::exp( -pow(dY,2.0) / m_dGaussScaleFactorY );
        const double dXFalloff = 
        -m_dScaleFactorFirst * std::exp( -pow(dX,2.0) / m_dGaussScaleFactorFirst ) +
        m_dScaleFactorSecond * std::exp( -pow(dX,2.0) / m_dGaussScaleFactorSecond );
        return dXFalloff * dYFalloff;
    }
    double MaxEdgeFunction( ) const
    {
        return m_dScaleFactorY * ( m_dScaleFactorFirst - m_dScaleFactorSecond );
    }
    
    
    Utils_EdgeFilter() { SetFilterParams( 0.0, 1.0, 0.25 ); }
    ~Utils_EdgeFilter() {}
};

#endif
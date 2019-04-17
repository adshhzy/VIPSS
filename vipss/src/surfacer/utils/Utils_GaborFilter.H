/*
 *  Utils_GaborFilter.H
 *  utils
 *
 *  Created by Cindy Grimm on 12/13/11.
 *  Copyright 2011 Washington University in St. Louis. All rights reserved.
 *
 */

#ifndef _UTILS_GABOR_FILTER_DEFS_H_
#define _UTILS_GABOR_FILTER_DEFS_H_

#include <WINSystemDefines.H>

/* Gabor filter function
 * dx : the horizontal distance between the focal point and center of the filter
 * dy : the vertical distance between the focal point and center of the filter
 * angle (theta) : the angle of the Gabor filter
 * frequency (phi) : the frequency of the Gabor filter
 * sigma : the variance of the normal distribution, or: the effective distance (i.e. size) of the filter.
 *
 * g(x,y) exp( 2 pi i phi( x cos theta + y sin theta) )
 * 
 * g(x,y) = 1 / (2 pi sigma^2)  e( -(x^2 + y^2) / (2 sigma^2) )
 *
 * G(x,y,S,F,W,P)=k*Gaussian(x,y,S)*(Sinusoid(x,y,F,W,P)-DC(F,S,P)), 
 * where: 
 * Gaussian(x,y,S)=exp(-pi*S^2*(x^2+y^2)) 
 * Sinusoid(x,y,F,W,P)=exp(j*(2*pi*F*(x*cos(W)+y*sin(W))+P))) 
 * DC(F,S,P)=exp(-pi*(F/S)^2+j*P)
 */
class Utils_GaborFilter {
private:
    double m_dTheta, m_dSigma, m_dFrequency;
    double m_dGaborScaleX, m_dGaborScaleY, m_dScaleFactor, m_dGaussScaleFactor; 
    
public:
    void SetFilterParams( const double in_dTheta, const double in_dSigma, const double in_dFrequency )
    {
        m_dTheta = in_dTheta;
        m_dSigma = in_dSigma < 0.0 ? 1.0 : in_dSigma;
        m_dFrequency = in_dFrequency;
        m_dGaborScaleX = std::cos(in_dTheta)* m_dFrequency; //Sometimes called 'k0x'
        m_dGaborScaleY = std::sin(in_dTheta)* m_dFrequency; //Sometimes called 'k0y'
        m_dScaleFactor = 1.0 / ( std::sqrt(2.0*M_PI*m_dSigma*m_dSigma) );
        m_dGaussScaleFactor = (2.0*m_dSigma*m_dSigma);
    }
    
    double operator()( const double in_dX, const double in_dY ) const
    {
        return (std::cos((m_dGaborScaleX * in_dX)+(m_dGaborScaleY * in_dY))
                *  m_dScaleFactor // g(x,y)            
                * std::exp(-1.0 * ((in_dX*in_dX)+(in_dY*in_dY)) / m_dGaussScaleFactor)
                
                );
    }
    std::pair<double,double> ComplexEval( const double in_dX, const double in_dY ) const
    {
        const double dGaussFactor = m_dScaleFactor * std::exp(-1.0 * ((in_dX*in_dX)+(in_dY*in_dY)) / m_dGaussScaleFactor);
        const double dExpVal = (m_dGaborScaleX * in_dX)+(m_dGaborScaleY * in_dY);
        std::pair<double,double> ptRet;
        ptRet.first = std::cos(dExpVal) * dGaussFactor;
        ptRet.second = std::sin(dExpVal) * dGaussFactor;
        
        return ptRet;
    }
    double MaxGaborFunction( ) const
    {
        return m_dScaleFactor;
    }
    
    
    Utils_GaborFilter() { SetFilterParams( 0.0, 1.0, 1.0 ); }
    ~Utils_GaborFilter() {}
};

#endif
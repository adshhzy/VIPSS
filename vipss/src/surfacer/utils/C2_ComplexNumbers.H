#ifndef _CS_COMPLEX_NUMBERS_DEFS_H
#define _CS_COMPLEX_NUMBERS_DEFS_H

#include <WINSystemDefines.H>
#include <utils/C2_PointTC.H>
#include <utils/C2_LFTTC.H>
#include <utils/Mesh_Array.H>
#include <utils/Rn_Defs.H>
//#include <complex>


///
template <class T>
inline 
void C2FindRoots( const Array<T> &in_adPoly, Array< C2PointTC<T> > &out_adRoots )
{
    int iHighestPower = 0;
    for ( int i = 0; i < in_adPoly.num(); i++ ) {
        if ( !RNIsZero( (T) in_adPoly[i] ) )
            iHighestPower = i;
    }
    
    out_adRoots.need( iHighestPower );

    if ( iHighestPower == 0 ) {
        return;
    } else if ( iHighestPower == 1 ) {
        // a x + b == 0  => x == - b / a
        out_adRoots[0] = C2PointTC<T>( -in_adPoly[0] / in_adPoly[1], 0.0 );
    } else if ( iHighestPower == 2 ) {
        const T dSqrt = pow(in_adPoly[1], 2) - 4.0 * in_adPoly[0] * in_adPoly[2];
        if ( dSqrt >= 0.0 ) {
            out_adRoots[0] = C2PointTC<T>( ( - in_adPoly[1] + sqrt( dSqrt ) ) / (2.0 * in_adPoly[2]), 0.0);
            out_adRoots[1] = C2PointTC<T>( ( - in_adPoly[1] - sqrt( dSqrt ) ) / (2.0 * in_adPoly[2]), 0.0);
        } else {
            out_adRoots[0] = C2PointTC<T>( - in_adPoly[1] / (2.0), -sqrt(-dSqrt) / (2.0 * in_adPoly[2]) );
            out_adRoots[1] = C2PointTC<T>( - in_adPoly[1] / (2.0),  sqrt(-dSqrt) / (2.0 * in_adPoly[2]) );
        }
    } else if ( iHighestPower == 3 ) {
        const T dA = in_adPoly[2] / in_adPoly[3];
        const T dB = in_adPoly[1] / in_adPoly[3];
        const T dC = in_adPoly[0] / in_adPoly[3];
        const T dR = ( 9.0 * dA * dB - 27.0 * dC - 2.0 * pow(dA, 3.0) );
        const T dT = ::pow((T)dA,(int)2) - dB * 3.0;
        const T dC0 = ::pow(2.0, 1.0 / 3.0);
        const T dC1 = ::pow(2.0, 2.0 / 3.0);
        const C2PointTC<T> dC2 = C2PointTC<T>(0,2) * dC0;
        const C2PointTC<T> dQ = pow( dT,3.0) * -4.0;
        const C2PointTC<T> dS = pow( sqrt(dQ + pow(-dR,2.0)) + dR, 1.0 / 3.0);
        out_adRoots[0] = ( ( (C2PointTC<T>) dT * 2.0 * dC0 ) / dS + 
                        dS * dC1 -
                        2.0 * dA ) * (1.0 / 6.0);

        out_adRoots[1] = ( -( dC2 * ( C2PointTC<T>(0,-1) + ::sqrt(3.0) ) * dT ) / dS + 
                           C2PointTC<T>(0,1) * dC1 * (C2PointTC<T>(0,1) + ::sqrt(3.0) ) * dS - 4.0 * dA) / 12.0;   

        out_adRoots[2] = ( ( dC2 * (C2PointTC<T>(0,1) + sqrt(3.0) ) * dT)/
                     dS - (C2PointTC<T>) dC1 * (C2PointTC<T>(0,1) * ::sqrt(3.0) + 1.0)*
                      dS - 4.0 * dA) / 12.0;
    }
}

#endif

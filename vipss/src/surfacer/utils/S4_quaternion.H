#ifndef _IG_QUATERN_H
#define _IG_QUATERN_H

#include <WINSystemDefines.H>
#include <utils/Rn_Defs.H>

/** \class S4QuaternionTC S4_quaternion.H utils/S4_quaternion.H
 \ingroup QuatGeometry
 \brief Quaternion

    An alternative to using Euler angles or a matrix for storing rotations. Supports
    algebraic operations on rotations
<br> Files: 
 - include/utils/S4_quaternion.H <br>
*/
template < class T = double >
class S4QuaternionTC : public R4VectorTC<T>
{
public:
    //**@name Access */
    //@{
    ///
    T Angle() const { return this->dw; }
    ///
    T X() const { return this->dx; }
    ///
    T Y() const { return this->dy; }
    ///
    T Z() const { return this->dz; }
    ///
    R3VectorTC<T> Axis() const { return R3VectorTC<T>( X(), Y(), Z() ); }
    //@}

    //**@name Converstion to/from rotation matrix */
    //@{
    ///
    R4MatrixTC<T> R4RotationMatrix( ) const;
    ///
    R3MatrixTC<T> R3RotationMatrix( ) const;
    //@}

    //**@name Making new rotations */
    //@{
    /// Inverse rotation (the Conjugate)
    S4QuaternionTC Conjugate() const ;
    /// Inverse q^{-1} = conj(a) / (a conj(a)) (i.e., Conjugate divided by length of vector)
    S4QuaternionTC Inverse() const ;
    /// returns a unit vector with same direction as this
    S4QuaternionTC Unit() const { return (*this) / ::Length(*this); }
    ///
    S4QuaternionTC Log( ) const ; 
    /// NOT communitive
    S4QuaternionTC & operator *= ( const S4QuaternionTC<T> & );
    ///
    S4QuaternionTC & operator /= ( const S4QuaternionTC<T> & );
    //@}

    /**@name Constructurs destructors */
    //@{
    ///
    S4QuaternionTC &operator=( const S4QuaternionTC<T> &in_quat ) { R4VectorTC<T>::operator=( in_quat ); return *this; }

    /// 
    S4QuaternionTC(const R4VectorTC<T> &vec): R4VectorTC<T>(vec) {};
    ///
    S4QuaternionTC(const S4QuaternionTC &qtn): R4VectorTC<T>(qtn) {};
    /// Axis of rotation, angle
    S4QuaternionTC(const R3VectorTC<T> &, const T & angle = 0);
    /// Rotation matrix
    S4QuaternionTC(const R3MatrixTC<T> &in_mat);
    /// x, y, z, angle
    S4QuaternionTC(const T & x=1, const T &y=0, const T &z=0, const T & angle = 0) : R4VectorTC<T>( x,y,z,angle) { }

    ///
    ~S4QuaternionTC() {};
    //@}
};

/**@name Quaternion operators */
//@{

/// NOT Communitive
template<class T>
inline S4QuaternionTC<T>
operator * ( const S4QuaternionTC<T> & q1, const S4QuaternionTC<T> & q2 )
{
    const T xx = q1.X() * q2.X();
    const T yy = q1.Y() * q2.Y();
    const T zz = q1.Z() * q2.Z();
    const T ww = q1.Angle() * q2.Angle();

    return S4QuaternionTC<T>(
        q1.Angle() * q2.X() + q1.X() * q2.Angle()   +   q1.Y() * q2.Z() - q1.Z() * q2.Y(),
        q1.Angle() * q2.Y() + q1.Y() * q2.Angle()   +   q1.Z() * q2.X() - q1.X() * q2.Z(),
        q1.Angle() * q2.Z() + q1.Z() * q2.Angle()   +   q1.X() * q2.Y() - q1.Y() * q2.X(),
        ww - xx - yy - zz );

}

/// NOT Communitive
template<class T>
inline R3VectorTC<T>
operator * ( const S4QuaternionTC<T> & q1, const R3VectorTC<T> & vec )
{
    const S4QuaternionTC<T> quat( vec[0], vec[1], vec[2], 0 );

    // Note: This is reversed in order to make RotationMatrix() * vec
    //       and this function return the same thing.
    // This also matches the sense of the R3Matrix::Rotation(axis, angle)
    // which is a rotation of the coordinate axis by the given angle (i.e,
    //  the point rotates the opposite direction)
    const S4QuaternionTC<T> quatRot = q1.Conjugate() * quat * q1;

    return quatRot.Axis();
}

/// NOT Communitive
template<class T>
inline S4QuaternionTC<T>
operator / ( const S4QuaternionTC<T> & q1, const S4QuaternionTC<T> & q2 )
{
    return q1 * q2.Inverse();
}

/// Equivalent of Lerp
template<class T>
inline S4QuaternionTC<T> SLerp( const S4QuaternionTC<T> &in_quat1, const S4QuaternionTC<T> &in_quat2, const T in_dT )
{
    const T dDot = Dot( in_quat1, in_quat2 ) / (Length( in_quat1 ) * Length( in_quat2 ));

    // sin will be zero
    //two quats are in the same direction
    if ( RNApproxEqual(  dDot , 1.0 , 1e-12 ) )
        return in_quat1 ;    
    
    //two quats differ by M_PI
    if ( RNApproxEqual(  dDot , -1.0, 1e-12 ) )
        return in_quat1 * (1.0 - in_dT) + in_quat2.Conjugate() * in_dT;

    const T dTheta = acos( dDot ) / 2.0;
    const T dW1 = sin( (1.0 - in_dT) * dTheta ) / sin( dTheta );
    const T dW2 = sin( in_dT * dTheta ) / sin( dTheta );

    if ( dDot > 0 ) {
        return in_quat1 * dW1 + in_quat2 * dW2;
    }
    return in_quat1 * dW1 + in_quat2.Conjugate() * dW2;
}

//@}

template<class T>
inline 
S4QuaternionTC<T> &S4QuaternionTC<T>::operator *= ( const S4QuaternionTC<T> & q2 )
{
    *this = (*this) * q2;
    return *this;
}

template<class T>
inline
S4QuaternionTC<T> &S4QuaternionTC<T>::operator /= ( const S4QuaternionTC<T> & q2 )
{
    *this = (*this) * q2.Inverse();
    return *this;
}


template < class T >
inline
S4QuaternionTC<T> S4QuaternionTC<T>::Conjugate() const 
{ 
    return S4QuaternionTC<T>( -this->dx, -this->dy, -this->dz, this->dw );
}

/// Inverse q^{-1} = conj(a) / (a conj(a))
template < class T >
inline
S4QuaternionTC<T> S4QuaternionTC<T>::Inverse() const 
{
    const S4QuaternionTC<T> quatConj = Conjugate();

  // const T dDot = R4VectorTC<T>::Length( *this );
  
    const T dDot = Length( *this ) ;

    return quatConj / dDot ;
//    return quatConj.R4VectorTC<T>::operator / ( dDot );
}

///logarithm of a quaternion : note this is only define for unit quaternions
template < class T >
inline
S4QuaternionTC<T> S4QuaternionTC<T>::Log() const 
{
   const S4QuaternionTC<T> log = Unit( ) ;
    
    const double dLen = Length( log.Axis () ) ;
    
    double dAng ;
    
    if( RNIsZero( log.Angle( ) ) )
    {
        dAng = M_PI * 0.5 ;
    }
    else
    {
        dAng = atan( dLen / log.Angle( ) ) ;
    }
    
    return S4QuaternionTC<T>( log.X() * dAng , log.Y() * dAng , log.Z() * dAng , 0.0 ) ;
    
}


template < class T >
inline S4QuaternionTC<T>::S4QuaternionTC(const R3VectorTC<T> & in_axis,
				                         const T & in_dAngle)
{
    const R3VectorTC<T> axis = ::Unit( in_axis );
    const T dSin = sin( in_dAngle / 2.0 );
    T dCos = cos( in_dAngle / 2.0 );

    if( RNIsZero( double (dCos) , 1e-5 ) )
        dCos = 0.0 ;

    this->dx = axis[0] * dSin;
    this->dy = axis[1] * dSin;
    this->dz = axis[2] * dSin;
    this->dw = dCos;
}

template < class T >
inline S4QuaternionTC<T>::S4QuaternionTC(const R3MatrixTC<T> & in_mat )
{
    //const T dTrace = in_mat.Trace( );
    const T dTrace = 1.0 + in_mat(0,0) + in_mat(1,1) + in_mat(2,2);
    if ( RNIsZero( dTrace ) ) {
        if ( in_mat(0,0) > in_mat(1,1) && in_mat(0,0) > in_mat(2,2) ) {
            const T dTraceSqrt = sqrt( 1.0 + in_mat(0,0) - in_mat(1,1) - in_mat(2,2) ) * 2.0;
            this->dx = 0.25 * dTraceSqrt;
            this->dy = ( in_mat(1,0) + in_mat(0,1) ) / dTraceSqrt;
            this->dz = ( in_mat(0,2) + in_mat(2,0) ) / dTraceSqrt;
            this->dw = ( in_mat(2,1) - in_mat(1,2) ) / dTraceSqrt;

        } else if ( in_mat(1,1) > in_mat(2,2) ) {
            const T dTraceSqrt = sqrt( 1.0 - in_mat(0,0) + in_mat(1,1) - in_mat(2,2) ) * 2.0;
            this->dx = ( in_mat(1,0) + in_mat(0,1) ) / dTraceSqrt;
            this->dy = 0.25 * dTraceSqrt;
            this->dz = ( in_mat(2,1) + in_mat(1,2) ) / dTraceSqrt;
            this->dw = ( in_mat(0,2) - in_mat(2,0) ) / dTraceSqrt;

        } else {
            const T dTraceSqrt = sqrt( 1.0 - in_mat(0,0) - in_mat(1,1) + in_mat(2,2) ) * 2.0;
            this->dx = ( in_mat(1,0) + in_mat(0,1) ) / dTraceSqrt;
            this->dy = ( in_mat(0,2) + in_mat(2,0) ) / dTraceSqrt;
            this->dz = 0.25 * dTraceSqrt;
            this->dw = ( in_mat(2,1) - in_mat(1,2) ) / dTraceSqrt;
        }
    } else {
        const T dTraceSqrt = sqrt( dTrace ) * 2.0;
        this->dx = -( in_mat(2,1) - in_mat(1,2) ) / dTraceSqrt;
        this->dy = -( in_mat(0,2) - in_mat(2,0) ) / dTraceSqrt;
        this->dz = -( in_mat(1,0) - in_mat(0,1) ) / dTraceSqrt;
        this->dw = 0.25 * dTraceSqrt;
    }
}

///
template < class T >
inline R4MatrixTC<T> S4QuaternionTC<T>::R4RotationMatrix( ) const
{
    const T qw = Angle();
    const T qx = X();
    const T qy = Y();
    const T qz = Z();

    const T qwqw = qw * qw;
    const T qxqx = qx * qx;
    const T qyqy = qy * qy;
    const T qzqz = qz * qz;

    const T qLenSq = qwqw + qxqx + qyqy + qzqz;

    R4MatrixTC<T> mat;

    mat(0,0) = 2.0 * (qLenSq / 2.0 - qyqy - qzqz) / qLenSq;
    mat(0,1) = (2.0 * (qx * qy + qw * qz)) / qLenSq;
    mat(0,2) = (2.0 * (qx * qz - qw * qy)) / qLenSq;
    mat(0,3) = 0.0;
    mat(1,0) = (2.0 * (qx * qy - qw * qz)) / qLenSq;
    mat(1,1) = (qLenSq - 2.0 * qxqx - 2.0 * qzqz) / qLenSq;
    mat(1,2) = (2.0 * (qw * qx + qy * qz)) / qLenSq;
    mat(1,3) = 0.0;
    mat(2,0) = (2.0 * (qw * qy + qx * qz)) / qLenSq;
    mat(2,1) = (2.0 * (qy * qz - qw * qx)) / qLenSq;
    mat(2,2) = (qLenSq - 2.0 * qxqx - 2.0 * qyqy) / qLenSq;
    mat(2,3) = 0.0;
    mat(3,0) = 0.0;
    mat(3,1) = 0.0;
    mat(3,2) = 0.0;
    mat(3,3) = 1.0;

    return mat;
}
///
template < class T >
inline R3MatrixTC<T> S4QuaternionTC<T>::R3RotationMatrix( ) const
{
    const T qw = Angle();
    const T qx = X();
    const T qy = Y();
    const T qz = Z();

    const T qwqw = qw * qw;
    const T qxqx = qx * qx;
    const T qyqy = qy * qy;
    const T qzqz = qz * qz;

    const T qLenSq = qwqw + qxqx + qyqy + qzqz;

    R3MatrixTC<T> mat;

    mat(0,0) = 2.0 * (qLenSq / 2.0 - qyqy - qzqz) / qLenSq;
    mat(0,1) = (2.0 * (qx * qy + qw * qz)) / qLenSq;
    mat(0,2) = (2.0 * (qx * qz - qw * qy)) / qLenSq;
    mat(1,0) = (2.0 * (qx * qy - qw * qz)) / qLenSq;
    mat(1,1) = (qLenSq - 2.0 * qxqx - 2.0 * qzqz) / qLenSq;
    mat(1,2) = (2.0 * (qw * qx + qy * qz)) / qLenSq;
    mat(2,0) = (2.0 * (qw * qy + qx * qz)) / qLenSq;
    mat(2,1) = (2.0 * (qy * qz - qw * qx)) / qLenSq;
    mat(2,2) = (qLenSq - 2.0 * qxqx - 2.0 * qyqy) / qLenSq;

    return mat;
}

typedef S4QuaternionTC<double> S4Quaternion;

#endif // IG_QUATERN_H

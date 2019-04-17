#ifndef __UTILS_COLOR_DEFS_H
#define __UTILS_COLOR_DEFS_H

#include <utils/Mesh_Array.H>
#include <utils/Rn_Defs.H>
#include <math.h>

/** \class UTILSColor Utils_Color.H utils/Utils_Color.H
 \ingroup UtilitiesColor
 \brief An RGB tripple.
 
    There are several default colors, or you can create your own with
    (0-1, 0-1, 0-1) RGB.

 <br>   Files: 
 - include/utils/Utils_Color.H
 - src/utils/utils/UTILSColor.cpp
*/

class UTILSColor {
public:
    /**@name Default colors */
    //@{
    typedef enum { 
        ///
        RED=0, 
        ///
        GREEN=1, 
        ///
        BLUE=2, 
        ///
        PURPLE=3, 
        ///
        YELLOW=4, 
        ///
        CYAN=5, 
        ///
        WHITE=6, 
        ///
        GREY=7,
        ///
        BLACK=8} Color;
    //@}

    /**@name Ramp types */
    //@{
    typedef enum { 
        ///
        INCREASING_RED = 0,
        ///
        DECREASING_RED, 
        ///
        INCREASING_GREEN, //=2
        ///
        DECREASING_GREEN, 
        ///
        INCREASING_BLUE, //=4
        ///
        DECREASING_BLUE,
        ///
        INCREASING_GREY, //=6
        ///
        DECREASING_GREY,
        ///
        INCREASING_RGB, //=8
        ///
        DECREASING_RGB,
        ///
        CONSTANT_LUMINANCE
    } RampType;
    //@}

protected:
    static const double sm_aadColors[ 9 ][ 3 ];
    double m_adC[ 3 ];

	static const R3Matrix rgbtolmsMatrix;
	static const R3Matrix lmstoatdMatrix;
	static const R3Matrix rgbtoatdMatrix;

	static const R3Matrix lmstorgbMatrix;
	static const R3Matrix atdtolmsMatrix;
	static const R3Matrix atdtorgbMatrix;

public:

    /**@name Accessors */
    //@{
    /// Pointer to rgb array
    const double *   Data( ) const { return m_adC; }
    ///
    double           operator[  ]( int in_iI ) const ;

    ///
    double &         operator[  ]( int in_iI ) ;

	///
	void             Print( ostream & os ) const;
    //@}

    /**@name Ramp colors */
    //@{
    ///
    static UTILSColor FromHSL( double in_dH, double in_dS, double in_dL );
    ///
    static UTILSColor Blend( const Array<UTILSColor> & in_aoColors, const Array<double> &in_adPerc );
    ///
    static UTILSColor Ramp( const double in_dValue, const double in_dMin, const double in_dMax, const RampType in_type );
    ///
    static UTILSColor RampMiddle( const double in_dValue, const double in_dMin, const double in_dMiddle, const double in_dMax, const RampType in_type );
    ///
    static UTILSColor LogRamp( const double in_dValue, const double in_dMin, const double in_dMax, const RampType in_type );
    /// Do a log ramp on left/right halves
    static UTILSColor LogRampMiddle( const double in_dValue, const double in_dMin, const double in_dMax, const RampType in_type );
    /// Go from number to dot color
    static UTILSColor Index( int );
    //@}

    ///
    WINbool         operator==( const UTILSColor & ) const;
    ///
    WINbool         operator==( Color ) const;

    /**@name Conversions */
    //@{
    ///
    void Clip( ) ;
    ///
	void ToCMY( double &, double &, double & ) const;
    ///
	void ToYIQ( double &, double &, double & ) const;
    ///
	void ToHSV( double &, double &, double & ) const;
    ///
	void ToHLS( double &, double &, double & ) const;
	///
	void ToXYZ( double &, double &, double & ) const;
	///
	void ToLUV( double &, double &, double & ) const;
	void ToLMS( double &, double &, double & ) const;
	///
	void ToATD( double &, double &, double & ) const;
	///
	void FromCMY( double , double , double  ) ;
    ///
	void FromYIQ( double , double , double  ) ;
    ///
	void FromHSV( double , double , double  ) ;
    ///
	void FromHLS( double , double , double  ) ;
	///
	void FromLUV( double , double , double  ) ;
	///
	void FromXYZ( double , double , double  ) ;
	///
	void FromLMS( double , double , double  ) ;
	///
	void FromATD( double , double , double  ) ;
	///
	void LUVtoXYZ( double , double , double, double &, double &, double &  ) const;
	///
	void XYZtoLUV( double , double , double, double &, double &, double &  ) const;
	///
	void LMStoATD( const double , const double , const double, double &, double &, double &  ) const;
	///
	void ATDtoLMS( const double , const double , const double, double &, double &, double &  ) const;
	//@}
    
    /**@name Constructors */
    //@{
    ///
    UTILSColor( );
    /// Construct with defined color
    inline UTILSColor( Color );
    /// RGB between 0 and 1
    UTILSColor( double, double, double );
    ///
    UTILSColor( const UTILSColor & );
    ///
    UTILSColor operator+( const UTILSColor & ) const;
    ///
    UTILSColor operator-( const UTILSColor & ) const;
    ///
    UTILSColor operator*( double ) const;
    ///
    UTILSColor & operator=( const UTILSColor & );
    ///
    inline UTILSColor & operator=( Color );
    ///
    ~UTILSColor( ) { }
    //@}

    /**@name Read/write */
    //@{
    ///
    void Write( ofstream & ) const;
    ///
    WINbool Read( ifstream & );
    ///
    void WriteBinary( ofstream & ) const;
    ///
    WINbool ReadBinary( ifstream & );
    ///
    static WINbool WriteBinary( const Array<UTILSColor> &, ofstream & );
    ///
    static WINbool ReadBinary( Array<UTILSColor> &, ifstream & );
    //@}

};

///
ostream & operator<<( ostream & os, const UTILSColor & c );

//------------------------------------------------------------------------------------------
inline double UTILSColor::operator[  ]( int in_iI ) const 
{ 
    ASSERT( in_iI >=0 && in_iI < 3 ); 
    return m_adC[ in_iI ]; 
}

//------------------------------------------------------------------------------------------
inline double & UTILSColor::operator[  ]( int in_iI ) 
{ 
    ASSERT( in_iI>=0 && in_iI < 3 ); 
    return m_adC[ in_iI ]; 
}

//------------------------------------------------------------------------------------------
inline
UTILSColor UTILSColor::operator+( const UTILSColor &in_col ) const
{
    return UTILSColor(in_col[0] + m_adC[0], in_col[1] + m_adC[1], in_col[2] + m_adC[2] );
}

//------------------------------------------------------------------------------------------
inline
UTILSColor UTILSColor::operator-( const UTILSColor &in_col ) const
{
    return UTILSColor(m_adC[0] - in_col[0], m_adC[1] - in_col[1], m_adC[2] - in_col[2] );
}

//------------------------------------------------------------------------------------------
inline 
UTILSColor UTILSColor::operator*( double in_dScl ) const
{
    return UTILSColor( m_adC[0] * in_dScl, m_adC[1] * in_dScl, m_adC[2] * in_dScl );
}

//------------------------------------------------------------------------------------------
inline
UTILSColor::UTILSColor( UTILSColor::Color _c ) 
{
    for ( int iI = 0; iI < 3; iI++ )
        m_adC[ iI ] = sm_aadColors[ _c ][ iI ];

}

//------------------------------------------------------------------------------------------
inline
UTILSColor & UTILSColor::operator=( UTILSColor::Color _c )
{ 
    for ( int iI = 0; iI < 3; iI++ )
        m_adC[ iI ] = sm_aadColors[ _c ][ iI ];

    return * this;

}//UTILSColor::operator=

#endif


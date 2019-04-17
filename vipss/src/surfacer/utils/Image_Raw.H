#ifndef __RAW_IMAGES_DEFS_H
#define __RAW_IMAGES_DEFS_H

#include <utils/Rn_Defs.H>
#include <utils/Utils_Color.H>
#include <utils/Rn_Polygon.H>


struct FIBITMAP;

/** Image manipulation.

    An image. See include/OpenGL/OGL_Images for OpenGL calls. Reads/writes bmps.

    <br>Files:
<dd><li> include/utils/Image_Raw.H
<dd><li> src/utils/image/IMAGERaw.cpp
*/
class IMAGERaw {
protected:
    typedef enum {
        RED = 0, 
        GREEN = 1,
        BLUE = 2,
        ALPHA = 3
    } MyRGBOffset;

    int           m_iXSize, m_iYSize, m_iStep;
    WINbool       m_bHasAlpha;
    Array<int>    m_aiXIndx;
    Array<float>  m_afImage;

public:
    /**@name Accessors */
    //@{
    ///
    WINbool HasAlpha() const { return m_bHasAlpha; }
    ///
    int     Width()  const { return m_iXSize; }
    ///
    int     Height() const { return m_iYSize; }
    ///
    int     NPixels() const { return m_iXSize * m_iYSize; }

    ///
    float     operator()( const int in_iX, const int in_iY, const int in_iRGB ) const;
    ///
    UTILSColor operator()( const R2Pt & ) const;
    ///
    UTILSColor operator()( const int in_iX, const int in_iY ) const;
    ///
    UTILSColor operator()( const int in_iPixel ) const;
    ///
    int        Index     ( const int in_iX, const int in_iY ) const;
    ///
    int        Index     ( const int in_iPixel ) const { return in_iPixel * m_iStep; }
    ///
    int        X         ( const int in_iIndex ) const;
    ///
    int        Y         ( const int in_iIndex ) const;
    ///
    float     Alpha( const int in_iX, const int in_iY ) const;
    ///
    float     Alpha( const int in_iPixel ) const;
    ///
    UTILSColor Interpolate( const R2Pt &in_pt ) const;
    ///
    Array<R2Pt> InterpolateIndex( const R2Pt &in_pt ) const;
    ///
    const float * Data() const { return &m_afImage[0]; }
    ///
    float * Data() { return &m_afImage[0]; }

    ///
    WINbool    IndexOk( const int in_iX, const int in_iY ) const;
    /// Does min/max
    UTILSColor Safe( const int in_iX, const int in_iY ) const;
    /// Does min/max
    float Safe( const int in_iX, const int in_iY, const int in_iRGB ) const;

    //@}

    /**@name Assignements */
    //@{
    ///
    void   SetSize( const int in_iX, const int in_iY, WINbool in_bHasAlpha = FALSE );

    ///
    float &operator()( const int in_iX, const int in_iY, const int in_iRGB ) ;
    ///
    void    Set       ( const int in_iPixel, const UTILSColor &in_pix );
    ///
    void    Set       ( const int in_iX, const int in_iY, const UTILSColor &in_pix );
    ///
    void    Set       ( const UTILSColor &in_pix );
    ///
    void    SetAlpha  ( const int in_iX, const int in_iY, const float in_dAlpha );
    ///
    void    SetAlpha  ( const int in_iPixel, const float in_dAlpha );
    ///
    void    SetAlpha  ( const float in_dAlpha );

    //@}

    /**@name Resampling */
    //@{
    ///
    IMAGERaw &operator=( const IMAGERaw & ) ;
    ///
    void Resample( IMAGERaw   & out_image, 
                   const R2Polygon & in_polyMine,
                   const R2Polygon & in_polyOut,
                   const float      in_dIgnoreAlpha = 1e30f ) const;
    ///
    void ResampleSimple( IMAGERaw   & out_image, 
                         const R2Polygon & in_polyMine,
                         const R2Polygon & in_polyOut,
                         Array<R2Pt_i>   & out_aipt) const;
    ///
    void ResampleEntire( IMAGERaw   & out_image, 
                         const R2Polygon & in_polyOut,
                         Array<R2Pt_i>   & out_aipt) const;
    ///
    void Resample( IMAGERaw &out_image, const R3Matrix &in_matMeToOut, const float in_dIgnoreAlpha = 1e30f ) const;
    //@}

    ///
    IMAGERaw();
    ///
    IMAGERaw( const int in_iX, const int in_iY, const WINbool in_bAlpha = FALSE );
    ///
    ~IMAGERaw();

    /**@name Read/write */
    //@{
	FIBITMAP *ConvertToFIB( ) const;
    /// Format determined by filename extension
    WINbool Read( const char *in_strFName );
    /// Format determined by filename extension
    WINbool Write( const char *in_strFName, const WINbool in_bUseAlpha = FALSE) const;
    //@}
};

inline int IMAGERaw::Index( const int in_iX, const int in_iY ) const
{
    ASSERT( in_iX >=0 && in_iX < Width() );
    return m_aiXIndx[ in_iY ] + m_iStep * in_iX;
}

inline float  IMAGERaw::operator()( const int in_iX, const int in_iY, const int in_iRGB ) const
{
    return m_afImage[ in_iRGB + Index(in_iX, in_iY) ];
}

inline float &IMAGERaw::operator()( const int in_iX, const int in_iY, const int in_iRGB ) 
{
    return m_afImage[ in_iRGB + Index(in_iX, in_iY) ];
}

inline UTILSColor IMAGERaw::operator()( const R2Pt &in_pt ) const
{
    const int iIndx =  Index( WINminmax( (int) (in_pt[0] * m_iXSize + 0.5), 0, m_iXSize - 1 ),
							  WINminmax( (int) (in_pt[1] * m_iYSize + 0.5), 0, m_iYSize - 1 ) );

    return UTILSColor( m_afImage[ RED + iIndx ],
                       m_afImage[ GREEN + iIndx ],
                       m_afImage[ BLUE + iIndx ] );
}

inline UTILSColor IMAGERaw::operator()( const int in_iX, const int in_iY ) const
{
    const int iIndx = Index(in_iX, in_iY);

    return UTILSColor( m_afImage[ RED + iIndx ],
                       m_afImage[ GREEN + iIndx ],
                       m_afImage[ BLUE + iIndx ] );
}

inline UTILSColor IMAGERaw::operator()( const int in_iIndx ) const
{
    const int iIndx = Index( in_iIndx );
    return UTILSColor( m_afImage[ RED + iIndx ],
                       m_afImage[ GREEN + iIndx ],
                       m_afImage[ BLUE + iIndx ] );
}

inline float IMAGERaw::Alpha( const int in_iX, const int in_iY ) const
{
    const int iIndx = Index(in_iX, in_iY);

    return m_afImage[ ALPHA + iIndx ];
}

inline float IMAGERaw::Alpha( const int in_iIndex ) const
{
    const int iIndx = Index( in_iIndex );
    return m_afImage[ ALPHA + iIndx ];
}

inline void IMAGERaw::Set( const int in_iX, const int in_iY, const UTILSColor &in_pix ) 
{
    const int iIndx = Index(in_iX, in_iY);

    m_afImage[ RED + iIndx ] = (float) in_pix[0];
    m_afImage[ GREEN + iIndx ] = (float) in_pix[1];
    m_afImage[ BLUE + iIndx ] = (float) in_pix[2];
}

inline void IMAGERaw::Set( const int in_iIndx, const UTILSColor &in_pix ) 
{
    const int iIndx = Index( in_iIndx );
    m_afImage[ RED + iIndx ] = (float) in_pix[0];
    m_afImage[ GREEN + iIndx ] = (float) in_pix[1];
    m_afImage[ BLUE + iIndx ] = (float) in_pix[2];
}

inline void IMAGERaw::SetAlpha( const int in_iX, const int in_iY, const float in_dAlpha ) 
{
    const int iIndx = Index(in_iX, in_iY);

    ASSERT( m_bHasAlpha == TRUE );
    m_afImage[ ALPHA + iIndx ] = in_dAlpha;
}

inline void IMAGERaw::SetAlpha( const int in_iIndx, const float in_dAlpha ) 
{
    ASSERT( m_bHasAlpha == TRUE );
    const int iIndx = Index( in_iIndx );
    m_afImage[ ALPHA + iIndx ] = in_dAlpha;
}

inline int IMAGERaw::X( const int in_iIndex ) const
{
    return (in_iIndex / m_iStep) % Width();
}

inline int IMAGERaw::Y( const int in_iIndex ) const
{
    return (in_iIndex / m_iStep) / Width();
}

inline WINbool IMAGERaw::IndexOk( const int in_iX, const int in_iY ) const
{
    if ( in_iX < 0 || in_iY < 0 )
        return FALSE;

    if ( in_iX >= Width() || in_iY >= Height() )
        return FALSE;

    return TRUE;
}

inline UTILSColor IMAGERaw::Safe( const int in_iX, const int in_iY ) const
{
    return (*this)( WINminmax( in_iX, 0, Width() - 1 ),
                    WINminmax( in_iY, 0, Height() - 1 ) );
}
/// Does min/max
inline float IMAGERaw::Safe( const int in_iX, const int in_iY, const int in_iRGB ) const
{
    return (*this)( WINminmax( in_iX, 0, Width() - 1 ),
                    WINminmax( in_iY, 0, Height() - 1 ), WINminmax(in_iRGB, 0, 2) );
}


#endif


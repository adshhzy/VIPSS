#ifndef __RAW_IMAGES_DEFS_1D_H
#define __RAW_IMAGES_DEFS_1D_H

#include <utils/Rn_Defs.H>
#include <utils/Utils_Color.H>
#include <utils/Rn_Polygon.H>

/** Image manipulation.

    A 1D image

    <br>Files:
<dd><li> include/utils/Image_Raw1D.H
<dd><li> src/utils/image/IMAGERaw1D.cpp
*/
class IMAGERaw1D {
protected:

    int           m_iXSize, m_iYSize;
    Array<int>    m_aiXIndx;
    Array<float>  m_afImage;

public:
    /**@name Accessors */
    //@{
    ///
    int     Width()  const { return m_iXSize; }
    ///
    int     Height() const { return m_iYSize; }
    ///
    int     NPixels() const { return m_iXSize * m_iYSize; }

    ///
    float     operator()( const int in_iX, const int in_iY ) const;
    ///
    float operator()( const R2Pt & ) const;
    ///
    float operator()( const int in_iPixel ) const { return m_afImage[in_iPixel]; }
    ///
    int        Index     ( const int in_iX, const int in_iY ) const;
    ///
    int        X         ( const int in_iIndex ) const;
    ///
    int        Y         ( const int in_iIndex ) const;
    ///
    float Interpolate( const R2Pt &in_pt ) const;
    ///
    Array<R2Pt> InterpolateIndex( const R2Pt &in_pt ) const;
    ///
    const float * Data() const { return &m_afImage[0]; }
    ///
    float * Data() { return &m_afImage[0]; }

    ///
    WINbool    IndexOk( const int in_iX, const int in_iY ) const;
    /// Does min/max
    float Safe( const int in_iX, const int in_iY ) const;

    //@}

    /**@name Assignements */
    //@{
    ///
    void   SetSize( const int in_iX, const int in_iY );

    ///
    float &operator()( const int in_iX, const int in_iY ) ;
    ///
    void    Set       ( const int in_iPixel, const float &in_pix );
    ///
    void    Set       ( const int in_iX, const int in_iY, const float &in_pix );
    ///
    void    Set       ( const float &in_pix );

    //@}

    /**@name Resampling */
    //@{
    ///
    IMAGERaw1D &operator=( const IMAGERaw1D & ) ;
    ///
    void Resample( IMAGERaw1D   & out_image, 
                   const R2Polygon & in_polyMine,
                   const R2Polygon & in_polyOut ) const;
    ///
    void ResampleSimple( IMAGERaw1D   & out_image, 
                         const R2Polygon & in_polyMine,
                         const R2Polygon & in_polyOut,
                         Array<R2Pt_i>   & out_aipt) const;
    ///
    void ResampleEntire( IMAGERaw1D   & out_image, 
                         const R2Polygon & in_polyOut,
                         Array<R2Pt_i>   & out_aipt) const;
    ///
    void Resample( IMAGERaw1D &out_image, const R3Matrix &in_matMeToOut ) const;
    //@}

    ///
    IMAGERaw1D();
    ///
    IMAGERaw1D( const int in_iX, const int in_iY);
    ///
    ~IMAGERaw1D();

    /**@name Read/write */
    //@{
    /// Format determined by filename extension
    WINbool Read( const char *in_strFName );
    /// Format determined by filename extension
    WINbool Write( const char *in_strFName ) const;
    //@}
};

inline int IMAGERaw1D::Index( const int in_iX, const int in_iY ) const
{
    ASSERT( in_iX >=0 && in_iX < Width() );
    return m_aiXIndx[ in_iY ] + in_iX;
}

inline float  IMAGERaw1D::operator()( const int in_iX, const int in_iY) const
{
    return m_afImage[ Index(in_iX, in_iY) ];
}

inline float &IMAGERaw1D::operator()( const int in_iX, const int in_iY) 
{
    return m_afImage[ Index(in_iX, in_iY) ];
}

inline float IMAGERaw1D::operator()( const R2Pt &in_pt ) const
{
    const int iIndx =  Index( WINminmax( (int) (in_pt[0] * m_iXSize + 0.5), 0, m_iXSize - 1 ),
							  WINminmax( (int) (in_pt[1] * m_iYSize + 0.5), 0, m_iYSize - 1 ) );

    return m_afImage[iIndx];
}


inline void IMAGERaw1D::Set( const int in_iX, const int in_iY, const float &in_pix ) 
{
    const int iIndx = Index(in_iX, in_iY);

    m_afImage[ iIndx ] = in_pix;
}

inline void IMAGERaw1D::Set( const int in_iIndx, const float &in_pix ) 
{
    m_afImage[ in_iIndx ] = in_pix;
}


inline int IMAGERaw1D::X( const int in_iIndex ) const
{
    return (in_iIndex) % Width();
}

inline int IMAGERaw1D::Y( const int in_iIndex ) const
{
    return (in_iIndex) / Width();
}

inline WINbool IMAGERaw1D::IndexOk( const int in_iX, const int in_iY ) const
{
    if ( in_iX < 0 || in_iY < 0 )
        return FALSE;

    if ( in_iX >= Width() || in_iY >= Height() )
        return FALSE;

    return TRUE;
}

inline float IMAGERaw1D::Safe( const int in_iX, const int in_iY ) const
{
    return (*this)( WINminmax( in_iX, 0, Width() - 1 ),
                    WINminmax( in_iY, 0, Height() - 1 ) );
}


#endif


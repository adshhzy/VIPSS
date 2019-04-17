#ifndef _UTILS_BITS_DEFS_H
#define _UTILS_BITS_DEFS_H

#include <WINSystemDefines.H>

#include <utils/Mesh_Array.H>

/* -----------------------  Constants  ------------------------------- */


/* -----------------------  Classes  ------------------------------- */

/** An array of bits.
 \ingroup UtilitiesFunc

  \brief This class has been replaced by the boost:dynamic_bitset<> class, located in externalLibraries/include/boost
 
   An array of bits which can be treated as an array of bools. Results in very compact storage of boolean arrays and very fast evaluation of bitwise operations, such as and-ing and or-ing of boolean arrays. <br>
 Files: <br>
  - include/utils/Utils_Bits.H <br>
  -  src/utils/utils/UTILBits.cpp <br>
  */
class UTILSBits  {
private:
    // 2^0, 2^1, ..., 2^32
    static unsigned long int   m_saiFlags[], m_siFlagsAll;
    int                        m_iNBits; // number of bits
    Array <unsigned long int>  m_aiData;
    
    int  NBitsPerInt( ) const { return sizeof( unsigned long int ) * 8; }
    
    inline WINbool CompareSizes( const UTILSBits & ) const;
    
public:
    /**@name Access to data */
    //@{
    /// Number of bits currently in use.
    inline int     NBits( ) const { return m_iNBits; }
    /// Number of bits set to TRUE.
    inline int     NBitsOn( ) const;
    ///
    inline int     AllOn( Array<int> &out_ai ) const;
    ///
    inline int     AllOn( Array<int> &out_ai, Array<int> &out_aiMap ) const;
    /// Will generate error if in_iIndx > m_iNBits
    inline void    SetBit( const int in_iIndx, const WINbool in_bValue );
    /// Expand the number of bits available if need be
    inline void    SetBitExpand( const int in_iIndx, const WINbool in_bValue );
    /// Set every bit to the input value
    inline UTILSBits &operator=( const WINbool );
    ///
    inline void       AddBit( const WINbool in_bValue );
    ///
    inline WINbool    operator[]( const int in_iIndx ) const;
    //@}

    /**@name Queries of which bits are on */
    //@{
    /** Which is the first TRUE bit?
        @return -1 if none set to TRUE */
    inline int FirstBitOn() const;
    /** Which is the first FALSE bit?
        @return -1 if none set to FALSE */
    inline int FirstBitOff() const;
    //@}

    /**@name Comparisions of bit arrays */
    //@{
    ///
    inline int        NOverlap( const UTILSBits & ) const;
    /// If a bit in b is on, then the corresponding bit in data is also on
    inline WINbool    IsSubset( const UTILSBits & ) const;
    /// Do the two arrays share bits ?
    inline WINbool    Overlap ( const UTILSBits & ) const;
    //@}

    /**@name Operations on bit arrays */
    //@{
    ///
    inline void    And     ( const UTILSBits &, int in_iStart = 0 ) ;
    ///
    inline void    Or      ( const UTILSBits &,  int in_iStart = 0) ;
    //@}


    /**@name Constructors and assignments */
    //@{
    ///
    inline UTILSBits &operator=( const UTILSBits &b );
    /// 
    inline UTILSBits( const int in_iNObjects, const WINbool in_bValue = FALSE );
    ///
    inline UTILSBits( const UTILSBits & );
    ///
    inline UTILSBits(  );
    ///
    inline UTILSBits( ifstream &in, const WINbool in_bBinary );
    //@}

    ///
    virtual ~UTILSBits() { }
    
    /**@name Read and Write */
    //@{
    ///
    void Write( ofstream & ) const;
    ///
    void Read ( ifstream & ) ;
    /// Safe for cross platform
    void WriteBinary( ofstream & ) const;
    /// Safe for cross platform
    void ReadBinary( ifstream & ) ;
    //@}

};

inline
WINbool UTILSBits::CompareSizes( const UTILSBits &b ) const
{ 
    return ( b.m_iNBits == m_iNBits ) ? TRUE : FALSE;
}

inline int UTILSBits::NBitsOn() const
{
    int n = 0;

    for (int i = 0; i < m_aiData.num( ); i++) {
        unsigned long int ui = m_aiData[i];
        if (( ui & m_siFlagsAll ) == m_siFlagsAll ) {
            n += NBitsPerInt( );
        } else {
            for (int j = 0; j < NBitsPerInt( ); j++) {
                n += ui % 2;
                ui = ui >> 1;
            }
        }
    }
    return n;
}

inline int UTILSBits::NOverlap( const UTILSBits &in_oBits ) const
{
    ASSERT( CompareSizes( in_oBits ) );
    int iNOverlap = 0;
    
    for (int i = 0; i < m_aiData.num( ); i++) {
        unsigned long int ui = m_aiData[i] & in_oBits.m_aiData[i];
        while ( ui ) {
            iNOverlap += ( ui % 2 );
            ui = ui >> 1;
        }
    }
    return iNOverlap;
}

inline int UTILSBits::FirstBitOn( ) const
{
    for (int i = 0; i < m_aiData.num(); i++) {
        if ( m_aiData[i] != 0 ) {
            for ( int j = 0; j < NBitsPerInt(); j++ ) {
                if ( m_saiFlags[j] & m_aiData[i] ) {
                    const int iN = i * NBitsPerInt() + j;
                    if ( iN >= m_iNBits ) return -1;
                    else                  return iN;
                }
            }
        }
    }
    return -1;
}


inline int UTILSBits::FirstBitOff( ) const
{
    if ( m_aiData.num() == 0 ) return -1;

    int i;
    for (i = 0; i < m_aiData.num(); i++) {
        if ( m_aiData[i] != m_siFlagsAll ) {
            for ( int j = 0; j < NBitsPerInt(); j++ ) {
                if ( ! ( m_saiFlags[j] & m_aiData[i] ) ) {
                    const int iN = i * NBitsPerInt() + j;
                    if ( iN >= m_iNBits ) return -1;
                    else                  return iN;
                }
            }
        }
    }
    const int iN = ( m_aiData.num() - 1 ) * NBitsPerInt();
    for ( i = 0; i < m_iNBits - iN; i++ ) {
        if ( ! ( m_saiFlags[i] & m_aiData.last() ) )
          return iN + i;
    }
    return -1;
}
inline void UTILSBits::And( const UTILSBits &in_oBits, const int ciStart )
{
    ASSERT( CompareSizes( in_oBits ) );

    for (int i = ciStart; i < m_aiData.num( ); i++)
      m_aiData[i] = m_aiData[i] & in_oBits.m_aiData[i];
}

inline void UTILSBits::Or( const UTILSBits &in_oBits, const int ciStart )
{
    ASSERT( CompareSizes( in_oBits ) );

   for (int i = ciStart; i < m_aiData.num( ); i++)
     m_aiData[i] = m_aiData[i] | in_oBits.m_aiData[i];
}

inline WINbool UTILSBits::IsSubset( const UTILSBits &in_oBits ) const
{
    if ( CompareSizes( in_oBits ) == FALSE ) return FALSE;

    for (int i = 0; i < m_aiData.num( ); i++)
      if (( m_aiData[i] & in_oBits.m_aiData[i] ) != in_oBits.m_aiData[i]) 
        return FALSE;
    
    return TRUE;
}

inline WINbool UTILSBits::Overlap( const UTILSBits &in_oBits ) const 
{
    if ( CompareSizes( in_oBits ) == FALSE ) return FALSE;
    
    for (int i = 0; i < m_aiData.num( ); i++)
      if ( ( m_aiData[i] & in_oBits.m_aiData[i] ) != 0 ) 
        return TRUE;

    return FALSE;
}

inline void UTILSBits::SetBit( const int in_iIndx, const WINbool m_bValue ) 
{
    ASSERT( in_iIndx < m_iNBits );

    const int iIndx = in_iIndx / NBitsPerInt( );
    const int iBitIndx = in_iIndx - ( iIndx * NBitsPerInt() );

    unsigned long int &ui = m_aiData[ iIndx ];
    
    if ( m_bValue == TRUE) ui = ui | m_saiFlags[ iBitIndx ];
    else                   ui = ui & ( m_siFlagsAll ^ m_saiFlags[ iBitIndx ] );
}

inline void UTILSBits::SetBitExpand( const int in_iIndx, const WINbool in_bValue ) 
{
    const int iIndx = in_iIndx / NBitsPerInt( );
    const int iBitIndx = in_iIndx - ( iIndx * NBitsPerInt() );

    m_iNBits = WINmax( in_iIndx+1, m_iNBits );
    m_aiData.needSet( iIndx+1, 0 );
    unsigned long int &ui = m_aiData[ iIndx ];
    
    if ( in_bValue == TRUE) ui = ui | m_saiFlags[ iBitIndx ];
    else                    ui = ui & ( m_siFlagsAll ^ m_saiFlags[ iBitIndx ] );
}

inline void UTILSBits::AddBit( const WINbool in_bValue ) 
{
    if ( m_iNBits % NBitsPerInt( ) == 0 )
      m_aiData += 0;

    m_iNBits++;
    SetBit( m_iNBits-1, in_bValue );
}

inline UTILSBits &UTILSBits::operator=( const WINbool in_bValue ) 
{
    if ( m_aiData.num() == 0 ) return *this;

    const int iN = sizeof( unsigned long int ) * m_aiData.num( );
    if ( in_bValue == FALSE) 
      memset( &m_aiData[0], 0, iN );
    else                    
      memset( &m_aiData[0], 255, iN );

    return *this;
}

inline UTILSBits &UTILSBits::operator=( const UTILSBits &in_oBits ) 
{
    m_aiData = in_oBits.m_aiData;
    m_iNBits = in_oBits.m_iNBits;

    return *this;
}

inline WINbool UTILSBits::operator[]( const int in_iIndx ) const 
{
    const int iIndx = in_iIndx / NBitsPerInt( );
    return 
      (( m_aiData[ iIndx ] >> ( in_iIndx % NBitsPerInt( ) ) ) & 1) == 1 ? TRUE : FALSE;
}

inline UTILSBits::UTILSBits( const int in_iN, const WINbool in_bValue )
: m_iNBits( in_iN ), m_aiData( (in_iN / NBitsPerInt( )) + 1 )
{
    (*this) = in_bValue;
}

inline UTILSBits::UTILSBits(  )
: m_iNBits( 0 ), m_aiData( 1 )
{
    (*this) = FALSE;
}

inline UTILSBits::UTILSBits( const UTILSBits &in_oBits )
{
    (*this) = in_oBits;
}

inline UTILSBits::UTILSBits( ifstream &in, const WINbool in_bBinary )  
{
    if ( in_bBinary ) ReadBinary( in );
    else              Read( in );
}

inline int UTILSBits::AllOn( Array<int> &out_ai ) const
{
    out_ai.need(0);

    const int iNBytes = m_iNBits / NBitsPerInt();
    for ( int i = 0; i < iNBytes; i++ ) {
        if ( m_aiData[i] == 0 ) continue;

        const int iByte = m_aiData[i];
        const int ciStart = i * NBitsPerInt();
        for (int j = 0; j < NBitsPerInt( ); j++) {
            const int iEval = ( iByte >> ( j ) ) & 1;

            if ( iEval ) out_ai += ciStart + j;
        }
    }

    const int iLeft = iNBytes * NBitsPerInt();
    for (int j = iLeft; j < m_iNBits; j++) {
        if ( (*this)[j] == TRUE ) out_ai += j;
    }

    return out_ai.num();
}

inline int UTILSBits::AllOn( Array<int> &out_ai, Array<int> &out_aiMap ) const
{
    out_ai.need(0);
    out_aiMap.need( m_iNBits );

    const int iNBytes = m_iNBits / NBitsPerInt();
    for ( int i = 0; i < iNBytes; i++ ) {
        if ( m_aiData[i] == 0 ) continue;

        const int iByte = m_aiData[i];
        const int ciStart = i * NBitsPerInt();
        for (int j = 0; j < NBitsPerInt( ); j++) {
            const int iEval = ( iByte >> ( j ) ) & 1;
            out_aiMap[ciStart + j] = out_ai.num();

            if ( iEval ) 
                out_ai += ciStart + j;
        }
    }

    const int iLeft = iNBytes * NBitsPerInt();
    for (int j = iLeft; j < m_iNBits; j++) {
        out_aiMap[j] = out_ai.num();
        if ( (*this)[j] == TRUE ) out_ai += j;
    }

    return out_ai.num();
}

#endif


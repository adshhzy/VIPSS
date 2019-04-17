/*************
 * UTILSElement routines
 *************/

/****************************************************/
inline
void UTILSSetUtilsElement( UTILSElement    &      in_oElem, 
                           int                    in_iId, 
                           UTILSElement::TypeEnum in_oType)
{
    in_oElem.m_iId = in_iId;
    in_oElem.m_oType = in_oType;
}

/****************************************************/
inline
void UTILSValidate( UTILSElement &in_oElem )
{
    in_oElem.m_bValid = TRUE;
}

/****************************************************/
inline
void UTILSInvalidate( UTILSElement &in_oElem )
{
    in_oElem.m_bValid = FALSE;
}

/****************************************************/
template<class Data>
inline
void  UTILSElementArray<Data>::UnIdElement( Data &in_oData ) 
{ 
    UTILSSetUtilsElement ( in_oData, -1, UTILSElement::NONE ); 
    UTILSInvalidate( in_oData );
}

/****************************************************/
inline
UTILSElement::UTILSElement( )  
{ 
    m_iId    = -1; 
    m_oType  = UTILSElement::NONE; 
    m_bValid = FALSE;
}

/****************************************************/
inline
UTILSElement::~UTILSElement( ) 
{ 
    m_iId    = -1; 
    m_oType  = UTILSElement::NONE; 
    m_bValid = FALSE;
}

/*************
 * UTILSElementArray access routines
 *************/

template< class Data >
inline
int UTILSElementArray<Data>::CGetAll( Array< const Data * > &out_aopData ) const
{
    out_aopData.need( N() );

    int iCount = 0;
    for (int i = 0; i < NMax(); i++)
        if ( m_aopData[i] && m_aopData[i]->UTILSValid() )
            out_aopData[iCount++] = m_aopData[i];

    ASSERT( iCount == out_aopData.num() );
    return iCount;

}

/****************************************************/
template< class Data >
inline
int UTILSElementArray<Data>::GetAll( Array< Data * > &out_aopData ) 
{
    out_aopData.need( N() );

    int iCount = 0;
    for (int i = 0; i < NMax(); i++)
        if ( m_aopData[i] && m_aopData[i]->UTILSValid() )
            out_aopData[iCount++] = m_aopData[i];

    ASSERT( iCount == out_aopData.num() );
    return iCount;

}

/****************************************************/
template<class Data>
inline
const Data &
UTILSElementArray<Data>::operator[ ]( int in_iID ) const
{ 
    return *m_aopData[ FindData( in_iID ) ] ; 
}

/****************************************************/
template<class Data>
inline
Data &
UTILSElementArray<Data>::operator[ ]( int in_iID ) 
{ 
    return *m_aopData[ FindData( in_iID ) ] ; 
}

/****************************************************/
template<class Data>
inline
int UTILSElementArray<Data>::Id( const Data &in_oData ) const 
{ 
    return in_oData.UTILSId( ); 
}

/****************************************************/
template<class Data>
inline
WINbool UTILSElementArray<Data>::Exists( int in_iID ) const
{
    if ( in_iID >= m_aopData.num() || in_iID == -1 ) return FALSE;

    if ( m_aopData[ in_iID ] == NULL ) return FALSE;

    ASSERT ( m_aopData[ in_iID ]->UTILSId( ) == in_iID );

    if ( m_aopData[ in_iID ]->UTILSValid() ) return TRUE;

    return FALSE;
}

/****************************************************/
template<class Data>
inline
int UTILSElementArray<Data>::FindData( int in_iID ) const
{
    ASSERT( m_aopData[ in_iID ]->UTILSId( ) == in_iID );

    return in_iID;
}



/*************
 * UTILSElementArray add/remove routines
 *************/

template<class Data>
inline
Data &UTILSElementArray<Data>::NewElement( ) 
{
    int   iN    = m_oBits.FirstBitOff();
    Data *oData = NULL;

    ASSERT( m_bAllocateData == TRUE );

    if ( iN == -1 ) {
        iN = m_oBits.NBits();
        m_oBits.AddBit( TRUE );
        oData = new Data;
        m_aopData += oData;
    } else {
        oData = m_aopData[iN];
    }

    m_oBits.SetBit(iN, TRUE);

    UTILSSetUtilsElement ( *oData, iN, m_oType );
    UTILSValidate        ( *oData );

    m_iCount++;

    return *oData;
}

/****************************************************/

template<class Data>
inline
WINbool UTILSElementArray<Data>::DeleteElement( Data &in_oData ) 
{
    int iN = in_oData.UTILSId();
    WINbool out_bRet = (in_oData == *m_aopData[iN]) ? TRUE : FALSE;

    ASSERT( m_bAllocateData == TRUE );

    if ( in_oData.UTILSValid() == FALSE )
      out_bRet = FALSE;

    UTILSInvalidate( in_oData );
    m_oBits.SetBit( iN, FALSE );

    m_iCount--;

    ASSERT( out_bRet );
    return out_bRet;
}

/****************************************************/

template<class Data>
inline
Data &UTILSElementArray<Data>::AddElement( Data *in_oElem ) 
{
    int   iN    = m_oBits.FirstBitOff();
    Data *oData = in_oElem;

    ASSERT( m_bAllocateData == FALSE );

    if ( iN == -1 ) {
        iN = m_oBits.NBits();
        m_oBits.AddBit( TRUE );
        m_aopData += oData;
    } else {
        m_aopData[iN] = oData;
    }

    m_oBits.SetBit(iN, TRUE);

    UTILSSetUtilsElement ( *oData, iN, m_oType );
    UTILSValidate        ( *oData );

    m_iCount++;

    return *oData;
}

/****************************************************/
template<class Data>
inline
Data &UTILSElementArray<Data>::AddElementWithDataId( Data *in_oElem ) 
{
    int   iId   = in_oElem->UTILSId();

    ASSERT( m_bAllocateData == FALSE );

    ASSERT( iId >= m_oBits.NBits() || m_oBits[iId] == FALSE );
    ASSERT( in_oElem->UTILSType() == m_oType );

    m_oBits.SetBitExpand( iId, TRUE );
    m_aopData.needSet( WINmax( iId+1, m_aopData.num() ), NULL );
    m_aopData[iId] = in_oElem;

    m_iCount++;

    return *in_oElem;
}

/****************************************************/
template<class Data>
inline
Data &UTILSElementArray<Data>::AddElementWithDataId( Data *in_oElem, int in_iId ) 
{
    ASSERT( m_bAllocateData == FALSE );

    ASSERT( in_iId >= m_oBits.NBits() || m_oBits[in_iId] == FALSE );

    UTILSSetUtilsElement ( *in_oElem, in_iId, m_oType );
    UTILSValidate        ( *in_oElem );

    m_oBits.SetBitExpand( in_iId, TRUE );
    m_aopData.needSet( WINmax( in_iId+1, m_aopData.num() ), NULL );
    m_aopData[in_iId] = in_oElem;

    m_iCount++;

    return *in_oElem;
}


/****************************************************/
/* Don't unset the in_oData type/valid/id since the calling routine may still
 * want it - use UnIdElement explicitly if unsetting is desired
 */
template<class Data>
inline
WINbool UTILSElementArray<Data>::RemoveElement( Data &in_oData ) 
{
    int iN = in_oData.UTILSId();
    WINbool out_bRet = (in_oData == *m_aopData[iN]) ? TRUE : FALSE;

    ASSERT( m_bAllocateData == FALSE );

    if ( in_oData.UTILSValid() == FALSE )
      out_bRet = FALSE;

    m_aopData[ iN ] = NULL;
    m_oBits.SetBit( iN, FALSE );

    m_iCount--;

    ASSERT( out_bRet );
    return out_bRet;
}


/****************************************************/
template<class Data>
inline
void UTILSElementArray<Data>::RemoveAll( ) 
{
    if ( m_bAllocateData == FALSE ) {
        for (int i = 0; i < m_aopData.num(); i++ ) {
            ASSERT( m_aopData[i] == NULL );
        }
    } else {
        for (int i = 0; i < m_aopData.num(); i++ ) {
            UTILSInvalidate( *m_aopData[i] );
        }
    }
    m_iCount = 0;
    m_oBits = FALSE;
}

/*************
 * Renumber
 *************/
/* Checks to see that the data is ok while it's renumbering.
   Returns the new id number (out_aiMap.num() == old size,
   out_aiMap[i] == -1 if element i was not valid).
   */
template<class Data>
inline
void UTILSElementArray<Data>::Renumber( Array<int> &out_aiMap ) 
{
    ASSERT( m_aopData.num() == m_oBits.NBits() );
    out_aiMap.need( m_aopData.num() );

    Array< Data * > aopData;
    Array< Data * > aopDead;
    for ( int i = 0; i < m_aopData.num(); i++ ) {
        if ( m_aopData[i] == NULL ) {
            ASSERT( m_bAllocateData == FALSE );
            ASSERT( m_oBits[i]      == FALSE );
            out_aiMap[i] = -1;
            continue;
        }
        ASSERT( m_aopData[i]->UTILSType() == m_oType );
        ASSERT( m_aopData[i]->UTILSId()   == i );

        if ( m_oBits[i] == FALSE ) {
            ASSERT( m_aopData[i]->UTILSValid() == FALSE );
            aopDead += m_aopData[i];
            out_aiMap[i] = -1;
        } else {
            ASSERT( m_aopData[i]->UTILSValid() == TRUE );
            UTILSSetUtilsElement( *m_aopData[i], aopData.num(), m_oType );
            out_aiMap[i] = aopData.num();

            aopData += m_aopData[i];
        }
    }

    if ( m_bAllocateData == TRUE ) {
        for ( int i = 0; i < aopDead.num(); i++ )
            delete aopDead[i];
    } else {
        for ( int i = 0; i < aopDead.num(); i++ ) 
            ASSERT( aopDead[i] == NULL );
    }
    UTILSBits oBits( aopData.num(), TRUE );
    
    m_oBits = oBits;
    m_aopData = aopData;

    m_aopData.squeeze();
}

/*************
 * UTILSElementArray constructor/destructor routines
 *************/
template<class Data>
inline
UTILSElementArray<Data>::~UTILSElementArray( ) 
{
    RemoveAll();
    m_oType = UTILSElement::NONE; 

    if ( m_bAllocateData == TRUE ) {
        for (int i = 0; i < m_aopData.num(); i++)
            delete m_aopData[i];
    }
}

/****************************************************/
template<class Data>
inline
UTILSElementArray<Data>::UTILSElementArray( UTILSElement::TypeEnum in_oType,
                                            WINbool in_bDoAlloc )
: m_oType( in_oType ), 
  m_oBits(0), 
  m_iCount(0), 
  m_bAllocateData( in_bDoAlloc )
{ 
}

/*************
 * UTILSElementArray read/write/compare routines
 *************/
template<class Data>
inline
WINbool
UTILSElementArray<Data>::operator==( const UTILSElementArray &in_oE ) const
{ 
    if ( m_oType  != in_oE.m_oType  ) return FALSE;
    if ( m_iCount != in_oE.m_iCount ) return FALSE;

    for ( int i = 0; i < m_aopData.num(); i++ ) {
        if ( m_aopData[i] && UTILSValid( m_aopData[i] ) ) {
            if ( (*in_oE[ m_aopData[i]->Id() ]) != (*m_aopData[i]) )
                return FALSE;
        }
    }
    return TRUE;
}

/****************************************************/
template<class Data>
inline
void
UTILSElementArray<Data>::Read( ifstream &in,
                               WINbool   in_bForceAlloc )
{
    char strTemp[255];
    in >> strTemp;
    ASSERT( !strcmp(strTemp, "UTILSElementArray") );
    ASSERT( in.good() );
    
    int iType = -1;
    in >> m_iCount >> iType;
    m_oType = ( UTILSElement::TypeEnum) iType;
    m_oBits.Read( in );

    if ( m_aopData.num() < m_oBits.NBits() ) {
        m_aopData.need( m_oBits.NBits() );
    }

    for ( int i = 0; i < m_oBits.NBits(); i++ ) {
        if ( m_bAllocateData == TRUE ) {
            m_aopData[i] = new Data;
            UTILSSetUtilsElement( *m_aopData[i], i, m_oType );
        }

        if ( m_oBits[i] == TRUE ) {
            if (  in_bForceAlloc == TRUE ) {
                m_aopData[i] = new Data;
                UTILSSetUtilsElement( *m_aopData[i], i, m_oType );
            }

            UTILSValidate( *m_aopData[i] );
            m_aopData[i]->Read(in);
            ASSERT( in.good() );
        } else if ( m_bAllocateData == FALSE ) {
            m_aopData[i] = NULL;
        }
    }
    ASSERT( in.good() );
    in.ignore(2, '\n');
    in >> strTemp;
    ASSERT( !strcmp(strTemp, "UTILSElementArrayEnd") );
    ASSERT( in.good() );
}

/****************************************************/
template<class Data>
inline
void
UTILSElementArray<Data>::ReadBinary(ifstream &in,
                                    WINbool   in_bForceAlloc)
{
    char strTemp[255];

    in >> strTemp;
    in.ignore( 1, '\n' );
    ASSERT( !strcmp(strTemp, "UTILSElementArrayBinary") );
    
    int iType = -1;
    UTILSIo::ReadBinary( m_iCount, in );
    UTILSIo::ReadBinary( iType, in );

    m_oType = ( UTILSElement::TypeEnum) iType;

    m_oBits.ReadBinary( in );

    if ( m_aopData.num() < m_oBits.NBits() ) {
        m_aopData.need( m_oBits.NBits() );
    }

    for ( int i = 0; i < m_oBits.NBits(); i++ ) {
        if ( m_bAllocateData == TRUE ) {
            m_aopData[i] = new Data;
            UTILSSetUtilsElement( *m_aopData[i], i, m_oType );
        }

        if ( m_oBits[i] == TRUE ) {
            if (  in_bForceAlloc == TRUE ) {
                m_aopData[i] = new Data;
                UTILSSetUtilsElement( *m_aopData[i], i, m_oType );
            }

            UTILSValidate( *m_aopData[i] );
            m_aopData[i]->ReadBinary(in);
        } else if ( m_bAllocateData == FALSE ) {
            m_aopData[i] = NULL;
        }
    }
    in >> strTemp;
    in.ignore( 2, '\n' );
    ASSERT( !strcmp(strTemp, "UTILSElementArrayBinaryEnd") );
}

/****************************************************/
/* Write out all valid elements */
template<class Data>
inline
void
UTILSElementArray<Data>::Write(ofstream &out) const
{
    out << " UTILSElementArray " << m_iCount << " " << m_oType << "\n";
    m_oBits.Write(out);

    for ( int i = 0; i < m_oBits.NBits(); i++ ) {
        if ( m_oBits[i] == TRUE ) {
            m_aopData[i]->Write(out);
        }
    }
    out << "\nUTILSElementArrayEnd\n";
}


/****************************************************/
/* Write out all valid elements */
template<class Data>
inline
void
UTILSElementArray<Data>::WriteBinary(ofstream &out) const
{
    out << "\nUTILSElementArrayBinary\n";
    UTILSIo::WriteBinary( m_iCount, out );
    UTILSIo::WriteBinary( (int) m_oType, out );

    m_oBits.WriteBinary(out);

    for ( int i = 0; i < m_oBits.NBits(); i++ ) {
        if ( m_oBits[i] == TRUE ) {
            m_aopData[i]->WriteBinary(out);
        }
    }
    out << "\nUTILSElementArrayBinaryEnd\n";
}

/****************************************************/
/* Check that everything's as it should be */
template<class Data>
inline
WINbool
UTILSElementArray<Data>::Check() const
{
    ASSERT ( m_oBits.NBits() == m_aopData.num() );

    for (int i = 0; i < m_aopData.num(); i++) {
        if ( m_oBits[i] == TRUE ) {
            ASSERT( m_aopData[i] );
            ASSERT( m_aopData[i]->UTILSId() == i );
        } else {
            if ( m_bAllocateData == FALSE ) {
                ASSERT( m_aopData[i] == NULL );
            } else {
                ASSERT( m_aopData[i]->UTILSId() == -1 );
            }
        }
    }
    return TRUE;
}

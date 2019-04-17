#ifndef _UTILS_ELEMENTS_DEFS_H
#define _UTILS_ELEMENTS_DEFS_H

#include <utils/Utils_Bits.H>
#include <utils/Utils_Io.H>

/** \defgroup UtilitiesFuncElements Array of elements with unique ids
 \ingroup UtilitiesFunc

  Not really used anymore.

 The purpose of these classes is to provide an array of elements
  where each element has a unique id. The elements are (in general)
  invalidated instead of deleted, with new elements coming from the
  invalidated ones. <br>

  Each element in the array is subclassed off of a UTILSElement. This
  gives each element an id, a type, and a valid/invalid flag. <br>

  The array can be formed in one of two ways; either the array
  allocates the data or the caller creates the data and the array just
  manages it. The latter allows for two arrays each pointing to the
  same data. This is useful if one array contains data elements whose
  type is inherited from the type in the other array. (See
  manifolds/MFMeshGenerator for an example of this.) <br>
<br>   Files: 
-     include/utils/Utils_Element.H
-     include/utils/Utils_Element_i.H
-     utils/utils/UTILSElement.cpp
  */
//@{

/** Base class for elements in the UTILSElement_array class.

<br>   Data: 
-     A unique id, assigned by the holding UTILSElement_array
-     A type, of the same type as the holding UTILSElement_array
-     A valid flag, indicates if this element is currently in use
                   or if it's been "killed"
 */
class UTILSElement {
public:
    /**@name Element type
       If the program is using multiple arrays, each array can have a
       different type. This prevents elements of one array from being
       plunked in another. A sort of meta-level type checking. <br>
       Add in types as needed. <br>
       When the holding array destroys the element (as opposed to just
       marking it as invalid) the type will be set to NONE. <br> */
    typedef enum {
        ///
        NONE = -1,
        /// PickMesh
        VERTEX,
        /// PickMesh
        EDGE,
        /// PickMesh
        FACE,
        /// MFManifold
        VERTEX_CHART,
        /// MFManifold
        EDGE_CHART,
        /// MFManifold
        FACE_CHART,
        /// MFEmbedding
        VERTEX_EMBEDDING,
        /// MFEmbedding
        EDGE_EMBEDDING,
        /// MFEmbedding
        FACE_EMBEDDING,
        /// MFEmbedViewDep
        COLOR,
        /// MFEmbedTexture
        TEXTURE,
        ///
        OTHER
    } TypeEnum;  
    
protected:
    WINbool          m_bValid;
    int              m_iId;
    TypeEnum         m_oType;
    
    // Lets UTILSElement_array manipulate the id, type, and valid flags
    friend void UTILSSetUtilsElement ( UTILSElement &, int, TypeEnum );
    friend void UTILSValidate        ( UTILSElement & );
    friend void UTILSInvalidate      ( UTILSElement & );

public:
    /**@name Accessor macros 
       UTILS is put in front to allow other classes to make their own
       Id() method. */
    //@{
    ///
    WINbool        UTILSValid( )  const { return m_bValid; }
    ///
    TypeEnum       UTILSType( )   const { return m_oType; }
    ///
    int            UTILSId( )     const { return m_iId; }
    //@}

    ///
    UTILSElement( );
    ///
    virtual ~UTILSElement( ) ;

    /**@name Default read, write 
       Id, valid and type are written by UTILSElementArray
       Override these to let the inheriting class write data */
    //@{
    ///
    virtual void Read(ifstream &in) { } 
    ///
    virtual void Write(ofstream &out) const { }
    ///
    virtual void ReadBinary(ifstream &in) { } 
    ///
    virtual void WriteBinary(ofstream &out) const { }
    //@}
};


/** Array of elements with unique ids.
  class Data must inherit from UTILSElement, i.e.
       class Data : public UTILSElement { }; 
<br>
 * The array can be constructed in one of two ways; either the data is
 * allocated by the array ( m_bAllocateData == TRUE ) in which case
 * the routines NewElement and DeleteElement are used, or the data is
 * allocated by the caller, in which case AddElement and RemoveElement
 * are used.
<br>

 The array has a type which is checked agains the elements which are
 passed in and returned; if all the UTILSElementArrays are created
 with different types this will prevent one element as masquerading as
 an element of another array.
<br>

 For all elements, m_aopData[ iId ]->UTILSId() == iId. This is
 checked.
<br>

 m_oBits is a bit array which keeps track of which elements are
 "dead", and hence can be re-used 
<br>
*/
template< class Data >
class UTILSElementArray {
protected:
    const WINbool                 m_bAllocateData;
    UTILSElement::TypeEnum        m_oType;    // All elements have this type
    UTILSBits                     m_oBits;    // Which items in m_aopData are in use
    Array< Data * >               m_aopData;  // The actual data
    int                           m_iCount;   // How many valid elements there are

    // To set the data in the individual elements
    friend void UTILSSetUtilsElement ( UTILSElement &, int in_iId, UTILSElement::TypeEnum );
    friend void UTILSValidate        ( UTILSElement & );
    friend void UTILSInvalidate      ( UTILSElement & );

    int FindData( int in_iID ) const;

    // Dummy routines
    UTILSElementArray ( const UTILSElementArray & );
    UTILSElementArray & operator=( const UTILSElementArray & )const ;

public:
    /**@name Accessor methods */
    //@{
    /// An actual count
    int                 N() const { return m_iCount; }
    /// An upper bound
    int                 NMax() const { return m_aopData.num(); }
    /// Does it exist and is it valid?
    WINbool             Exists      ( int in_iId ) const;
    /// Returns data with id in_iId
    const Data         &operator[ ] ( int in_iID ) const;
    /// Returns data with id in_iId
          Data         &operator[ ] ( int in_iID ) ;

    /// Get all valid/existing elements
    int                 CGetAll( Array< const Data * > & ) const;
    /// Get all valid/existing elements
    int                 GetAll ( Array< Data * > & ) ;

    /// Get the id of the data
    int                 Id          ( const Data &in_oData ) const;
    //@}

    /**@name Adding and deleting data with allocation 
       m_bAllocateData == TRUE <br>
       UTILSElement creates and deletes data */
    //@{
    ///
    Data               &NewElement( );
    ///
    WINbool             DeleteElement( Data & );
    //@}

    /**@name Adding and deleting without allocation
       Element created and deleted by calling application. <br>
       m_bAllocateData == FALSE <br>
       Add in at the id location given by the data element or by the id */
    //@{
    /// Use the id in the data element
    Data               &AddElementWithDataId( Data *  );
    /// Use the given id (sets id)
    Data               &AddElementWithDataId( Data *, int in_iId  );

    /// Add in at the first available location (sets id)
    Data               &AddElement( Data * );
    ///
    WINbool             RemoveElement( Data & );
    // Set id to -1, valid to false, and type to NONE
    void                UnIdElement( Data &in_oData ) ;
    //@}

    /// Will "delete" if m_bAllocateData is TRUE
    void                RemoveAll();

    /** Clean up. m_aopData[i]->UTILSId() == i after this call.
         Also deletes unused elements */
    void                Renumber( Array<int> &out_aiMap );

    ///
    WINbool             operator==( const UTILSElementArray & )const ;

    /// See above for affect of allocating or not
    UTILSElementArray( UTILSElement::TypeEnum in_oType, WINbool in_bAllocate );
    /// Will actually delete elements if in_bAllocate == TRUE
    ~UTILSElementArray( );

    /**@name Read, write */
    //@{
    ///
    void Read (ifstream &in, WINbool in_bForceAlloc);    
    ///
    void Write(ofstream &out) const;
    ///
    void ReadBinary (ifstream &in, WINbool in_bForceAlloc);    
    ///
    void WriteBinary(ofstream &out) const;
    //@}

    WINbool Check() const;
};
    
//@}

#include "Utils_Element_i.H"

#endif

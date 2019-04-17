// This may look like C code, but it is really -*- C++ -*-
// Copyright Microsoft Corporation.  Written by Hugues Hoppe.

#ifndef _DEFS_MESH_ARRAY_H
#define _DEFS_MESH_ARRAY_H

#include <WINSystemDefines.H>

/** \class SArray Mesh_Array.H utils/Mesh_Array.H

  \ingroup MeshesHelper

  \brief Array whose size is defined once.

    SArray<T> supports classes that have no public operator=() <br>

  This is essentially the stl vector class, but with a fixed, allocate once size.
  If you need to be able to re-allocate the size, use the Array class.

  Use: Either specify the size at creation:

  SArray array(num elems);

  or in a one-time init call

  SArray array;

  array.init(num elems)

  The only reason to use this one (vs Array) is if the class you are allocating does
  not have a copy constructor or operator= method.

File: include/Mesh_Array.H
*/
template<class T>
class SArray {
public:
    /// Initialize. If size is zero, can resize once.
    inline SArray(const int psize=0) : size(psize), a(size?new T[size]:0) { }
    ///
    inline ~SArray()                    { delete[] a; }

    /// One-time late allocation of array
    inline void init(const int s)             { ASSERT(!a); size=s; a=new T[size]; }
    /// Debug check
    inline void ok(const int i) const         { ASSERT(a && i>=0 && i<size); }
    /// Return a pointer to the first element of the list
    inline operator const T*() const    { return a; }
    /// Return a pointer to the first element of the list
    inline operator T*()                { return a; }
    /// Checks if index is ok
    inline const T& operator[](const int i) const     { ok(i); return a[i]; }
    /// Checks if index is ok
    inline       T& operator[](const int i)           { ok(i); return a[i]; }
    /// Number of elements in the list
    inline int num() const              { return size; }

protected:
    /// must come before a
    int size;                   
    /// Alloc'd array
    T* a;
private:
    inline SArray<T>& operator=(SArray<T>&) { return *this; } // disable
    inline SArray(SArray<T>&) { }         // disable
};

/** \class Array Mesh_Array.H utils/Mesh_Array.H
  \ingroup MeshesHelper
  \brief Basically the vector class in stl; a re-sizable list of elements

  For Array<T>, T must have a public operator=() and copy constructor. <br>

  It is more efficient, if you know the size of the array in advance, to allocate
  it at initialization rather than use Add() or operator += to add elements one
  at a time. Or use preallocate() to allocate approximately how much memory you will need.

  The class allocates slightly more memory than needed; this reduces the number
  of times malloc is called.

  Only use clearcompletely() if you really want all of the memory to go away, i.e.,
  you allocated something big and no longer need the memory. Otherwise, use need(0) to
  resize the array without unallocating the memory.

File: include/Mesh_Array.H
*/
template<class T>
class Array {
protected:
    /// How much extra to allocate
    const double resize_factor;
    /// must come before a; actual size of memory allocated
    int size;                   
    /// Data
    T* a;
    /// Size of array (n <= size)
    int n;
    
    /// Resize the arrays, if need be
    inline void resize(const int nsize);
    /// Resize the arrays, if need be
    inline void resize_aux(const int ne);
    
public: 
    /**@name To allow ARRAY_BIT_COPY(T). Don't use directly. */
    //@{
    /// Do not use
    inline void copy_raw_bits(T* na, const T* oa, const int n);
    /// Do not use
    inline void copy_array(T* na, const T* oa, const int n);
    //@}

 public:
    /**@name Constructors and assignments */
    //@{
    /// Default size is zero.
    inline Array(const int ne=0);
    ///
    inline Array(const Array<T>& ar);
    ///
    inline Array<T>& operator=(const Array<T>& ar);
    ///
    inline virtual ~Array()                     { delete[] a; }
    //@}

    /**@name Accessors */
    //@{
    /// Number of elements
    inline int  num() const              { return n; }
    /// Get direct access to the allocated memory. Use with caution.
    inline operator const T*() const    { ASSERT(a); return a; }
    /// Get direct access to the allocated memory. Use with caution.
    inline operator T*()                { ASSERT(a); return a; }
    /// Is the index 0 <= i < n? assert false if not
    inline void ok(const int i) const         
    { 
#if defined(DEBUG) || defined(_DEBUG)
        if ( a == NULL )
            TRACE("Empty\n");
        ASSERT(a); 
        if ( i < 0 || i >= n )
            TRACE("Overlow %d %d\n", i, n);
        ASSERT(i>=0 && i<n); 
#endif
    }
    /// Is the member[i] valid? false if NULL
    inline bool Ok(const int i) const         
    { 
#if defined(DEBUG) || defined(_DEBUG)
        if ( a == NULL )
            TRACE("Empty\n");
        ASSERT(a); 
        if ( i < 0 || i >= n )
            TRACE("Overlow %d %d\n", i, n);
        ASSERT(i>=0 && i<n); 
        if ( a[i] == NULL ){
            TRACE("Null member\n");
            return FALSE;
        }
#endif
        return TRUE;
    }
    /// To prevent compiler problems with .net
    inline const T& operator[](const size_t i) const { ok(i); return a[i]; }
    ///
    inline const T& operator[](const int i) const { ok(i); return a[i]; }
    ///
    inline       T& operator[](const int i)       { ok(i); return a[i]; }
    ///
    inline       T& operator[](const size_t i)    { ok(i); return a[i]; }
    /// operator[] with a wrap (doesn't work with greatly negative numbers)
    inline const T& wrap(const int i) const { ASSERT(a); return a[(i+n)%n]; }
    /// operator[] with a wrap (doesn't work with greatly negative numbers)
    inline      T& wrap(const int i)        { ASSERT(a); return a[(i+n)%n]; }
    /// operator[] with a clamp
    inline const T& clamp(const int i) const { return operator[]( (i<0) ? 0 : ((i>=n) ? n-1 : i) ); }
    /// operator[] with a clamp
    inline      T& clamp(const int i)        { return operator[]( (i<0) ? 0 : ((i>=n) ? n-1 : i) ); }
    ///
    inline const T& last() const      { ASSERT(a); return operator[](n-1); }
    ///
    inline       T& last()            { ASSERT(a); return operator[](n-1); }
    //@}

    /**@name Add and delete elements */
    //@{
    /// Really wipe; use need(0) if you are still using the Array later
    inline void clearcompletely()       { delete[] a,a=0; n=size=0; }
    /// Resizes memory storage if array is much smaller than array size
    inline void init(const int ne);
    /// Alloc at least e+1, COPY if too small
    inline void access(const int e);
    /// Allocate exactly ne space, COPY if too small
    inline void need(const int ne);  
    /// Allocate at least ne space, COPY if too small
    inline void needExpand(const int ne);  
    /// allocate at least ne space, COPY if too small, set new to input
    inline void needSet(const int ne, const T&);  
    /// set item ne, alloc and COPY if too small
    inline void set(const int ne, const T& );  
    /// set item ne
    inline void fill( const T& );  
    /// Use memset to fill the array
    inline void memfill( unsigned char in_oE ) ;
    /// add ne items, returns previous size
    inline int  add(const int ne);   
    /// subtract ne items off the end of the list
    inline void sub(const int ne);
    /// Add items to end
    inline void add( const Array<T> & ); 
    /// Add items, first item is at pts[iAt]
    inline void add( const Array<T> &, const int iAt ); 
    /// Add items to end
    inline void addItem( const T & ); 
    /// Add items, first item is at pts[iAt]
    inline void addItem( const T &, const int iAt ); 
    /// Delete items from iMe on, iN of them
    inline void del( const int iMe, const int iN = 1 );
    /// Adds an element. Do not use a+=a[..]!!
    inline Array<T>& operator+=(const T& e);
    /// Resize allocated space to be number of elements in array
    inline void squeeze()               { if (n<size) resize(n); }
    /// Make some space
    inline void preallocate(const int s);
    //@}

    /**@name Reorder items */
    //@{
    ///
    inline void swap(const int i, const int j)      { ok(i); ok(j); T e=a[i]; a[i] = a[j]; a[j]=e; }
    /// Reverses the elements in the list
    inline void reverse();
    /// Order items. Uses passed-in comparision function to determine order.
    inline void sort( int comp(const T&, const T&));
    /** Sort items by the values given.
        values should be the same size as the array.
        Sorts the values list, too.
        @param io_adValues assigns a value to each item. Is sorted along with array. */
    inline void sort(Array<double> &io_adValues); // sorts values, too
    ///
    inline WINbool unique() const;
    //@}

    /**@name Finding items in the list */
    //@{
    ///
    inline WINbool contains(const T& e) const;
    ///
    inline int index(const T& e, const int iAt = 0) const;
    //@}

    /**@name Queries about ordering of list */
    //@{
    ///
    inline int min_i(int comp(const T&, const T&)) const;
    ///
    inline int max_i(int comp(const T&, const T&)) const;
    ///
    inline const T& min_e(int comp(const T&, const T&)) const;
    ///
    inline const T& max_e(int comp(const T&, const T&)) const;

    /// Use the values in the in_adValues array to determine order
    inline int min_i(const Array<double> &in_adValues) const;
    /// Use the values in the in_adValues array to determine order
    inline int max_i(const Array<double> &in_adValues) const;
    /// Use the values in the in_adValues array to determine order
    inline const T&  min_e(const Array<double> &in_adValues) const;
    /// Use the values in the in_adValues array to determine order
    inline const T&  max_e(const Array<double> &in_adValues) const;
    //@}

    /**@name Read and write binary. Use at own risk - assumes class type T has no pointers */
    //@{
    /// Not valid if cross-platform (big vs little indian)
    inline void ReadBinary( ifstream &in ) {
        int iN;
        in.read( (char *) &iN, sizeof(int) );
        need(iN);
        in.read( (char *) a, iN * sizeof(T) );
    }
    /// Not valid if cross-platform (big vs little indian)
    inline void WriteBinary( ofstream &out ) const {
        out.write( (const char *) &n, sizeof(int) );
        out.write( (const char *) a, n * sizeof(T) );
    }

};

/** \brief Operations on all items in the list 
  \ingroup MeshesHelper

  It would be nice to specialize copy_array for all T of type 
   pointer but I cannot find a way to do this.

   template<class T*>
   void Array<T*>::copy_array(T** na, const T** oa, const int n) \
   { copy_raw_bits(na,oa,n); }

*/
  //@{
///Perform a mapping of items in the array using a function pointer.
///This should work just like Scheme's or Ruby's map function.
#define ARRAY_BIT_COPY(T) \
EmptyTemplate \
inline void Array<T>::copy_array(T* na, const T* oa, const int in_iN) \
{ copy_raw_bits(na,oa,in_iN); }

///
ARRAY_BIT_COPY(int)
///
ARRAY_BIT_COPY(float)
///
ARRAY_BIT_COPY(double)


//@}
//-------------------------------------------------------


#if defined(linux)
extern void *memcpy (void *__restrict __dest,
		     __const void *__restrict __src, size_t __n)
     __THROW __nonnull ((1, 2));

#endif

/* Copy N bytes of SRC to DEST.  */
template<class T>
inline
void Array<T>::copy_raw_bits(T* na, const T* oa, const int in_iN)
{
    // problem: does not call T::operator=()
  if (in_iN) memcpy((void*)na,(const void*)oa, (size_t) ( in_iN*sizeof(T) ));
}


template<class T>
inline void Array<T>::copy_array(T* na, const T* oa, const int in_iN)
{
  for (int i = 0; i < in_iN; i++) { na[i]=oa[i]; } 
}

template<class T>
inline
Array<T>& Array<T>::operator=(const Array<T>& ar)
{
    if (&ar==this) return *this;
    need(ar.num());
    copy_array(a,ar.a,n);
    return *this;
}


template<class T>
inline 
Array<T>::Array(const int ne) : resize_factor(1.5), a(NULL), n(0), size(0)
{ 
   size = ne;
   a = size ? new T[size] : 0;
   n = ne; 
}

template<class T>
inline
Array<T>::Array(const Array<T>& ar) : resize_factor(1.5), size(0), a(0), n(0)
{
    *this=ar;
}

template<class T>
inline
void Array<T>::resize(const int nsize)
{
    size=nsize;
    ASSERT( n<=size );
    T* na=size?new T[size]:0;
    copy_array(na,a,n);
    delete[] a,a=na;
}

template<class T>
inline
void Array<T>::resize_aux(const int ne)
{
    resize(WINminmax( int(n*resize_factor)+3, ne, 2 * ne ));
}

template<class T>
inline
void Array<T>::preallocate(const int s)
{
    if (n<s) resize(s);
    else ASSERT(n==s);
}

template<class T>
inline
WINbool Array<T>::contains(const T& e) const
{
    for (int i = 0; i < n; i++) { 
        if (a[i] == e) 
          return TRUE; 
    } 
    return FALSE;
}

template<class T>
inline
int Array<T>::index(const T& e, const int iAt) const
{
    for (int i = iAt; i < n; i++) { 
        if (a[i] == e) 
          return i; 
    } 
    return -1;
}

// allocate ne, CLEAR if too small
template<class T>
inline void Array<T>::init(const int ne) 
{  
    if (ne>size || ne < 1.5 * size) { 
        delete[] a;
        a = new T[ne];
        size = ne; 
    } 
    n = ne;
}

// allocate at least e+1, COPY if too small
template<class T>
inline void Array<T>::access(const int e) { 
    int ne=e+1;
    if (ne>size) resize_aux(ne);
    if (ne>n) n=ne;
}

// allocate ne, COPY if too small
template<class T>
inline void Array<T>::need(const int ne) {   
    if (ne>size) resize(ne);
    n=ne;
}

// allocate ne, COPY if too small
template<class T>
inline void Array<T>::needExpand(const int ne) {   
    if (ne>size) resize_aux(ne);
    n=ne;
}

// allocate ne, COPY if too small, set new to in
template<class T>
inline void Array<T>::needSet(const int ne, const T &in_oNew) 
{   
    int iN = num();
    if (ne>size) {
        resize_aux(ne);
    }
    n = WINmax( iN, ne );
    for (int i = iN; i < num(); i++)
        (*this)[i] = in_oNew;
}

// allocate ne, COPY if too small, set new to in,
// ne to last
template<class T>
inline void Array<T>::set(const int ne, const T &in_oE) 
{   
    needExpand( ne+1 );

    (*this)[ne] = in_oE;
}

// fill array with items
template<class T>
inline void Array<T>::fill(const T &in_oE) 
{   
    for ( int i = 0; i < num(); i++ )
        (*this)[i] = in_oE;
}

// fill array with items, using memcopy
template<class T>
inline void Array<T>::memfill( unsigned char in_oE ) 
{   
    memset( a, in_oE, sizeof( T ) * num() );
}

// ret: previous num()
template<class T>
inline int Array<T>::add(const int ne) 
{    
    int cn = n; 
    needExpand(n+ne); 
    return cn;
}

// ret: previous num()
template<class T>
inline void Array<T>::addItem(const T& e) 
{    
    // add() may allocate a[]!!
    int i=add(1); 
    a[i]=e; 
}

template<class T>
inline void Array<T>::sub(const int ne) 
{
    n-=ne;
}

// do not use a+=a[..]!!
template<class T>
inline Array<T>& Array<T>::operator+=(const T& e) 
{ 
    // add() may allocate a[]!!
    int i=add(1); 
    a[i]=e; 
    return *this; 
}

template<class T>
inline void Array<T>::reverse() 
{
    for(int i = 0; i < n/2; i++) { swap(i, n-1-i); } 
}

template<class T>
inline void Array<T>::sort(int comp(const T&, const T&))
{
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (comp(a[j], a[i]) < 0)
              swap(i, j);
        }
    }
}

template<class T>
inline
WINbool Array<T>::unique() const
{
    for (int i = 0; i < n; i++) {
        if (index(a[i],i+1) != -1)
            return FALSE;
    }
    return TRUE;
}

template<class T>
inline int Array<T>::min_i(int comp(const T&, const T&)) const
{
    int iMin = 0;
    for (int i = 0; i < n; i++) {
        if (comp(a[i], a[iMin]) < 0) iMin = i;
    }
    return iMin;
}

template<class T>
inline int Array<T>::max_i(int comp(const T&, const T&)) const
{
    int iMax = 0;
    for (int i = 0; i < n; i++) {
        if (comp(a[i], a[iMax]) > 0) iMax = i;
    }
    return iMax;
}

template<class T>
inline const T& Array<T>::min_e(int comp(const T&, const T&)) const
{
    return (*this)[min_i(comp)];
}


template<class T>
inline const T& Array<T>::max_e(int comp(const T&, const T&)) const
{
    return (*this)[max_i(comp)];
}


template<class T>
inline void Array<T>::sort(Array<double> &values)
{
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (values[j] < values[i]) {
                swap(i, j);
                values.swap(i, j);
            }
        }
    }
}

template<class T>
inline int Array<T>::min_i(const Array<double> &values) const
{
    int m_i = 0;
    for (int i = 0; i < n; i++) {
        if (values[i] < values[m_i]) m_i = i;
    }
    return m_i;
}

template<class T>
inline int Array<T>::max_i(const Array<double> &values) const
{
    int m_i = 0;
    for (int i = 0; i < n; i++) {
        if (values[i] > values[m_i]) m_i = i;
    }
    return m_i;
}

template<class T>
inline const T& Array<T>::max_e(const Array<double> &values) const
{
    return (*this)[max_i(values)];
}

template<class T>
inline const T& Array<T>::min_e(const Array<double> &values) const
{
    return (*this)[min_i(values)];
}

// Add points, first point at pts[iAt]
template<class T>
inline void Array<T>::add( const Array<T> &at, const int iAt )
{
    Array<T> aNew;

    int i;
    for ( i = 0; i < iAt; i++ )
        aNew += (*this)[i];

    for ( i = 0; i < at.num(); i++ )
        aNew += at[i];

    for ( i = iAt; i < num(); i++ )
        aNew += (*this)[i];

    (*this) = aNew;
}

// Add points, first point at pts[iAt]
template<class T>
inline void Array<T>::addItem( const T &e, const int iAt )
{
    Array<T> aNew;

    int i;
    for ( i = 0; i < iAt; i++ )
        aNew += (*this)[i];

    aNew += e;

    for ( i = iAt; i < num(); i++ )
        aNew += (*this)[i];

    (*this) = aNew;
}

// Add points, first point at pts[iAt]
template<class T>
inline void Array<T>::add( const Array<T> &at )
{
    for ( int i = 0; i < at.num(); i++ )
        (*this) += at[i];
}

//Delete items from iMe on, iN of em
template<class T>
inline void Array<T>::del( const int iMe, const int iN )
{
    int iStart = iMe;
    int iOff = iN;

    if ( (iMe + iN) - num() > 0 ) {
        iStart = 0;
        iOff = iN - (num() - iMe);
    }
    for ( int i = iStart; i < num(); i++ ) {
        if ( iOff + i < num() ) {
            (*this)[i] = (*this)[iOff + i];
        } else {
            break;
        }
    }

    sub( iN );
}


#endif


// -*-c++-*-
#ifndef IG_MATRIX_H
#define IG_MATRIX_H

// $Id: Rn_Matrix.H,v 1.17 2011/03/11 16:55:07 rsowell Exp $

// CwMtx matrix and vector math library
// Copyright (C) 1999-2001  Harry Kuiper
// Copyright (C) 2000  Will DeVore(template conversion)

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

// 5/25/2002 - Made changes to default includes --Cindy Grimm

#include <WINSystemDefines.H>


// CWMatrix class

// This library was designed to mirror as closely as possible the
// notation used in mathematical writing. A matrix is indexed:
// matMatrixName[row][col].

// CAUTION!!!

// This matrix library was implemented with emphasis on
// speed. Consequently no attempts were made to trap and report
// errors. It is left entirely to the user to write code that does not
// cause errors within the matrix routines.

namespace CwMtx
{
  // matrix status values
  enum { N_NOTALLOCATED, N_ALLOCATED, N_MAPPED };

  // Classes that create unity and zero "objects".  The ones directly
  // below work only for template arguments that are basic (numerical)
  // types that can be initialised by a literal 0 or 1.
  //template <class T> class CWTUnity
  //{
  //public:
  //  operator T() { return 1; }
  //};

  //template <class T> class CWTZero
  //{
  //public:
  //  operator T() { return 0; }
  //};

  // This template defaults to double. Most of the time this template
  // will be working with math functions that only work with
  // doubles. For example, the transcendental function sin(x) takes
  // and returns a double which would force the compiler to convert
  // back and forth from some other data type to a double.

  // prefix mat
  template < class T >
  class CWTMatrix
  {
  public:
    typedef T element;

    // creates a matrix, does NOT allocate rows and columns
    CWTMatrix();
    // creates a matrix, allocates rows and columns
    CWTMatrix(unsigned, unsigned);
    CWTMatrix(const CWTMatrix &);
    // sub-matrix mapped into another
    CWTMatrix(const CWTMatrix &, unsigned, unsigned, unsigned, unsigned);
	
    // removes matrix elements from free store
    ~CWTMatrix() { Deallocate(); };

    // allocates rows and colums
    void Dimension(unsigned, unsigned);
    // maps matrix into another
    void MapInto(const CWTMatrix&, unsigned, unsigned, unsigned, unsigned);
    // reverses the effect of Dimension() and MapInto()
    void Deallocate();

    int GetStatus() const { return m_nMatStatus; };
    unsigned GetRows() const { return m_crow; };
    unsigned GetCols() const { return m_ccol; };

    // basic matrix operations

    // returns a row of modifyable elements
    T* operator [](unsigned irow) { return m_rgrow[irow]; };
    // returns a row of non-modifyable elements
    const T* operator [](unsigned irow) const { return m_rgrow[irow]; };

    CWTMatrix operator +(const CWTMatrix &) const;
    CWTMatrix operator -(const CWTMatrix &) const;
    CWTMatrix operator -() const;
    CWTMatrix operator *(const T &) const;
    CWTMatrix operator *(const CWTMatrix &) const;

    // Interesting note here. Because we are defining the "/" operator
    // we can't expect to use it inside the definition unless we
    // safely restrict it to operating on constants. Without the
    // parens the operator does the "* 1" first then does the
    // "/value" which then leads to calling the "/" etc...  With the
    // parens, the intended scalar operation, "1/value", occurs
    // first then the Matrix/Scalar operations occur.  If the parens
    // are missing in the operator below, an ifinite loop
    // occurs. -----------------------------------------------V----------V
    CWTMatrix operator /(const T &value) const
    {
      return (*this)*static_cast<const T &>(1/value);
    }

    // not inherited
    CWTMatrix & operator =(const CWTMatrix &);
    CWTMatrix & operator +=(const CWTMatrix &);
    CWTMatrix & operator -=(const CWTMatrix &);
    CWTMatrix & operator *=(const T &);
    CWTMatrix & operator /=(const T &value)
    {
      return (*this) *= static_cast<const T &>(1/value);
    }

    int operator ==(const CWTMatrix &) const;
    int operator !=(const CWTMatrix &mat) const { return !( (*this) == mat ); }

    // stores CWMatrix + CWMatrix in this
    void StoreSum(const CWTMatrix &, const CWTMatrix &);
    // stores CWMatrix*CWMatrix in this
    void StoreProduct(const CWTMatrix &, const CWTMatrix &);
    // stores transpose of CWMatrix in this
    void StoreTranspose(const CWTMatrix &);
    // stores transpose of CWMatrix in this
    void StoreInverse( const CWTMatrix<T> &, const bool do_not_check_ident = false );
    // stores CWMatrix at indicated position in this
    void StoreAtPosition(unsigned, unsigned, const CWTMatrix &);
    // fills the whole array with a value.
    void Fill(const T &);
    ///
    void SetIdentity();
    /// the matrix's determinant
    T Determinant();

    void InterchangeRows(unsigned, unsigned);
    void AddRowToRow(unsigned, unsigned, const T &);
    void MultiplyRow(unsigned, const T &);

    void Print() const;

  private:
    // initializes data members
    void Initialize();

    /// calculate the determinant of a subset of the matrix
    T Determinant( bool * useRows, bool * useCols, int nuse );

    // we keep the data structures used for CWMatrix implementation
    // private

    // row count
    unsigned int m_crow;
    // column count
    unsigned int m_ccol;
    // an array of rows (stored on free store)
    T **m_rgrow;
    // matrix status
    int m_nMatStatus;
  };

  // Templates to create self-dimensioning CWTMatrix classes - or one
  // of its derived classes - using the syntax of a default
  // constructor.  This facility is required for using matrices as
  // elements of matrices since these will always be created by a call
  // to the default constructor.

  template <class T, unsigned crow, unsigned ccol>
  class CWTMat: public T
  {
  public:
    CWTMat(): T(crow, ccol) {}

    T & operator =(const T &mtx) { return T::operator=(mtx); }
  };

  // NOTE: There exists no unity matrix for a non-square matrix!

  // Zero matrix.  NOTE: A zero matrix can only be constructed for a
  // matrix of known dimensions.  Hence the use of CWTMat<T,n,m>.
  //template <class T, unsigned crow, unsigned ccol>
  //class CWTZero< CWTMat< CWTMatrix<T>, crow, ccol> > :
  //  public CWTMat< CWTMatrix<T>, crow, ccol>
  //{
  //public:
  //  CWTZero() { Fill(CWTZero<T>()); }
  //};

  //
  // Constructors
  //

  template < class T >
  inline CWTMatrix<T>::CWTMatrix()
  {
    Initialize();
  }

  // creates a matrix, allocates rows and columns
  template < class T >
  inline CWTMatrix<T>::CWTMatrix(unsigned crow, unsigned ccol)
  {
    Initialize();
    Dimension(crow, ccol);
  }

  template < class T >
  inline CWTMatrix<T>::CWTMatrix(const CWTMatrix<T> &mat)
  {
    Initialize();

    if (mat.m_nMatStatus == N_NOTALLOCATED)
      {
	// input matrix not allocated, so there's nothing to copy
	return;
      }
    else
      {
	// copy contents of input matrix
	(*this) = mat;
      }
  }

  // Mapped matrix constructor
  template < class T >
  inline CWTMatrix<T>::CWTMatrix(const CWTMatrix<T> &mat,
				 unsigned irowStart,
				 unsigned icolStart,
				 unsigned irowEnd,
				 unsigned icolEnd)
  {
    Initialize();
    MapInto(mat, irowStart, icolStart, irowEnd, icolEnd);
  }

  //
  // Private Methods
  //

  // initial values for matrix attributes
  template < class T >
  inline void CWTMatrix<T>::Initialize()
  {
    m_crow       = 0;
    m_ccol       = 0;
    m_rgrow      = NULL;
    m_nMatStatus = N_NOTALLOCATED;
  }

  //
  // User Methods
  //

  template < class T >
  void CWTMatrix<T>::Dimension(unsigned crowInit, unsigned ccolInit)
  {
    if (m_nMatStatus != N_NOTALLOCATED)
      {
	Deallocate();
      }

    m_crow = crowInit;
    m_ccol = ccolInit;

#ifdef CC_CWTMTX_ASSUME_BASIC_TYPES
    // NOTE: CWTMatrix normally uses the standard C++ new() operator
    // to allocate memory for matrix elements.  Using malloc(3) can
    // save time because it can allocate all required memory in a
    // single call.  However, malloc(3) only works for template
    // arguments that are C++ basic types because basic types are no
    // classes and thus need not be initialised by a constructor.  To
    // use malloc(3) #define CC_CWTMTX_ASSUME_BASIC_TYPES in all
    // source files that use CWTMatrix templates before you #include
    // them.

    // Allocate space for row pointers and the rows themselves using
    // ANSI C malloc(3) function.
    m_rgrow = reinterpret_cast<T **>(malloc(m_crow*sizeof(T *)
					    + m_crow*m_ccol*sizeof(T)));
    T *ptTmp = reinterpret_cast<T *>(&(m_rgrow[m_crow]));
#else
    // Allocate space for row pointers and the rows themselves using
    // ANSI C++ new() operator.
    m_rgrow = new T*[m_crow];
    T *ptTmp = new T[m_crow*m_ccol];
#endif

    // make row pointers point to start of each row
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	m_rgrow[irow] = &(ptTmp[irow*m_ccol]);
      }

    SetIdentity();
    m_nMatStatus = N_ALLOCATED;
  }

  // maps a matrix into another matrix, deallocates first if neccesary
  // allocates space on the free store for a matrix, deallocates first
  // if neccesary
  template < class T >
  void CWTMatrix<T>::MapInto(const CWTMatrix<T> &mat,
			     unsigned irowStart,
			     unsigned icolStart,
			     unsigned irowEnd,
			     unsigned icolEnd )
  {
    if (m_nMatStatus != N_NOTALLOCATED)
      {
	Deallocate();
      }

    // calculate columns
    m_crow = irowEnd - irowStart + 1;

    // calculate rows
    m_ccol = icolEnd - icolStart + 1;

    // allocate space for row pointers
#ifdef CC_CWTMTX_ASSUME_BASIC_TYPES
    m_rgrow = reinterpret_cast<T **>(malloc(m_crow*sizeof(T *)));
#else
    m_rgrow = new T*[m_crow];
#endif

    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	// get values for row pointers
	m_rgrow[irow] = &mat.m_rgrow[irow + irowStart][icolStart];
      }

    m_nMatStatus = N_MAPPED;
  }

  // deallocates a matrix' elements from the free store
  template < class T >
  void CWTMatrix<T>::Deallocate()
  {
    // has to take account of the state the matrix is in
    switch (m_nMatStatus)
      {
      case N_NOTALLOCATED:
	// nothing needs to be deallocated
	break;

      case N_MAPPED:
	// delete the array of row pointers
#ifdef CC_CWTMTX_ASSUME_BASIC_TYPES
	free(m_rgrow);
#else
	delete [] m_rgrow;
#endif
	break;

      case N_ALLOCATED:
	//   delete the contents of the matrix completely
#ifdef CC_CWTMTX_ASSUME_BASIC_TYPES
	// Deallocate space for row pointers and the rows themselves
	// using ANSI C free(3) function.
	free(m_rgrow);
#else
	// Deallocate space for row pointers and the rows themselves
	// using ANSI C++ delete() operator.
	delete [] *m_rgrow;
	delete [] m_rgrow;
#endif
	break;
      };

    // bring this matrix in initialized state
    Initialize();
  }

  template < class T >
  CWTMatrix<T> 
  CWTMatrix<T>::operator +(const CWTMatrix<T> &mat) const
  {
    // copy this and add argument
    return CWTMatrix( *this ) += mat;
  }

  template < class T >
  CWTMatrix<T> 
  CWTMatrix<T>::operator -(const CWTMatrix<T> &mat) const
  {
    // copy this and subtract argument
    return CWTMatrix( *this ) -= mat;
  }

  template < class T >
  CWTMatrix<T> CWTMatrix<T>::operator -() const
  {
    return (*this)*static_cast<const T &>(0 - 1);
  }

  template < class T >
  CWTMatrix<T> CWTMatrix<T>::operator *(const T &value) const
  {
    // copy this and multiply by argument
    return CWTMatrix( *this ) *= value;
  }

  template < class T >
  CWTMatrix<T>
  CWTMatrix<T>::operator *(const CWTMatrix<T> &mat) const
  {
    // create result matrix
    CWTMatrix matResult( m_crow , mat.m_ccol );
    // store (*this)*mat in result matrix
    matResult.StoreProduct( *this , mat );
    return matResult;
  }

  // assignment operator
  template < class T >
  CWTMatrix<T> & 
  CWTMatrix<T>::operator =(const CWTMatrix<T> &mat)
  {
    if (m_nMatStatus == N_NOTALLOCATED)
      {
	// if this is not allocated, we'll allocate it on the fly
	Dimension(mat.m_crow, mat.m_ccol);
      }

    // Copy the values from mat to this, excess values in mat are
    // ignored. If mat has too few elements, garbage will be copied
    // into remaining elements of this
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] = mat.m_rgrow[irow][icol];
	  }
      }

    return *this;
  }

  template < class T >
  CWTMatrix<T> & 
  CWTMatrix<T>::operator +=(const CWTMatrix<T> &mat)
  {
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] += mat.m_rgrow[irow][icol];
	  }
      }

    return *this;
  }

  template < class T >
  CWTMatrix<T> & 
  CWTMatrix<T>::operator -=(const CWTMatrix<T> &mat)
  {
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] -= mat.m_rgrow[irow][icol];
	  }
      }

    return *this;
  }

  template < class T >
  CWTMatrix<T> & CWTMatrix<T>::operator *=(const T &value)
  {
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] *= value;
	  }
      }

    return *this;
  }

  template < class T >
  int CWTMatrix<T>::operator ==(const CWTMatrix<T> &mat) const
  {
    if ((m_crow == mat.m_crow) && (m_ccol == mat.m_ccol))
      {
	for (unsigned irow = 0; irow < m_crow; ++irow)
	  {
	    for (unsigned icol = 0; icol < m_ccol; ++icol)
	      {
		if (m_rgrow[irow][icol] != mat.m_rgrow[irow][icol])
		  {
		    return 0;
		  }
	      }
	  }

	return 1;
      }
    else
      {
	return 0;
      }
  }

  // stores mat1 + mat2 in this
  template < class T >
  void CWTMatrix<T>::StoreSum(const CWTMatrix<T> &mat1,
			      const CWTMatrix<T> &mat2)
  {
    // NOTE: it is assumed that this has correct dimensions,
    // i.e. CWMatrix(mat1.m_crow, mat2.m_ccol)

      ASSERT( mat1.GetCols() == mat2.GetCols() );
      ASSERT( mat1.GetRows() == mat2.GetRows() );
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] = 
	      mat1.m_rgrow[irow][icol] + mat2.m_rgrow[irow][icol];
	  }
      }
  }

  // stores mat1*mat2 in this
  template < class T >
  void CWTMatrix<T>::StoreProduct(const CWTMatrix<T> &mat1,
				  const CWTMatrix<T> &mat2)
  {
    // NOTE: it is assumed that this has correct dimensions,
    // i.e. CWMatrix(mat1.m_crow, mat2.m_ccol)

      ASSERT( mat1.GetCols() == mat2.GetRows() );
    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] = 0.0;

	    for (unsigned icol2 = 0; icol2 < mat1.m_ccol; ++icol2)
	      {
		m_rgrow[irow][icol] += 
		  mat1.m_rgrow[irow][icol2]*mat2.m_rgrow[icol2][icol];
	      }
	  }
      }
  }

  // mat should fit in this
  template < class T >
  void CWTMatrix<T>::StoreAtPosition(unsigned irowStart,
				     unsigned icolStart,
				     const CWTMatrix<T> &mat)
  {
    for (unsigned irow = 0; irow < mat.m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < mat.m_ccol; ++icol)
	  {
	    m_rgrow[irow + irowStart][icol + icolStart] =
	      mat.m_rgrow[irow][icol];
	  }
      }
  }

  // stores transpose of mat in this
  template < class T >
  void CWTMatrix<T>::StoreTranspose(const CWTMatrix<T> &mat)
  {
    // NOTE: it is assumed that this has correct dimensions,
    // i.e. CWMatrix(mat.m_ccol, mat.m_crow)

    for (unsigned irow = 0; irow < m_crow; ++irow)
      {
	for (unsigned icol = 0; icol < m_ccol; ++icol)
	  {
	    m_rgrow[irow][icol] = mat.m_rgrow[icol][irow];
	  }
      }
  }

  template < class T >
  inline void CWTMatrix<T>::Fill(const T &elemFill)
  {
    unsigned iEnd = m_crow*m_ccol;

    for (unsigned i = 0; i < iEnd; ++i)
      {
	(*m_rgrow)[i] = elemFill;
      }
  }

  template < class T >
  inline void CWTMatrix<T>::SetIdentity( )
  {
    for ( unsigned int iR = 0; iR < m_crow; iR++ ) {
        memset( m_rgrow[iR], 0, sizeof(T) * m_ccol );
        if ( iR < m_ccol )
            (*this)[iR][iR] = 1.0;
    }    
  }

  //for now, this only works for square matrices
  template < class T >
  inline T CWTMatrix<T>::Determinant( )
  {
      //return nan for non-square matrices
      if( m_crow != m_ccol ) 
          //return nan(""); //I hope this gives you NaN; it gives me NaN
          return sqrt(-1.0); //I hope this gives you NaN; it gives me NaN

      double det=1;
      //make the matrix upper-triangular via the Gauss-Jordan method
      CWTMatrix<T> workmat(*this);

      for (unsigned int iCol = 0; iCol < m_ccol; ++iCol) {
          int iRowSwap = -1;
          for ( unsigned int i = iCol; i < m_crow; i++ ) {
              if ( iRowSwap==-1 || 
                    fabs(workmat[i][iCol])>fabs(workmat[iRowSwap][iCol])) {
                  iRowSwap = i;
                  break;
              }
          }
          //if the entire column is zero, the determinant is 0
          if ( iRowSwap == -1 ) {
              return 0;
          }

          //swap rows
          if( iRowSwap != iCol ){
              T tt;
              for( unsigned int tc=0; tc<m_ccol; ++tc ){
                  //swap mat
                  tt = workmat[iRowSwap][tc];
                  workmat[iRowSwap][tc] = workmat[iCol][tc];
                  workmat[iCol][tc] = tt;
              }
          }
  
          //scale row
          const double dScl = 1.0 / workmat[iCol][iCol];
          workmat.MultiplyRow( iCol, dScl );
          det /= dScl;
          if( !RNApproxEqual( workmat[iCol][iCol], 1.0, 1e-6 ) )
              cout << "determinant: triangularization failure at " << iCol << ", " << iCol << " -> " << workmat[iCol][iCol] << endl;
          
          //subtract row from all others
          for ( unsigned int iR = iCol+1; iR < m_crow; iR++ ) {
              const double dAdd = -workmat[iR][iCol];
              workmat.AddRowToRow( iCol, iR, dAdd );
              if( !RNApproxEqual( workmat[iR][iCol], 0.0, 1e-6 ) ){
                  //if dAdd is very big, then rounding errors can keep the matrix element from going to 0.  Therefore, we should repeat in that case
                  if( fabs(dAdd) > 1e3 ){
                      cout << "failure at " << iCol << " " << iR << " -> " << workmat[iR][iCol] << "  -- repeating!" << endl;
                      iR--;
                      continue;
                  }
                  //however, if that wasn't the case, we have a problem
                  ASSERT(FALSE);
              }
          }
      }
      //determinant code -- because it's upper-triangular, only the
      //elements on the diagonal contribute to the determinant
      for( unsigned int i=0; i<m_crow; i++ )
          det *= workmat[i][i];
      return det;
  }

  template < class T >
  void CWTMatrix<T>::InterchangeRows(unsigned irow1, unsigned irow2)
  {
    T* rowSav      = m_rgrow[irow1];

    // switch row pointers
    m_rgrow[irow1] = m_rgrow[irow2];
    m_rgrow[irow2] = rowSav;
  }

  template < class T >
  void CWTMatrix<T>::MultiplyRow(unsigned irow, const T &value)
  {
    for (unsigned icol = 0; icol < m_ccol; ++icol)
      {
	m_rgrow[irow][icol] *= value;
      }
  }

  template < class T >
  void CWTMatrix<T>::AddRowToRow(unsigned irowSrc,
				 unsigned irowDest,
				 const T &value)
  {
    for (unsigned icol = 0; icol < m_ccol; ++icol)
      {
	m_rgrow[irowDest][icol] += m_rgrow[irowSrc][icol]*value;
      }
  }

  //
  // template functions designed work with the base matrix class.
  //

  template < class T >
  inline CWTMatrix<T> operator *(const T &value,
				 const CWTMatrix<T> &mat)
  {
    return mat*value;
  }

  template < class T >
  CWTMatrix<T> transpose(const CWTMatrix<T> &mat)
  {
    CWTMatrix<T> matTranspose( mat.GetCols(), mat.GetRows() );
    matTranspose.StoreTranspose(mat);
    return matTranspose;
  }

  template < class T >
  CWTMatrix<T> inverse( const CWTMatrix<T> &mat )
  {
      CWTMatrix<T> matInverse( mat.GetRows(), mat.GetCols() );
      matInverse.StoreInverse( mat);
      return matInverse;
  }
  template < class T >
  inline ostream & operator <<(ostream &os, const CWTMatrix<T>& mtx)
  {
    os << "{" ;

    for (unsigned i = 0; i < mtx.GetRows(); i++)
      {
	if (i != 0)
	  {
	    cout << ",";
	  }

	for (unsigned j = 0; j < mtx.GetCols(); j++)
	  {
	    os  << " " << mtx[i][j];
	  }
      }

    os << " }";

    return os;
  }

  template < class T >
  inline void CWTMatrix<T>::Print() const
  {
      for ( unsigned irow = 0; irow < GetRows(); irow++ ) {
          for (unsigned icol = 0; icol < GetCols(); ++icol) {
	          TRACE("%f ", (*this)[irow][icol]);
          }
          TRACE("\n");
      }
      TRACE("\n");
  }

  typedef CWTMatrix<double> RNMatrix;
}


#endif // IG_MATRIX_H


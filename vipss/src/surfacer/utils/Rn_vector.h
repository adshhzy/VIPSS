// -*-c++-*-
#ifndef IG_VECTOR_H
#define IG_VECTOR_H

// $Id: Rn_vector.h,v 1.2 2004/05/07 20:24:51 cmg Exp $

// CwMtx matrix and vector math library
// Copyright (C) 1999-2001  Harry Kuiper
// Copyright (C) 2000  Will DeVore (template conversion)

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
#include <WINSystemDefines.h>

#ifndef IG_MATRIX_H
#include <utils/Rn_Matrix.H>
#endif

namespace CwMtx
{

  template < class T = double >
  class CWTVector : public CWTMatrix<T>
  {
  public:
    CWTVector(): CWTMatrix<T>() {};
    CWTVector(unsigned crowInit): CWTMatrix<T>(crowInit, 1) {};
    CWTVector(const CWTMatrix<T>& mat): CWTMatrix<T>(mat) {};
    CWTVector(const CWTVector& vec): CWTMatrix<T>(vec) {};
    // mapping into a matrix
    CWTVector(const CWTMatrix<T>& , unsigned, unsigned, unsigned);
    // mapping into a vector
    CWTVector(const CWTVector& , unsigned, unsigned);

    ~CWTVector() {};

    void MapInto(const CWTMatrix<T> &, unsigned, unsigned, unsigned);
    void MapInto(const CWTVector &, unsigned, unsigned);
    void Dimension(unsigned crowInit)
    {
      CWTMatrix<T>::Dimension(crowInit, 1);
    }

    T & operator [](unsigned irow)
    {
      return this->CWTMatrix<T>::operator[](irow)[0];
    }

    const T & operator [](unsigned irow) const
    {
      return this->CWTMatrix<T>::operator[](irow)[0];
    }

    CWTVector operator +(const CWTVector &) const;
    CWTVector operator -(const CWTVector &) const;
    CWTVector operator -() const;
    CWTVector operator *(const T &) const;
    // CWTVector*CWTVector, inner product
    T operator *(const CWTVector &) const;
    CWTVector operator /(const T &value) const
    {
      return (*this)*static_cast<const T &>(CWTUnity<T>()/value);
    }

    // not inherited
    CWTVector & operator =(const CWTVector &vec);
    CWTVector & operator +=(const CWTVector &vec);
    CWTVector & operator -=(const CWTVector &vec);
    CWTVector & operator *=(const T &value);
    CWTVector & operator /=(const T &value);

    // CWTVector norm
    T operator !() const { return (*this).Norm(); };

    void StoreAtRow(unsigned, const CWTVector &);
    // returns vector norm (length)
    T Norm() const;
    // returns a unit vector with same direction as this
    CWTVector Unit() const { return (*this)/Norm(); }

    // make this a unit vector
    void MakeUnit() { (*this) /= Norm(); }
  };

  template <class T, unsigned crow>
  class CWTVec: public T
  {
  public:
    CWTVec(): T(crow) {}

    T & operator =(const T &mtx) { return T::operator=(mtx); }
  };

  // NOTE: There exists no unity vector for a general vector!

  // Zero matrix.
  //template <class T, unsigned crow>
  //class CWTZero< CWTVec<CWTVector<T>, crow> >:
  //  public CWTVec<CWTVector<T>, crow>
  //{
  //public:
  //  CWTZero() { Fill(CWTZero<T>()); }
  //};

  //
  // Constructors
  //

  // mapping into a vector
  template < class T >
  inline CWTVector<T>::CWTVector(const CWTVector<T> &vec,
				 unsigned irowStart,
				 unsigned irowEnd)
    : CWTMatrix<T>(vec, irowStart, 0, irowEnd, 0)
  {
  }

  // mapping into a matrix
  template < class T >
  inline CWTVector<T>::CWTVector(const CWTMatrix<T> &mat,
				 unsigned irowStart,
				 unsigned icolStart,
				 unsigned irowEnd)
    :
    CWTMatrix<T>(mat, irowStart, icolStart, irowEnd, icolStart)
  {
  }

  //
  // User Methods
  //

  template < class T >
  inline void CWTVector<T>::MapInto(const CWTMatrix<T> &mat,
				    unsigned irowStart,
				    unsigned icol,
				    unsigned irowEnd)
  {
    CWTMatrix<T>::MapInto(mat, irowStart, icol, irowEnd, icol);
  }

  template < class T >
  inline void CWTVector<T>::MapInto(const CWTVector &vec,
				    unsigned irowStart,
				    unsigned irowEnd)
  {
    CWTMatrix<T>::MapInto(vec, irowStart, 0, irowEnd, 0);
  }

  // not inherited
  template < class T >
  inline CWTVector<T> & CWTVector<T>::operator =(const CWTVector<T> &vec)
  {
    return static_cast<CWTVector &>(CWTMatrix<T>::operator=(vec));
  }

  template < class T >
  inline CWTVector<T> & CWTVector<T>::operator +=(const CWTVector<T> &vec)
  {
    return static_cast<CWTVector &>(CWTMatrix<T>::operator+=(vec));
  }

  template < class T >
  inline CWTVector<T> & CWTVector<T>::operator -=(const CWTVector<T> &vec)
  {
    return static_cast<CWTVector &>(CWTMatrix<T>::operator-=(vec));
  }

  template < class T >
  inline CWTVector<T> & CWTVector<T>::operator *=(const T &value)
  {
    return static_cast<CWTVector &>(CWTMatrix<T>::operator*=(value));
  }

  template < class T >
  inline CWTVector<T> & CWTVector<T>::operator /=(const T &value)
  {
    return (*this) *= static_cast<const T &>(CWTUnity<T>()/value);
  }

  template < class T >
  inline void CWTVector<T>::StoreAtRow(unsigned irowStart,
				       const CWTVector<T> &vec)
  {
    CWTMatrix<T>::StoreAtPosition(irowStart, 0, vec);
  }

  template < class T >
  CWTVector<T> CWTVector<T>::operator +(const CWTVector<T> &vec) const
  {
    return CWTVector<T>(*this) += vec;
  }

  template < class T >
  CWTVector<T> CWTVector<T>::operator -(const CWTVector<T> &vec) const
  {
    return CWTVector<T>(*this) -= vec;
  }

  template < class T >
  CWTVector<T> CWTVector<T>::operator -() const
  {
    return (*this)*static_cast<const T &>(CWTZero<T>() - CWTUnity<T>());
  }

  template < class T >
  CWTVector<T> CWTVector<T>::operator *(const T &value) const
  {
    return CWTVector<T>(*this) *= value;
  }

  template < class T >
  T CWTVector<T>::operator *(const CWTVector<T> &vec) const
  {
    T elemResult = 0;

    for (unsigned irow = 0; irow < (*this).GetRows(); ++irow)
      {
	     elemResult += (*this)[irow]*vec[irow];
      }

    return elemResult;
  }

  // length of vector
  template < class T >
  T CWTVector<T>::Norm() const
  {
    T elemResult = 0;

    elemResult = (*this)*(*this);

    return sqrt( elemResult );
  }

  //
  // template functions designed work with the vector class.
  //

  template < class T >
  inline CWTVector<T> operator *(const T &value, const CWTVector<T> &vec)
  {
    return vec*value;
  }

  // matrix*vector must yield a vector
  template < class T >
  CWTVector<T> operator *(const CWTMatrix<T> &mat, const CWTVector<T> &vec)
  {
    CWTVector<T> vecResult(vec.GetRows());
    vecResult.StoreProduct(mat, vec);
    return vecResult;
  }

  // norm computation as a function
  template < class T >
  inline T Length(const CWTVector<T> &vec)
  {
    return vec.Norm();
  }

  typedef CWTVector<double> RNVector;

}


#endif // IG_VECTOR_H


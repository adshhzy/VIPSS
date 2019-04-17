// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef HashStruct_h
#define HashStruct_h

#if 0
{
	HashStruct hs<segment>;
	ForHashStruct(hs,segment,s) { do(s); } EndFor;
}
#endif

#include <utils/Mesh_Map.H>


/** \class BHashStruct Mesh_HashStruct.H utils/Mesh_HashStruct.H
   \ingroup MeshesHash

  \brief Basic hash table, user supplied has function

  This should be a templated class, but uses Univ (void *) instead.

 Maintains a set of structures.
 The structures are hashed according to a user-provided function hashf.
  Note: hashf can return any value including zero.
 Also, equality of two structures is determined by user function cmpf.
  Note: cmpf should return: does_not_match (as does strcmp)!
  Note for hacker: when cmpf is invoked, its first arg is st (below)
   Thus, st (structure template) can be fake and cmpf can be asymmetric.
 Note: BHashStruct should be empty prior to destruction.
*/
class BHashStruct : public BMap {
  public:
      /// Hash function
	typedef Univ (*HASHF)(Univ);
    /// Comparison function
	typedef int (*CMPF)(Univ, Univ);
    /// Need both a hash function and comparison function
	BHashStruct(HASHF phashf, CMPF pcmpf);
    ///
	~BHashStruct();
    /// s!=0, s must be new/not in table
	void enter(Univ s);	
    /// return: is_new y/n
	int add(Univ s);	
    /// return s found or 0
	Univ retrieve(Univ st) const; 
    /// return wasfound
	int remove(Univ st);
    /// die if empty, used for iteration ??
	Univ getone() const;	      
    /// die if empty
	Univ removeone();
  private:
	HASHF hashf;
	CMPF cmpf;
	Univ replace(Univ k, Univ v); // disable ?? permitted
};

/** \class BHashStructIter Mesh_HashStruct.H utils/Mesh_HashStruct.H
   \ingroup MeshesHash

  \brief Basic hash tableiterator
  */
class BHashStructIter : public BMapIter {
 public:
    inline BHashStructIter(const BHashStruct& s) : BMapIter(s) { }
    inline BHashStructIter(const BHashStruct& s, MeshRandom& r)
        : BMapIter(s,r) { }
};

//----------------------------------------------------------------------------

inline Univ BHashStruct::getone() const {
	BMapIter si(*this); return si.value();
}
inline Univ BHashStruct::removeone() {
	Univ e=getone(); remove(e); return e;
}

//----------------------------------------------------------------------------

template<class T> class HashStructIter;

/** \class HashStruct Mesh_HashStruct.H utils/Mesh_HashStruct.H
   \ingroup MeshesHash

  \brief Hash table. Template wrapper around HashStruct.

  See BHashStruct for method descriptions
  GNUG chokes on Conv<T*>, so I did it manually here
  */
template<class T>
class HashStruct : public BHashStruct {
  public:
      ///
	typedef Univ (*HASHF)(const T*);
      ///
	typedef int (*CMPF)(const T*, const T*);
      ///
	HashStruct(HASHF phashf, CMPF pcmpf) :
      ///
	BHashStruct((BHashStruct::HASHF)phashf,(BHashStruct::CMPF)pcmpf) { }
      ///
	~HashStruct() { }
      ///
	inline void enter(T* e) { BHashStruct::enter(Univ(e)); }
      ///
	inline T* retrieve(const T* e) const
	{ return (T*)(void*)(BHashStruct::retrieve(Conv<T*>::e((T*)e))); }
      ///
	inline int add(T* e) { return BHashStruct::add(Univ(e)); }
      ///
	inline int remove(const T* e)
	{ return BHashStruct::remove(Conv<T*>::e((T*)e)); }
      ///
	inline T* getone() const { return (T*)(void*)(BHashStruct::getone()); }
      ///
	inline T* removeone() { return (T*)(void*)(BHashStruct::removeone()); }
// 	typedef HashStructIter<T> Iter;
};

/** \class HashStructIter Mesh_HashStruct.H utils/Mesh_HashStruct.H
   \ingroup MeshesHash

  \brief Template wrapper around Basic hash tableiterator
  */
template<class T>
class HashStructIter : public BMapIter {
  public:
	HashStructIter(const HashStruct<T>& s) : BMapIter(s) { }
	HashStructIter(const HashStruct<T>& s, MeshRandom& r) : BMapIter(s,r) { }
	~HashStructIter() { }
	inline T* operator()() const { return (T*)(void*)(BMapIter::value()); }
  private:
	void key() const;	// disable
	void value() const;	// disable
};

#define ForHashStruct(S,T,V) \
{ for (HashStructIter<T> zz(S);zz;zz.next()) { T* V=zz();
					       
#endif


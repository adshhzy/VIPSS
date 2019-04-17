// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Hash_h
#define Hash_h

#if 0
{
	Map<Edge,Vertex> mev;
	ForMapKeyValue(mev,Edge,e,Vertex,v) { do(e,v); } EndFor;
	ForMapKey(mev,Edge,Vertex,e) { do(e); } EndFor;
	ForMapValue(mev,Edge,Vertex,v) { do(v); } EndFor;
}
#endif

#include <utils/Mesh_Hh.H>
#include <utils/Mesh_Pool.H>


/** \class BMap Mesh_Map.H utils/Mesh_Map.H
   \ingroup MeshesHash

  \brief Mapping between elements. Essentially a hash table. 

    This class forms the base class for the templated one below.
    */
class BMap {
  public:
	BMap();
	~BMap();
	void clear();
	void enter(Univ k, Univ v); // k must be new
	int contains(Univ k) const;
	Univ retrieve(Univ k, int& present) const;
	Univ retrieve(Univ k) const; // ret 0 if absent (die if v==0)
	Univ get(Univ k) const;	     // die if absent
	Univ remove(Univ k);	     // ret prevv or 0(ambiguous)
	Univ replace(Univ k, Univ v); // ret prevv or 0(ambiguous)
	int num() const;
	int empty() const;
	void OK() const;
	static const BMap EMPTY;
// not recommented:
	Univ* specialretrieveaddr(Univ k) const; // ret 0 if absent
	Univ* specialretrieveenteraddr(Univ k);
	int specialadd(Univ k, Univ v);	// ret was_there
	void specialslowenter(Univ k, Univ v);
	POOLALLOCATION(BMap);
#if !( defined(__DECCXX) || defined(SOLCC42))
  protected:
#endif
	struct Node {
		POOLALLOCATION(Node);
		Univ k;
		Univ v;
		Node* n;
	};
	// {b,isize,inum} could be made Array<Node*> but for efficiency
	Node** b;		// is zero if empty set!
	int hashk(Univ k) const;
	void removeaux(int buckn, Node* n, Node* last);
  private:
    friend class BMapIter;
	int isize;		// 0 if !b
	int inum;
	int fbuckn;		// <= first non-zero index, 0 if !b ??
	void quickenter(Univ k, Univ v);
	Node* find(Univ k) const;
	void resize(int newsize);
	DISABLECOPY(BMap);
};

class MeshRandom;

/** \class BMapIter Mesh_Map.H utils/Mesh_Map.H
   \ingroup MeshesHash

  \brief Iterator for a hash table

    This class is templated in MapIter.
    */
class BMapIter {
  public:
	BMapIter(const BMap& hh);
	BMapIter(const BMap& hh, MeshRandom& r);
	~BMapIter();
	operator void*() const;
	void next();
	Univ key() const;
	Univ value() const;
  private:
	const BMap& h;
	int bn;
	const BMap::Node* n;
	void findrealfbuckn();
	void advance();
	void advancewrap();
	// shallow copy is safe
};

//----------------------------------------------------------------------------

inline int BMap::hashk(Univ k) const { return Conv<unsigned>::d(k)%isize; }

// no inum++, no resize, no fbuckn
inline void BMap::quickenter(Univ k, Univ v)
{
	register int buckn=hashk(k);
	register Node* n=new Node;
	n->k=k; n->v=v; n->n=b[buckn]; b[buckn]=n;
}

inline void BMap::enter(Univ k, Univ v)
{
	inum++;
	if (!b) resize(5);
	else if (inum>isize*3) resize((isize-1)*5+3);
	quickenter(k,v);
	fbuckn=0;
}

inline BMap::Node* BMap::find(Univ k) const
{
	if (!b) return 0;
	for (Node* n=b[hashk(k)];n;n=n->n)
		if (n->k==k) return n;
	return 0;
}

inline int BMap::contains(Univ k) const
{
	return find(k)?1:0;
}

inline Univ BMap::retrieve(Univ k, int& present) const
{
	Node* n=find(k); present=n?1:0; return n?n->v:0;
}

inline Univ BMap::retrieve(Univ k) const
{
	Node* n=find(k);
	// return n?assertv(n->v):0;
	return n?n->v:0;
}

inline Univ BMap::get(Univ k) const
{
	// return assertv(find(k))->v;
	return find(k)->v;
}

inline Univ BMap::replace(Univ k, Univ v)
{
	Node* n=find(k);
	Univ ov;
	return n?(ov=n->v,n->v=v,ov):0;
}

inline int BMap::num() const { return inum; }
inline int BMap::empty() const { return !inum; }

inline Univ* BMap::specialretrieveaddr(Univ k) const
{
	Node* n=find(k); return n?&n->v:0;
}

inline Univ* BMap::specialretrieveenteraddr(Univ k)
{
	if (!b) { specialslowenter(k,0); return &find(k)->v; }
	int buckn=hashk(k);
     Node* n = NULL;

	for (n=b[buckn];n;n=n->n)
		if (n->k==k) return &n->v;
	n=new Node; n->k=k; n->v=0; n->n=b[buckn]; b[buckn]=n;
	fbuckn=0;
	if (++inum>isize*3) { resize((isize-1)*5+3); return &find(k)->v; }
	return &n->v;
}

inline int BMap::specialadd(Univ k, Univ v)
{
	if (!b) { specialslowenter(k,v); return 1; }
	int buckn=hashk(k);
     Node* n = NULL;

	for (n=b[buckn];n;n=n->n)
		if (n->k==k) return 0;
	n=new Node; n->k=k; n->v=v; n->n=b[buckn]; b[buckn]=n;
	fbuckn=0;
	if (++inum>isize*3) resize((isize-1)*5+3);
	return 1;
}


// Hand optimize this one
inline void BMapIter::advance()
{
	// assertx(!n);  too expensive
	// for (++bn;bn<h.isize;bn++) if (n=h.b[bn]) break;
	register int lbn=bn;
	register const int lisize=h.isize;
	// cast added for DECCXX
	register const BMap::Node** const lb=(const BMap::Node**)h.b;
	register const BMap::Node* ln=0;
	for (++lbn;lbn<lisize;lbn++) {
		ln = lb[lbn];
		if ( ln ) break;
	}
	bn=lbn; n=ln;
}

inline BMapIter::~BMapIter() { }
inline BMapIter::operator void*() const { return (void*)n; }
inline Univ BMapIter::key() const { return n->k; }
inline Univ BMapIter::value() const { return n->v; }

inline void BMapIter::next()
{
	n=n->n;
	if (!n) advance();
}

//----------------------------------------------------------------------------

template<class K, class V> class MapIter;

/** \class MMap Mesh_Map.H utils/Mesh_Map.H
   \ingroup MeshesHash

  \brief Mapping between elements. Essentially a hash table. 

    This class is a templated version of BMap.

  Warning: do not use with objects that have expensive copy constructors
    */
template<class K, class V>
class MMap : public BMap {
  public:
	MMap() { }
	virtual ~MMap() { }
	inline void enter(K k, V v)
	{ BMap::enter(Conv<K>::e(k),Conv<V>::e(v)); }
	inline int contains(K k) const
	{ return BMap::contains(Conv<K>::e(k)); }
	inline V retrieve(K k, int& present) const
	{ return Conv<V>::d(BMap::retrieve(Conv<K>::e(k),present)); }
	inline V retrieve(K k) const
	{ return Conv<V>::d(BMap::retrieve(Conv<K>::e(k))); }
	inline V get(K k) const
	{ return Conv<V>::d(BMap::get(Conv<K>::e(k))); }
	inline V remove(K k)
	{ return Conv<V>::d(BMap::remove(Conv<K>::e(k))); }
	inline V replace(K k, V v)
	{ return Conv<V>::d(BMap::replace(Conv<K>::e(k),Conv<V>::e(v))); }
// 	typedef MapIter<K,V> Iter;
};

/** \class MapIter Mesh_Map.H utils/Mesh_Map.H
   \ingroup MeshesHash

  \brief Iterator for a hash table

    This class is a templated version of BMapIter.

  Warning: do not use with objects that have expensive copy constructors
    */
template<class K, class V>
class MapIter : public BMapIter {
  public:
	inline MapIter(const MMap<K,V>& hh) : BMapIter(hh) { }
	inline MapIter(const MMap<K,V>& hh, MeshRandom& r) : BMapIter(hh,r) { }
	inline ~MapIter() { }
	inline K key() const { return Conv<K>::d(BMapIter::key()); }
	inline V value() const { return Conv<V>::d(BMapIter::value()); }
};

#define ForMapKeyValue(S,T1,V1,T2,V2) \
{ for (MapIter<T1,T2> zz1(S);zz1;zz1.next()) \
  { T1 V1=zz1.key(); T2 V2=zz1.value();
#define ForMapKey(S,T1,T2,V) { for (MapIter<T1,T2> zz2(S);zz2;zz2.next()) \
			       { T1 V=zz2.key();
#define ForMapValue(S,T1,T2,V) { for (MapIter<T1,T2> zz3(S);(zz3);zz3.next()) \
				 { T2 V=zz3.value();
				   
#endif


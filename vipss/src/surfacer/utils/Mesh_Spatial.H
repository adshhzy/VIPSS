// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Spatial_h
#define Spatial_h

#include <utils/Rn_Defs.H>
#include <utils/Mesh_Map.H>
#include <utils/Mesh_Stack.H>
#include <utils/Mesh_Set.H>
#include <utils/Mesh_Pqueue.H>


class Bbox;

class Spatial {			// abstract class
  public:
	Spatial(int pgn);
	virtual ~Spatial();
  protected:
    friend class SpatialSearch;
	int gn;			// grid size
	double gni;		// 1/gn
	int indexinbounds(int i) const;
	int indicesinbounds(int ci[3]) const;
	int double2index(double fd) const;
	double index2double(int i) const;
	void point2indices(const R3Pt& p, int ci[3]) const;
	void indices2point(int ci[3], R3Pt& p) const;
	void indices2bbox(int ci[3], Bbox& bb) const;
	int encode(int ci[3]) const;
	void decode(int en, int ci[3]) const;
	// for SpatialSearch:
	// Add elements from cell ci[3] to priority queue with priority equal
	// to distance from pcenter squared.  May use set to avoid duplication.
	virtual void addcell(int ci[3], Pqueue<Univ>& pq, const R3Pt& pcenter,
			     Set<Univ>& set) const=0;
	// Refine distance estimate of first entry in pq (optional)
	virtual void pqrefine(Pqueue<Univ>& pq, const R3Pt& pcenter) const;
	virtual Univ pqid(Univ pqe) const=0; // given pq entry, return id
  private:
	DISABLECOPY(Spatial);
};

class PointSpatial : public Spatial {
  public:
	PointSpatial(int pgn);
	~PointSpatial();
	void clear();
	// id!=0
	void enter(Univ id, const R3Pt* pp); // note: pp not copied!
	void remove(Univ id, const R3Pt* pp);	// must exist, else die
  protected:
	void addcell(int ci[3], Pqueue<Univ>& pq, const R3Pt& pcenter,
		     Set<Univ>& set) const;
	Univ pqid(Univ pqe) const;
  private:
	struct Node {
		Node(Univ pid, const R3Pt* pp) : id(pid), p(pp) { }
		Univ id;
		const R3Pt* p;
	};
	MMap<int,Stack<Node*>*> map; // encoded cube index -> Stack
};

class ObjectSpatial : public Spatial {
  public:
	typedef double (*DISTF)(const R3Pt& p, Univ id);
	ObjectSpatial(int pgn, DISTF papproxf2, DISTF pexactf2);
	~ObjectSpatial();
	void clear();
	// id!=0
	// Enter an object that comes with a containment function: the
	// function returns true if the object lies within a given bounding
	// box.  A starting point is also given.
	void enter(Univ id, const R3Pt& startp,
		       WINbool (*fcontains)(const Bbox& bb));
	// Find the objects that could possibly intersect the segment (p1,p2).
	// The objects are not returned in the exact order of intersection!
	// However, once should_stop is set (ftest's return), the procedure
	// will keep calling ftest with all objects that could be closer.
	void searchsegment(const R3Pt& p1, const R3Pt& p2,
			   int (*ftest)(Univ id)) const;
  protected:
	void addcell(int ci[3], Pqueue<Univ>& pq, const R3Pt& pcenter,
		     Set<Univ>& set) const;
	void pqrefine(Pqueue<Univ>& pq, const R3Pt& pcenter) const;
	Univ pqid(Univ pqe) const;
  private:
	MMap<int,Stack<Univ>*> map; // encoded cube index -> Stack
	DISTF approxf2;
	DISTF exactf2;
};

// Search from a point.
class SpatialSearch {
  public:
	// pmaxdis is only a request, you may get objects that lie farther
	SpatialSearch(const Spatial& psp, const R3Pt& pp,
		      double pmaxdis=10.);
	~SpatialSearch();
	int done();
	Univ next(double* dis2=0); // ret id
  private:
    friend class Spatial;
	const Spatial& sp;
	const R3Pt pcenter;
	double maxdis;
	Pqueue<Univ> pq;	// pq of entries by distance
	int ssi[2][3];		// search space indices (extents)
	double disbv2;		// distance to search space boundary
	int axis;		// axis to expand next
	int dir;		// direction in which to expand next (0,1)
	Set<Univ> setevis;	// may be used by addcell()
	int ncellsv;
	int nelemsv;
	void getclosestnextcell();
	void expandsearchspace();
	void consider(int ci[3]);
	DISABLECOPY(SpatialSearch);
};

#endif


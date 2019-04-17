// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Bbox_h
#define Bbox_h

#include <utils/Rn_Defs.H>

/** \class Bbox Mesh_BBox.H utils/Mesh_BBox.H
 \ingroup MeshesHelper
 \brief Axis-aligned bounding boxes, defined by min/max points.

  Can do boolean operations on these bounding boxes.
  */
class Bbox {
  public:
      ///
	Bbox() { }
    /// Initialize with given end points
	Bbox(const R3Pt& pmin, const R3Pt& pmax);
    ///
	~Bbox() { }
    /// Make zero size
	void clear();
    /// Make infinitely big
	void infinite();
    /// return min (0) or max (1) bounding box points.
	R3Pt& operator[](int i) { return p[i]; }
    /// return min (0) or max (1) bounding box points.
	const R3Pt& operator[](int i) const { return p[i]; }
    /// Make this box contain the input box
	void takeunion(const Bbox& bb);
    /// Expand box (if necessary) so it contains this point
	void takeunion(const R3Pt& pp);
    /// Only keep the bounding box contained inside of both this one and the input one
	void intersect(const Bbox& bb);
    /// Does this bounding box contain the input one? 
	int inside(const Bbox& bb) const;
	/// uniform scaling into unit cube, centered on x & y, rest at z=0
	R4Matrix getFrameToCube() const;
	/// non-uniform scaling (?) into unit cube, centered on x & y, rest at z=0
	R4Matrix getFrameToSmallCube() const;
  private:
      /// p[0] is min point, p[1] is max point
	R3Pt p[2];
	// shallow copy ok
};

#endif


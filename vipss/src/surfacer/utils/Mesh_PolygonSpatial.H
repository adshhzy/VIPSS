// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef PolygonSpatial_h
#define PolygonSpatial_h

#include <utils/Mesh_Spatial.H>
#include <utils/Rn_Polygon.H>

class PolygonSpatial : public ObjectSpatial {
  public:
	PolygonSpatial(int pgn);
	~PolygonSpatial();
	// clear() inherited from ObjectSpatial, does not delete Polygons!
	// deleteclear() must be implemented by user
	void enter(const R3Polygon &poly); // poly is not copied!
	int firstAlongSegment(const R3Pt& p1, const R3Pt& p2,
			      const R3Polygon*& poly, R3Pt& pint) const;
};

#endif


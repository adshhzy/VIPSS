#ifndef polygon_mapping_DEFS
#define polygon_mapping_DEFS

#include <utils/Rn_Polygon.H>

/* -----------------------  Constants  ------------------------------- */

/* -----------------------  Classes  ------------------------------- */

/** \class UTILSPolygon Utils_Polygon_mapping.H utils/Utils_Polygon_mapping.H
  \ingroup UtilitiesGeom

  \brief Map a quadrilateral to the unit square and back again.

    This map can be combined to take any non-coplanar quadrilateral to another non-coplanar quadrilateral, mapping the edges of one to the egdes of the other. Preserves lines. <br>

The mapping is a 3X3 matrix, which takes (x,y,1) to (s,t,w). To go from point to point, divide s,t by w. 

  For Barycentric coordinates, see Utils_Barys.H

<br>   Files: 
- include/utils/Utils_Polygon_mapping.H
- utils/utils/UTILSPolygon.cpp
*/
class UTILSPolygon  {
private:
    UTILSPolygon &operator=(UTILSPolygon &);
    static R3Matrix &MakeMatrix(const R2Polygon &p);
    
public:
    /** Find the matrix that takes a quadrilateral to the unit square.
        The unit square is (0,0) to (1,1). The polygon must be 4-sided */
    static R3Matrix Polygon_to_unit( const R2Polygon &p );
    /** Find the matrix that takes the unit square to a quadrilateral.
        The unit square is (0,0) to (1,1). The polygon must be 4-sided */
    static R3Matrix Unit_to_polygon( const R2Polygon &p );
    /** Find the matrix that takes the quadrilateral to a quadrilateral.
        Combines polygon->unit->polygon. */
    static R3Matrix Polygon_to_polygon( const R2Polygon & in_polyFrom, 
                                       const R2Polygon & in_polyTo );

    /** Map the 2D vector through the 3X3 matrix.
        Turns the vector into an homogenous point. */
    static R2Vec Map(const R3Matrix &m, const R2Vec &v);
    /** Map the 2D point through the 3X3 matrix.
        Turns the point into an homogenous point. */
    static R2Pt  Map(const R3Matrix &m, const R2Pt &v);

    ///
    UTILSPolygon() {}
    ///
    virtual ~UTILSPolygon() {}

    ///
    static WINbool Test();
};

#endif


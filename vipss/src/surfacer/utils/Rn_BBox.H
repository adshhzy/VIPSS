#ifndef _RN_DEFS_BBOX_H
#define _RN_DEFS_BBOX_H

#include <utils/Rn_Defs.H>
#include <utils/Mesh_Array.H>


/** \class RNBBoxTC RN_BBox.H utils/RN_BBox.H
  \ingroup LinearGeometry
  \brief  An ND bouinding box.

 Defined by center, extent. <br>
 Files: <br>
   - include/utils/Rn_BBox.H <br>
   - include/utils/Rn_BBox_i.H <br>
   -  src/utils/geometry/Rn_BBox.cpp <br>
  */
template <class Point, class Vector, class Matrix>
class RNBBoxTC {
protected:
    /// Center of bounding box
    Point           m_ptCenter;
    /// Size/width of box in each direction
    SArray<double>  m_adExtent;
    /// Matrix to take unit cube, centered at origin, to bbox
    Matrix m_mat, m_matInv;

    /// Location of corners in space (derived)
    SArray<Point>   m_aptCorners; 

public:
    /**@name Access to data */
    //@{
    /** Dimension.(Depends on template type) */
    inline int     Dim() const ;
    /// 1D: two corners, 2D: four corners, etc.
    inline int     NCorners() const { return m_aptCorners.num(); }
    /// Get corner
    inline const Point & Corner ( const int in_i ) const ;
    /// Get Center of bounding box
    inline const Point & Center ( ) const ;
    /// Get width/height etc. (x,y,z...)
    inline double Extent ( const int in_i ) const ;
    /// Matrix that transforms cube at origin to bounding box
    inline Matrix Transform() const { return m_mat; }
    /// Matrix that transforms bbox to cube at origin
    inline Matrix InvTransform() const { return m_matInv; }
    /// Extents in all dimensions
    inline Point Scale() const;
    //@}

    /**@name Changing data */
    //@{
    /// Move the bounding box
    inline void ApplyTransform( const Matrix & );
    /// width/height of box
    inline void SetExtents( const Point & in_ptExtents );
    /// Reset the matrix to the identity
    inline void ResetMatrix( );
    //@}

    /**@name Inside/outside/closest routines */
    //@{
    /// Is it on or inside the boundary?
    inline WINbool Inside( const Point & ) const;
    /// Is it within epsilon of the boundary?
    inline WINbool IsBoundary( const Point & ) const;
    /// Project onto cube
    inline Point   Project( const Point & ) const;
    /// Intersect ray with cube
    inline WINbool IntersectRay( const Point &, const Vector &, Array<Point> &out_apt ) const;
    //@}

    /**@name Comparisons with other bboxes */
    //@{
    /// Does an approximate (RNEpsilon_d) comparison
    inline WINbool   operator==( const RNBBoxTC<Point, Vector, Matrix> & ) const;
    //@}
    
    
    /**@name Construction and assignment routines */
    //@{
    ///
    inline RNBBoxTC &operator=(const RNBBoxTC<Point, Vector, Matrix> &p);
    /// 
    inline RNBBoxTC( const RNBBoxTC<Point, Vector, Matrix> &in_bbox ) ;
    /// Default constructor is unit square/cube
    inline RNBBoxTC( );
    /// 
    inline RNBBoxTC( const Point &in_ptCenter, const Point &in_ptScl );
    //@}

    ///
    virtual inline ~RNBBoxTC() {}
    
    /**@name Read, write, print */
    //@{
    /// Debugging print.
    inline void    Print() const ;
    /// Writes points to file.
    inline void    Write ( ofstream &out ) const ;
    /// Reads points from file.
    inline WINbool Read  ( ifstream &in );
    //@}

    /// Check closest point and intersection routines
    static WINbool Test();
};


#include "Rn_BBox_i.H"


/** \addtogroup LinearGeometry */
//@{
/// 1D bbox (line)
typedef RNBBoxTC<R1Pt, R1Vec, R2Matrix> R1BBox;
/// 2D bbox (square)
typedef RNBBoxTC<R2Pt, R2Vec, R3Matrix> R2BBox;
/// 3D bbox (cube)
typedef RNBBoxTC<R3Pt, R3Vec, R4Matrix> R3BBox;
//@}
#endif


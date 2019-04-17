#if _MSC_VER > 1000
#pragma once
#endif
// _MSC_VER > 1000

#ifndef __BASETYPE_PROJECTIVE_H
#define __BASETYPE_PROJECTIVE_H


/** \defgroup ProjectiveGeometry Homogenous 2D and 3D Points and Pluecker lines

\ingroup SimpleGeometry  

Note: operators such as addition and Length are written as global methods, so they don't appear in the class definition but as functions in the documentation

  I don't know if these classes are actually used anywhere...

  Copyright (C) 1998 Microsoft Corporation.  All Rights Reserved. <br>
<br>
  File:       include/utils/Rn_Projective.H, include/utils/Rn_Projective_i.H<br>
<br>
  Content:     2D and 3D projective - space classes <br>
  -          R2HPointTC   - 2D Projective point <br>
  -          R3HPointTC   - 3D Projective point <br>
  -          R3HPlaneTC   - 3D Projective plane <br>
  -          R3HLineTC    - 3D Pluecker line <br>
 <hr>
*/
//@{
//@}


/** 
\class R2HPointTC Rn_Projective.H utils/Rn_Projective.H
\ingroup ProjectiveGeometry
\brief  An R2 homogenous point */
template<class Coord>
class R2HPointTC
{
public:
    ///
    R2HPointTC() { }
    ///
    R2HPointTC( const Coord& _u, const Coord& _v, const Coord& _w ):
        v(_u, _v, _w) { }
    ///
    inline R2HPointTC( const R2PointTC<Coord>& p );
    ///
    inline operator R2PointTC<Coord> ();

public:
    ///
    Coord u, v, w;
};


/** 
\class R3HPointTC Rn_Projective.H utils/Rn_Projective.H
\ingroup ProjectiveGeometry
\brief  An R3 homogenous point */
template<class Coord>
class R3HPointTC
{
public:
    ///
    R3HPointTC() { }
    ///
    inline R3HPointTC( const Coord& _x, const Coord& _y, const Coord& _z, const Coord& _w );
    ///
    inline R3HPointTC( const R3PointTC<Coord>& p );
    ///
    inline operator R3PointTC<Coord> ();

public:
    ///
    Coord x, y, z, w;
};


/** 
\class R3HPlaneTC Rn_Projective.H utils/Rn_Projective.H
\ingroup ProjectiveGeometry
\brief  An R3 homogenous plane */
template<class Coord>
class R3HPlaneTC
{
public:
    ///
    R3HPlaneTC() { }
    ///
    inline R3HPlaneTC( const Coord& _x, const Coord& _y, const Coord& _z, const Coord& _w );
    ///
    inline R3HPlaneTC( const R3PointTC<Coord>& p1,
                       const R3PointTC<Coord>& p2,
                       const R3PointTC<Coord>& p3);


public:
    ///
    Coord x, y, z, w;
};


/** 
\class R3HLineTC Rn_Projective.H utils/Rn_Projective.H
\ingroup ProjectiveGeometry
\brief  R3HLineTC: Pluecker line class */
template<class Coord>
class R3HLineTC
{
public:
    ///
    R3HLineTC() { }
    ///
    inline R3HLineTC( const R3HPointTC<Coord>&, const R3HPointTC<Coord>& );
    ///
    inline R3HLineTC( const R3HPlaneTC<Coord>&, const R3HPlaneTC<Coord>& );


public:
    ///
    Coord pl[6];
};



#include "Rn_Projective_i.H"

#endif
// #ifndef __BASETYPE_PROJECTIVE_H


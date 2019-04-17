#ifndef _UTILS_GEOM_ARRAY_H
#define _UTILS_GEOM_ARRAY_H

/**@name Geometric utilities
  \ingroup UtilitiesGeom
  \brief Min, max, and average, standard deviation of arrays of points/vectors. 
   
   A set of templated methods for dealing with arrays of points or vectors.
   For example, average the list, find the bounding box, etc. */
//@{
/// Find the min of each of the elemnets of the array. Assumes WINmin can sort the elements of geom.
template<class geom> inline 
geom GEOMMin( const Array<geom> &in_a )
{
    geom g;

    for (int iD = 0; iD < g.Dim(); iD++ )
        g[iD] = 1e30;

    for ( int i = 0; i < in_a.num(); i++ ) {
        for (int iD = 0; iD < g.Dim(); iD++ )
            g[iD] = WINmin( g[iD], in_a[i][iD] );
    }

    return g;
}

/// Find the max of each of the elemnets of the array. Assumes WINmin can sort the elements of geom.
template<class geom> inline 
geom GEOMMax( const Array<geom> &in_a )
{
    geom g;

    for (int iD = 0; iD < g.Dim(); iD++ )
        g[iD] = -1e30;

    for ( int i = 0; i < in_a.num(); i++ ) {
        for (int iD = 0; iD < g.Dim(); iD++ )
            g[iD] = WINmax( g[iD], in_a[i][iD] );
    }

    return g;
}

/// Find the average of the elemnets of the array.
template<class geom> inline 
geom GEOMAvg( const Array<geom> &in_a )
{
    geom g;

    for (int iD = 0; iD < g.Dim(); iD++ )
        g[iD] = 0.0;

    for ( int i = 0; i < in_a.num(); i++ ) {
        for (int iD = 0; iD < g.Dim(); iD++ )
            g[iD] += in_a[i][iD];
    }

    if ( in_a.num() == 0 )
        return g;

    const double dDiv = 1.0 / (double) in_a.num();
    for (int iD = 0; iD < g.Dim(); iD++ )
        g[iD] *= dDiv;

    return g;
}


/// Find the average and the standard deviation of the array elements
template<class geom> inline 
geom GEOMAvgSD( const Array<geom> &in_a, geom &out_geomSD )
{
    const geom geomAvg = GEOMAvg( in_a );

    for (int iD = 0; iD < out_geomSD.Dim(); iD++ )
        out_geomSD[iD] = 0.0;

    for ( FORINT i = 0; i < in_a.num(); i++ ) {
        for (int iD = 0; iD < out_geomSD.Dim(); iD++ )
            out_geomSD[iD] += ( in_a[i][iD] - geomAvg[iD] ) * ( in_a[i][iD] - geomAvg[iD] );
    }

    if ( in_a.num() < 2 )
        return geomAvg;

    const double dDiv = 1.0 / (double) ( in_a.num() - 1.0 );
    for (int iD = 0; iD < out_geomSD.Dim(); iD++ )
        out_geomSD[iD] *= dDiv;

    return geomAvg;
}


/// Find the average, min, max and the standard deviation of the array elements
template<class geom> inline 
void GEOMAvgSDMinMax( const Array<geom> &in_a, 
                      geom & out_geomAvg, 
                      geom & out_geomSD,
                      geom & out_geomMin,
                      geom & out_geomMax )
{
    const geom geomAvg = GEOMAvg( in_a );

    for (int iD = 0; iD < out_geomMin.Dim(); iD++ ) {
        out_geomMin[iD] = 1e30;
        out_geomMax[iD] = -1e30;
        out_geomSD[iD] = 0.0;
    }

    for ( FORINT i = 0; i < in_a.num(); i++ ) {
        for (int iD = 0; iD < out_geomSD.Dim(); iD++ ) {
            out_geomSD[iD] += ( in_a[i][iD] - geomAvg[iD] ) * ( in_a[i][iD] - geomAvg[iD] );
            out_geomMin[iD] = WINmin( out_geomMin[iD], in_a[i][iD] );
            out_geomMax[iD] = WINmax( out_geomMax[iD], in_a[i][iD] );
        }
    }

    if ( in_a.num() < 2 )
        return geomAvg;

    const double dDiv = 1.0 / (double) ( in_a.num() - 1.0 );
    for (int iD = 0; iD < out_geomSD.Dim(); iD++ )
        out_geomSD[iD] *= dDiv;

    return geomAvg;
}

//@}
#endif


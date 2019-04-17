#ifndef BARY_POLY_H
#define BARY_POLY_H

#include <utils/Pm_MeshLite.H>

WINbool CalculateWeights( const PMeshLite &in_mesh, const R3Pt &in_pt, Array<double> &out_adWeights );

/// Sorts points by angle in tangent plane first
double LaplacianWeights( const R3Pt &in_pt, const R3Vec &in_vecNorm, const Array<R3Pt> &in_apt, Array<double> &out_adWeights );

/// Assumes points ordered
double LaplacianWeights( const R3Pt &in_pt, const Array<R3Pt> &in_apt, Array<double> &out_adWeights );

// Project onto plane by distributing angle error and keeping length
void ProjectPlane( const R3Pt &in_pt, const Array<R3Pt> &in_apt, Array<R2Pt> &out_apt2D );

/// Uses projection to sort
double FloaterWeights( const R3Pt &in_pt, const R3Vec &in_vecNorm, const Array<R3Pt> &in_apt, Array<double> &out_adWeights );

/// Assumes points are already sorted
double FloaterWeights( const R3Pt &in_pt, const Array<R3Pt> &in_apt, Array<double> &out_adWeights );

namespace Strain {

/// Find the symmetric matrix that aligns the two points
R3Matrix MatrixAlign(  const Array<R3Pt> &in_apt1, const Array<R3Pt> &in_apt2 );
    
void DiscreteStrain( const Array<R3Pt> &in_apt1, const Array<R3Pt> &in_apt2, const R3Vec &in_vecNormPts2,
                                    R3Matrix &out_matF, R3Matrix &out_matR1To2 );
    
/// Calculate the max,min strain of the deformation between two sets of points
pair<double, double> DiscreteStrain( const Array<R3Pt> &in_apt1, const Array<R3Pt> &in_apt2, const R3Vec &in_vecNormPts2,
                                    R3Matrix &out_matF  );

/// Same as above, but calculates the normal for you.
pair<double, double> DiscreteStrain( const Array<R3Pt> &in_apt1, const Array<R3Pt> &in_apt2, R3Matrix &out_matF  );

/// Same as above, but uses vecDs and vecDt to determine normal
/// Returns the 2D tensor in out_matFij
pair<double, double> DiscreteStrain( const Array<R3Pt> &in_apt1, const Array<R3Pt> &in_apt2, const R3Vec &in_vecDs, const R3Vec &in_vecDt, R3Matrix &out_matFij );


/// Calculates strains based on fitting a surface to the points and using derivatives
/// Assumes points are sampled at same s,t values
/// Returns the 2D tensor in out_matFij
pair<double, double> AnalyticStrain( const Array<R3Pt> &in_apt1, const Array<R3Pt> &in_apt2 );

/// Works from the raw derivatives - assumes using the same parameter space, i.e., f(s,t) = f'(s,t) everywhere
/// Returns the 2D tensor in out_matFij
pair<double, double> AnalyticStrain( R2Matrix &out_matFij, const R3Vec &in_vecDs1, const R3Vec &in_vecDt1, const R3Vec &in_vecDs2, const R3Vec &in_vecDt2 );

void AnalyticStrain3D(const R3Vec &vSA, const R3Vec &vTA, const R3Vec &vSB, const R3Vec &vTB,
                                   double &out_E1, double &out_E2, double &out_e1, double &out_e2);
}
#endif

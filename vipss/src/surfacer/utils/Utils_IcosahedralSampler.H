#ifndef UTILS_ICOSAHEDRALSAMPLER_DEFS_H
#define UTILS_ICOSAHEDRALSAMPLER_DEFS_H

#include <utils/Rn_Defs.H>

/** \class UTILSIcosahedralSampler Utils_IcosahedralSampler.H utils/Utils_IcosahedralSampler.H
\ingroup UtilitiesGeom
\brief Samples evenly spaced on an icosahedran.

  You can subsample by spliting each triangle. 

  Files: <br>
  - include/utils/Utils_IcosahedralSampler.H
  - src/utils/utils/UTILSIcosahedralSampler.cpp */
class UTILSIcosahedralSampler {
private:
    SArray<R3Vec>   m_avecVertices;
    SArray<R3Pt_i>  m_aiTriangles; 	// indices of triangles vertices
    SArray<R2Pt_i>  m_aiEdges; 	// indices of edge vertices

    R3Vec RectangularCoordsOf(double in_dPhi, double in_dTheta );
    void TriangleEdges(const R3Pt_i &in_aiTri, SArray<R2Pt_i> &out_aiEdges);
    void GetEdgeSamples(int in_iN, Array<R3Vec> &);
    void GetFaceSamples(int in_iN, Array<R3Vec> &);
    
public:
    /** Get samples as unit vectors.
      @param in_iN number of times to divide edge of icosohedran.
      @param out_avecSamples The samples */
    void GetSamples(int in_iN, Array<R3Vec> &out_avecSamples );

    ///
    UTILSIcosahedralSampler() ;
    ///
    ~UTILSIcosahedralSampler() {}
    
};

#endif


#ifndef _PLY_PMESH_LITE_PLY_H__
#define _PLY_PMESH_LITE_PLY_H__

#include <utils/Pm_MeshLite.H>

extern void PLYReadPlyToPMeshLite( const char * in_strFName, PMeshLite &out_mesh, const WINbool in_bReverse );
extern void PLYWritePMeshLiteToPly( const PMeshLite &in_mesh, const char * in_strFName );

#endif

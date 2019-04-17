/***** implicit.h */

/* header file for implicit surface polygonizer, implicit.c */

#ifndef IMPLICIT_HDR
#define IMPLICIT_HDR

#include <utils/Rn_Defs.H>

typedef struct {                   /* surface vertex */
    R3Pt  position;
    R3Vec normal;        /* position and surface normal */
} VERTEX;

typedef struct {                   /* list of vertices in polygonization */
    int count, max;                /* # vertices, max # allowed */
    VERTEX *ptr;                   /* dynamically allocated */
} VERTICES;



#ifdef __cplusplus
extern "C" {
#endif

bool polygonize (
    double (*function)(const R3Pt &in_pt),
    double size,
	int bounds,
    const R3Pt &in_ptStart,
	int (*triproc)(int i1, int i2, int i3, VERTICES vertices),
	void (*vertproc)(VERTICES vertices)	
	);

/* see implicit.c for explanation of arguments */

#ifdef __cplusplus
}
#endif

#endif




#include "Polygonizer.h"

//<rts> #include <misc/Srf_Sweep.H>
//#include <fitting/Fit_RBFHermite.H>

/***** implicit.c */

/*
 * ANSI C code from the article
 * "An Implicit Surface Polygonizer"
 * by Jules Bloomenthal
 * in "Graphics Gems IV", Academic Press, 1994
 *

 * implicit.c
 *     an implicit surface polygonizer, translated from Mesa to ANSI C
 *     applications should call polygonize()
 *
 * to compile a test program for ASCII output:
 *     cc implicit.c -o implicit -lm
 *
 * to compile a test program for display on an SGI workstation:
 *     cc -DSGIGFX implicit.c -o implicit -lgl_s -lm
 
 * to compile for subroutine use without main:
 *     cc -DNOMAIN -c implicit.c
 *
 * Authored by Jules Bloomenthal, Xerox PARC.
 * Copyright (c) Xerox Corporation, 1991.  All rights reserved.
 * Permission is granted to reproduce, use and distribute this code for
 * any and all purposes, provided that this notice appears in all copies
 *
 * Last modified 11 Jan 95 by Jules Bloomenthal and Paul Heckbert
 */

#include <math.h>
#ifdef WIN32
#include <malloc.h>
#elif  defined(__linux__)
#include <stdlib.h>
#else
#include <sys/malloc.h>
#endif
#include <stdio.h>
#include <sys/types.h>

#define RES     10 /* # converge iterations    */

#define L       0  /* left direction:   -x, -i */
#define R       1  /* right direction:  +x, +i */
#define B       2  /* bottom direction: -y, -j */
#define T       3  /* top direction:    +y, +j */
#define N       4  /* near direction:   -z, -k */
#define F       5  /* far direction:    +z, +k */
#define LBN     0  /* left bottom near corner  */
#define LBF     1  /* left bottom far corner   */
#define LTN     2  /* left top near corner     */
#define LTF     3  /* left top far corner      */
#define RBN     4  /* right bottom near corner */
#define RBF     5  /* right bottom far corner  */
#define RTN     6  /* right top near corner    */
#define RTF     7  /* right top far corner     */

/* the LBN corner of cube (i, j, k), corresponds with location
 * (start.x+(i-.5)*size, start.y+(j-.5)*size, start.z+(k-.5)*size) */

#define RAND()    ((rand()&32767)/32767.)  /* random number, 0--1 */
#define HASHBIT   (5)
#define HSIZE     (size_t)(1<<(3*HASHBIT)) /* hash table size (32768) */
#define MASK      ((1<<HASHBIT)-1)
#define HASH(i,j,k) \
((((((i)&MASK)<<HASHBIT)|((j)&MASK))<<HASHBIT)|((k)&MASK))
#define BIT(i, bit) (((i)>>(bit))&1)
#define FLIP(i,bit) ((i)^1<<(bit)) /* flip the given bit of i */

double RNEpsilon_d = 1e-15;
float RNEpsilon_f = 1e-6;


typedef struct {                   /* test the function for a signed value */
    R3Pt p;                       /* location of test */
    double value;                   /* function value at p */
    int ok;                        /* if value is of correct sign */
} TEST;

typedef struct cornerlist {        /* list of corners */
    int i, j, k;                   /* corner id */
    double value;                   /* corner value */
    struct cornerlist *next;       /* remaining elements */
} CORNERLIST;

typedef struct {                   /* partitioning cell (cube) */
    int i, j, k;                   /* lattice location of cube */
    CORNERLIST *corners[8];        /* eight corners */
} CUBE;

typedef struct cubes {             /* linked list (stack) of cubes */
    CUBE cube;                     /* a single cube */
    struct cubes *next;            /* remaining elements */
} CUBES;

typedef struct centerlist {        /* list of cube locations */
    int i, j, k;                   /* cube location */
    struct centerlist *next;       /* remaining elements */
} CENTERLIST;

typedef struct edgelist {          /* list of edges */
    int i1, j1, k1, i2, j2, k2;    /* edge corner ids */
    int vid;                       /* vertex id */
    struct edgelist *next;         /* remaining elements */
} EDGELIST;

typedef struct intlist {           /* list of integers */
    int i;                         /* an integer */
    struct intlist *next;          /* remaining elements */
} INTLIST;

typedef struct intlists {          /* list of list of integers */
    INTLIST *list;                 /* a list of integers */
    struct intlists *next;         /* remaining elements */
} INTLISTS;

typedef struct process {           /* parameters, function, storage */
    double (*function)
      (const R3Pt &in_pt);         /* implicit surface function */
    int (*triproc)(int i1, int i2,
      int i3, VERTICES vertices);  /* triangle output function */
    double size, delta;             /* cube size, normal delta */
    int bounds;                    /* cube range within lattice */
    R3Pt start;                   /* start point on surface */
    CUBES *cubes;                  /* active cubes */
    VERTICES vertices;             /* surface vertices */
    CENTERLIST **centers;          /* cube center hash table */
    CORNERLIST **corners;          /* corner value hash table */
    EDGELIST **edges;              /* edge and vertex id hash table */
} PROCESS;


void freeprocess (PROCESS *p);

int dotet (CUBE *cube, int c1, int c2, int c3, int c4, PROCESS *p);

int setcenter(CENTERLIST *table[], int i, int j, int k);

int vertid (CORNERLIST *c1, CORNERLIST *c2, PROCESS *p);


char *mycalloc (int nitems, int nbytes);
CORNERLIST *setcorner (PROCESS *p, int i, int j, int k);

void converge ( const R3Pt &in_p1, const R3Pt &p2, double v,
                double (*function)(const R3Pt &in_pt),
                R3Pt &p);

TEST find (int sign, PROCESS *p, const R3Pt &in_pt);


#ifndef NOMAIN


/**** A Test Program ****/

/* sphere: an inverse square function (always positive) */

double sphere (const R3Pt &in_pt) {
    double rsq = in_pt[0] * in_pt[0] + in_pt[1] * in_pt[1] + in_pt[2] * in_pt[2];
    return 1.0/(rsq < 0.00001? 0.00001 : rsq);
}

/* blob: a three-pole blend function, try size = .1 */

double blob (const R3Pt &in_pt) {
    return 4.0 -sphere( in_pt + R3Vec(0.5, -0.5,-0.5) )
               -sphere( in_pt + R3Vec(-0.5, 0.5,-0.5) )
               -sphere( in_pt + R3Vec(-0.5, -0.5,0.5) );
}
//#define SGIGFX
#ifdef SGIGFX /**************************************************************/
#include <GL/gl.h>
#include "gl.h"
#include "gl/device.h"

static int shadepolys = 1, nvertices = 0, ntriangles = 0;


/* triangle: called by polygonize() for each triangle; set SGI lines */

int triangle (int i1, int i2, int i3, VERTICES vertices) {
    int i, ids[3];
    ids[0] = i3;
    ids[1] = i2;
    ids[2] = i1;
    if (shadepolys) bgnpolygon(); else bgnclosedline();
    for (i = 0; i < 3; i++) {
        VERTEX *vv = &vertices.ptr[ids[i]];
        if (shadepolys) n3f (&vv->normal.x);
        v3f (&vv->position.x);
        } 
    if (shadepolys) endpolygon(); else endclosedline();
    ntriangles++;
    nvertices = vertices.count;
    return 1;
}


/* main: call polygonize() with blob function, display on SGI */

void PolygonizeMainDisplay (int ac, char *av[]) {
    char *err;
    int sp = shadepolys = ac < 2 || strcmp(av[1], "-lines");
    double size = shadepolys? 0.2 : 0.1;

    keepaspect(1, 1);
    winopen("implicit");
    doublebuffer ();

    if (shadepolys) {
        Matrix ident = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
        static double material[] = {
            AMBIENT, .05, .05, .05,
            DIFFUSE, 1., 1., 1.,
            SPECULAR, .5, .5, .5,
            SHININESS, 10,
            LMNULL,
            };
        static double lightModel[] = {
            AMBIENT, .05, .05, .05,
            LOCALVIEWER, 1.,
            LMNULL,
            };
        static double light[] = {
            LCOLOR, 1., 1., 1.,
            POSITION, 0., 0., 1., 0.,
            LMNULL,
            };
        RGBmode ();
        gconfig ();
        mmode (MVIEWING);
        loadmatrix (ident);
        lsetdepth (getgdesc (GD_ZMIN), getgdesc (GD_ZMAX));
        zbuffer (TRUE);
        backface (TRUE);
        czclear(0x404040, getgdesc (GD_ZMAX));
        lmdef (DEFMATERIAL, 1, 0, material);
        lmdef (DEFLIGHT, 1, 0, light);
        lmdef (DEFLMODEL, 1, 0, lightModel);
        lmbind (MATERIAL, 1);
        lmbind (LMODEL, 1);
        lmbind (LIGHT0, 1);
    }
    else {
        gconfig();
        color(7);
        clear();
    }
    swapbuffers();
    perspective(450, 1.0/1.0, 0.1, 10.0);
    translate(0.0, 0.0, -3.0);

    makeobj(1);
    if (shadepolys) {
        RGBcolor (0, 255, 0);
        lmcolor (LMC_DIFFUSE);
    }
    if ((err = polygonize(blob, size, 20, 0.,0.,0., triangle))
         != NULL) {
             fprintf(stderr, "%s\n", err);
             exit(1);
    }
    closeobj();
    printf ("%d vertices, %d triangles\n", nvertices, ntriangles);

    pushmatrix();
    for (;;) { /* spin the object */
        reshapeviewport();
        if (shadepolys) {
            czclear(0x404040, getgdesc (GD_ZMAX));
            lmcolor(LMC_DIFFUSE);
        }
        else {color(7); clear(); color(100);}
        rot(0.8, 'x');
        rot(0.3, 'y');
        rot(0.1, 'z');
        callobj(1);
        swapbuffers();
    }
}

#else /***********************************************************************/

int gntris;          /* global needed by application */
VERTICES gvertices;  /* global needed by application */



Array< int > g_aiFaces0, g_aiFaces1, g_aiFaces2;
R3Pt m_ptMin;
R3Pt m_ptMax;


/* triangle: called by polygonize() for ea. triangle; write to stdout */

int triangle (int i1, int i2, int i3, VERTICES vertices) {
    gvertices = vertices;
    gntris++;
    fprintf(stdout, "%d %d %d\n", i1, i2, i3);



    g_aiFaces0 += i1;

    g_aiFaces1 += i2;

    g_aiFaces2 += i3;


    return 1;
}

int TriProc(int in_i1, int in_i2, int in_i3, VERTICES vs) {
	const R3Pt pt = vs.ptr[in_i1].position;
	
	bool bOutside = false;
	for ( int j = 0; j < 3; j++ ) {
		const double dWidth = m_ptMax[j] - m_ptMin[j];
		if ( pt[j] < m_ptMin[j] - dWidth * 0.1 ) {
			bOutside = true;
		}
		if ( pt[j] > m_ptMax[j] + dWidth * 0.1 ) {
			bOutside = true;
		}
	}
	if ( bOutside == false ) {
		//s_srfData->m_afaceSurface.push_back( R3Pt_i( in_i1, in_i2, in_i3 ) );
		g_aiFaces0 += in_i1;
		g_aiFaces1 += in_i2;
		g_aiFaces2 += in_i3;
	}
	
	return 1;
}


/* main: call polygonize() with blob function
 * write points-polygon formatted data to stdout */


FILE *g_out;
/*<rts>
void PolygonizeMainFile (FILE *out, double in_fStep, int in_iNCubes ) {

    char *err;

    

    gntris = 0;

    R3Pt ptOut;



    g_out = out;

    SRFSweep::FirstPt(ptOut);

    if ((err = polygonize(SRFSweep::DistFunc, in_fStep, in_iNCubes, ptOut, triangle))

         != NULL) {

             fprintf(out, "%s\n", err);

             exit(1);

    }



}
void PolygonizeMainFile (FILE *out, double in_fStep, int in_iNCubes,  R3Pt in_pt, R3Pt in_ptMin, R3Pt in_ptMax) {

	char *err;
	//gntris = 0;
	//R3Pt ptOut;
	g_out = out;
	m_ptMin = in_ptMin;
	m_ptMax = in_ptMax;

	//SRFSweep::FirstPt(ptOut);

	if ((err = polygonize(FITRbfHermite::DistFunc, in_fStep, in_iNCubes, in_pt, TriProc))
		!= NULL) {
			fprintf(out, "%s\n", err);
			exit(1);
	}

}*/

void Write ( ) 

{

    for (int i = 0; i < gvertices.count; i++) {

        VERTEX v;

        v = gvertices.ptr[i];

		//cout << "i+1 = " << i+1 << "\n";
		//cout << v.position << "\n";
		//cout << v.normal << "\n";

        fprintf(g_out, "Vertex %d %f %f %f normal={%f %f %f}\n",

            i+1,

            v.position[0], v.position[1], v.position[2],

            v.normal[0],   v.normal[1],   v.normal[2]);

    }

    for (FORINT i = 0; i < g_aiFaces0.num(); i++) {

		fprintf(g_out, "Face %d  %d %d %d\n", i+1, g_aiFaces0[i]+1, g_aiFaces2[i]+1, g_aiFaces1[i]+1 );

    }

}



#endif /*********** endif for SGIGFX *****************************************/
#endif /*********** endif for NOMAIN *****************************************/


/**** An Implicit Surface Polygonizer ****/



void testface ( int i, int j, int k,

                CUBE *old,

                int face, int c1, int c2, int c3, int c4,   PROCESS *p);


/* polygonize: polygonize the implicit surface function
 *   arguments are:
 *       double function (const R3Pt &in_pt)
 *           the implicit surface function given an arbitrary point
 *           return negative for inside, positive for outside
 *       double size
 *           width of the partitioning cube
 *       int bounds
 *           max. range of cubes (+/- on the three axes) from first cube
 *       double x, y, z
 *           coordinates of a starting point on or near the surface
 *           may be defaulted to 0., 0., 0.
 *       int triproc (i1, i2, i3, vertices)
 *               int i1, i2, i3 (indices into the vertex array)
 *               VERTICES vertices (the vertex array, indexed from 0)
 *           called for each triangle
 *           the triangle coordinates are (for i = i1, i2, i3):
 *               vertices.ptr[i].position.x, .y, and .z
 *           vertices are ccw when viewed from the out (positive) side
 *               in a left-handed coordinate system
 *           vertex normals point outwards
 *           return 1 to continue, 0 to abort
 *   returns error or NULL
 */

bool polygonize (
    double (*function)(const R3Pt &in_pt),
    double size,
    int bounds,
    const R3Pt &in_pt,
    int (*triproc)(int i1, int i2, int i3, VERTICES vertices),
	void (*vertproc)(VERTICES vertices))
    {
    int n;
    PROCESS p;
    TEST in, out;
    
    p.function = function;
    p.triproc = triproc;
    p.size = size;
    p.bounds = bounds;
    p.delta = size/(double)(RES*RES);

    /* allocate hash tables: */
    p.centers = (CENTERLIST **) mycalloc(HSIZE, sizeof(CENTERLIST *));
    p.corners = (CORNERLIST **) mycalloc(HSIZE,   sizeof(CORNERLIST *));
    p.edges   = (EDGELIST   **) mycalloc(2*HSIZE, sizeof(EDGELIST *));

    p.vertices.count = p.vertices.max = 0; /* no vertices yet */
    p.vertices.ptr = NULL;
    
    /* find point on surface, beginning search at (x, y, z):  */
    srand(1);
    in = find(1, &p, in_pt);
    out = find(0, &p, in_pt);
    if (!in.ok || !out.ok) {
        freeprocess(&p);
        if (!in.ok) printf ("in not ok\n");
        if (!out.ok) printf ("out not ok\n");
        cerr << "ERR: polyganizer can't find starting point\n";
		return false;
    }
    converge(in.p, out.p, in.value, p.function, p.start);

    /* push initial cube on stack: */
    p.cubes = (CUBES *) mycalloc(1, sizeof(CUBES)); /* list of 1 */
    p.cubes->cube.i = p.cubes->cube.j = p.cubes->cube.k = 0;
    p.cubes->next = NULL;

    /* set corners of initial cube: */
    for (n = 0; n < 8; n++)
        p.cubes->cube.corners[n] = \
            setcorner(&p, BIT(n,2), BIT(n,1), BIT(n,0));

    setcenter(p.centers, 0, 0, 0);

    while (p.cubes != NULL) { /* process active cubes till none left */
        CUBE c;
        CUBES *temp = p.cubes;
        c = p.cubes->cube;

        /* decompose into tetrahedra and polygonize: */
        if (!(dotet(&c, LBN, LTN, RBN, LBF, &p) &&
              dotet(&c, RTN, LTN, LBF, RBN, &p) &&
              dotet(&c, RTN, LTN, LTF, LBF, &p) &&
              dotet(&c, RTN, RBN, LBF, RBF, &p) &&
              dotet(&c, RTN, LBF, LTF, RBF, &p) &&
              dotet(&c, RTN, LTF, RTF, RBF, &p)))
         {
			 freeprocess(&p);
			 cerr << "ERR: polyganizeraborted";
			 return false;
         }

        /* pop current cube from stack */
        p.cubes = p.cubes->next;
        free((char *) temp);
        /* test six face directions, maybe add to stack: */
        testface(c.i-1, c.j, c.k, &c, L, LBN, LBF, LTN, LTF, &p);
        testface(c.i+1, c.j, c.k, &c, R, RBN, RBF, RTN, RTF, &p);
        testface(c.i, c.j-1, c.k, &c, B, LBN, LBF, RBN, RBF, &p);
        testface(c.i, c.j+1, c.k, &c, T, LTN, LTF, RTN, RTF, &p);
        testface(c.i, c.j, c.k-1, &c, N, LBN, LTN, RBN, RTN, &p);
        testface(c.i, c.j, c.k+1, &c, F, LBF, LTF, RBF, RTF, &p);
    }

    gvertices = p.vertices;
	vertproc( gvertices );

	//cout << "Starting to write\n";
    //Write();
    freeprocess(&p);

    return NULL;
}

/* freeprocess: free all allocated memory */

void freeprocess (PROCESS *p) {
    int index;
    CORNERLIST *l, *lnext;
    CENTERLIST *cl, *clnext;
    EDGELIST *edge, *edgenext;

    for (index = 0; index < HSIZE; index++)
        for (l = p->corners[index]; l; l = lnext) {
            lnext = l->next;
            free((char *) l);           /* free CORNERLIST */
        }
    for (index = 0; index < HSIZE; index++)
        for (cl = p->centers[index]; cl; cl = clnext) {
            clnext = cl->next;
            free((char *) cl);          /* free CENTERLIST */
        }
    for (index = 0; index < 2*HSIZE; index++)
        for (edge = p->edges[index]; edge; edge = edgenext) {
            edgenext = edge->next;
            free((char *) edge);        /* free EDGELIST */
        }
    free((char *) p->edges);            /* free array EDGELIST ptrs */
    free((char *) p->corners);          /* free array CORNERLIST ptrs */
    free((char *) p->centers);          /* free array CENTERLIST ptrs */
    if (p->vertices.ptr)
        free((char *) p->vertices.ptr); /* free VERTEX array */
}

/* testface: given cube at lattice (i, j, k), and four corners of face,
 * if surface crosses face, compute other four corners of adjacent cube
 * and add new cube to cube stack */

void testface (
    int i, int j, int k,
    CUBE *old,
    int face, int c1, int c2, int c3, int c4,   PROCESS *p)
    {
    CUBE cubeNew;
    CUBES *oldcubes = p->cubes;
    static int facebit[6] = {2, 2, 1, 1, 0, 0};
    int n, pos = old->corners[c1]->value > 0.0 ? 1 : 0;
    int bit = facebit[face];

    /* test if  no surface crossing, cube out of bounds, or prev. visited? */
    if ((old->corners[c2]->value > 0) == pos &&
        (old->corners[c3]->value > 0) == pos &&
        (old->corners[c4]->value > 0) == pos) return;
    if (abs(i) > p->bounds || abs(j) > p->bounds || abs(k) > p->bounds)
        return;
    if (setcenter(p->centers, i, j, k)) return;

    /* create new cube: */
    cubeNew.i = i;
    cubeNew.j = j;
    cubeNew.k = k;
    for (n = 0; n < 8; n++) cubeNew.corners[n] = NULL;
    cubeNew.corners[FLIP(c1, bit)] = old->corners[c1];
    cubeNew.corners[FLIP(c2, bit)] = old->corners[c2];
    cubeNew.corners[FLIP(c3, bit)] = old->corners[c3];
    cubeNew.corners[FLIP(c4, bit)] = old->corners[c4];
    for (n = 0; n < 8; n++)
        if (cubeNew.corners[n] == NULL) cubeNew.corners[n] =
            setcorner(p, i+BIT(n,2), j+BIT(n,1), k+BIT(n,0));

    /*add cube to top of stack: */
    p->cubes = (CUBES *) mycalloc(1, sizeof(CUBES));
    p->cubes->cube = cubeNew;
    p->cubes->next = oldcubes;
}


/* setpoint: set point location given lattice location */

void setpoint (R3Pt &out_pt, int i, int j, int k, PROCESS *p) 

{
    out_pt[0] = p->start[0]+((double)i-0.5) * p->size;
    out_pt[1] = p->start[1]+((double)j-0.5) * p->size;
    out_pt[2] = p->start[2]+((double)k-0.5) * p->size;
}


/* setcorner: return corner with the given lattice location
   set (and cache) its function value */

CORNERLIST *setcorner (PROCESS *p, int i, int j, int k) {
    /* for speed, do corner value caching here */
    int index = HASH(i, j, k);
    CORNERLIST *l = p->corners[index];
    R3Pt pt;
    
    for (; l != NULL; l = l->next)
        if (l->i == i && l->j == j && l->k == k) return l;
            
    setpoint (pt, i, j, k, p);
    l = (CORNERLIST *) mycalloc(1, sizeof(CORNERLIST));
    l->i = i; l->j = j; l->k = k;
    l->value = p->function(pt);
    l->next = p->corners[index];
    p->corners[index] = l;
    return l;
}


/* find: search for point with value of given sign (0: neg, 1: pos) */

TEST find (int sign, PROCESS *p, const R3Pt &in_pt)
    {
    int i;
    TEST test;
    double range = p->size;
    test.ok = 1;
    for (i = 0; i < 10000; i++) {

        const R3Vec vec(range*(RAND()-0.5), range*(RAND()-0.5), range*(RAND()-0.5));

        test.p = in_pt + vec;


        test.value = p->function(test.p);
        if (sign == (test.value > 0.0)) return test;
        range = range*1.0005; /* slowly expand search outwards */
    }
    test.ok = 0;
    return test;
}


/**** Tetrahedral Polygonization ****/


/* dotet: triangulate the tetrahedron
 * b, c, d should appear clockwise when viewed from a
 * return 0 if client aborts, 1 otherwise */

int dotet (CUBE *cube, int c1, int c2, int c3, int c4, PROCESS *p) {
    CORNERLIST *a = cube->corners[c1];
    CORNERLIST *b = cube->corners[c2];
    CORNERLIST *c = cube->corners[c3];
    CORNERLIST *d = cube->corners[c4];
    int index = 0, apos, bpos, cpos, dpos, e1, e2, e3, e4, e5, e6;
    if (apos = (a->value > 0.0)) index += 8;
    if (bpos = (b->value > 0.0)) index += 4;
    if (cpos = (c->value > 0.0)) index += 2;
    if (dpos = (d->value > 0.0)) index += 1;
    /* index now 4-bit number equal to one of the 16 possible cases */
    if (apos != bpos) e1 = vertid(a, b, p);
    if (apos != cpos) e2 = vertid(a, c, p);
    if (apos != dpos) e3 = vertid(a, d, p);
    if (bpos != cpos) e4 = vertid(b, c, p);
    if (bpos != dpos) e5 = vertid(b, d, p);
    if (cpos != dpos) e6 = vertid(c, d, p);
    /* 14 productive tet. cases (0000 and 1111 do not yield polygons */
    switch (index) {
        case 1:  return p->triproc(e5, e6, e3, p->vertices);
        case 2:  return p->triproc(e2, e6, e4, p->vertices);
        case 3:  return p->triproc(e3, e5, e4, p->vertices) &&
                        p->triproc(e3, e4, e2, p->vertices);
        case 4:  return p->triproc(e1, e4, e5, p->vertices);
        case 5:  return p->triproc(e3, e1, e4, p->vertices) &&
                        p->triproc(e3, e4, e6, p->vertices);
        case 6:  return p->triproc(e1, e2, e6, p->vertices) &&
                        p->triproc(e1, e6, e5, p->vertices);
        case 7:  return p->triproc(e1, e2, e3, p->vertices);
        case 8:  return p->triproc(e1, e3, e2, p->vertices);
        case 9:  return p->triproc(e1, e5, e6, p->vertices) &&
                        p->triproc(e1, e6, e2, p->vertices);
        case 10: return p->triproc(e1, e3, e6, p->vertices) &&
                        p->triproc(e1, e6, e4, p->vertices);
        case 11: return p->triproc(e1, e5, e4, p->vertices);
        case 12: return p->triproc(e3, e2, e4, p->vertices) &&
                        p->triproc(e3, e4, e5, p->vertices);
        case 13: return p->triproc(e6, e2, e4, p->vertices);
        case 14: return p->triproc(e5, e3, e6, p->vertices);
    }
    return 1;
}


/**** Storage ****/


/* mycalloc: return successful calloc or exit program */

char *mycalloc (int nitems, int nbytes) {
   char *ptr = (char *) calloc(nitems, nbytes);
   if (ptr != NULL) return ptr;
   fprintf(stderr, "can't calloc %d bytes\n", nitems*nbytes);
   exit(1);
}


/* setcenter: set (i,j,k) entry of table[]
 * return 1 if already set; otherwise, set and return 0 */

int setcenter(CENTERLIST *table[], int i, int j, int k) {
    int index = HASH(i, j, k);
    CENTERLIST *centerListNew, *l, *q = table[index];
    for (l = q; l != NULL; l = l->next)
        if (l->i == i && l->j == j && l->k == k) return 1;
    centerListNew = (CENTERLIST *) mycalloc(1, sizeof(CENTERLIST));
    centerListNew->i = i; centerListNew->j = j; centerListNew->k = k; centerListNew->next = q;
    table[index] = centerListNew;
    return 0;
}


/* setedge: set vertex id for edge */

void setedge (
    EDGELIST *table[],
    int i1, int j1, int k1, int i2, int j2, int k2, int vid)
    {
    unsigned int index;
    EDGELIST *edgeListNew;
    if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
        int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
    }
    index = HASH(i1, j1, k1) + HASH(i2, j2, k2);
    edgeListNew = (EDGELIST *) mycalloc(1, sizeof(EDGELIST));
    edgeListNew->i1 = i1; edgeListNew->j1 = j1; edgeListNew->k1 = k1;
    edgeListNew->i2 = i2; edgeListNew->j2 = j2; edgeListNew->k2 = k2;
    edgeListNew->vid = vid;
    edgeListNew->next = table[index];
    table[index] = edgeListNew;
}


/* getedge: return vertex id for edge; return -1 if not set */

int getedge (EDGELIST *table[],
             int i1, int j1, int k1, int i2, int j2, int k2)
    {
    EDGELIST *q;
    if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
        int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
    }
    q = table[HASH(i1, j1, k1)+HASH(i2, j2, k2)];
    for (; q != NULL; q = q->next)
        if (q->i1 == i1 && q->j1 == j1 && q->k1 == k1 &&
            q->i2 == i2 && q->j2 == j2 && q->k2 == k2)
            return q->vid;
    return -1;
}


/**** Vertices ****/
int addtovertices (VERTICES *vertices, VERTEX v);

void vnormal (const R3Pt &in_point, PROCESS *p, R3Vec &out_v);



/* vertid: return index for vertex on edge:
 * c1->value and c2->value are presumed of different sign
 * return saved index if any; else compute vertex and save */

int vertid (CORNERLIST *c1, CORNERLIST *c2, PROCESS *p) {
    VERTEX v;
    R3Pt a, b;
    int vid =
        getedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k);
    if (vid != -1) return vid;                /* previously computed */
    setpoint (a, c1->i, c1->j, c1->k, p);
    setpoint (b, c2->i, c2->j, c2->k, p);
    converge (a, b, c1->value, p->function, v.position); /* posn.  */
    vnormal(v.position, p, v.normal);                     /* normal */
    vid = addtovertices(&p->vertices, v);                   /* save   */
    setedge(p->edges, c1->i, c1->j, c1->k, c2->i, c2->j, c2->k, vid);
    return vid;
}


/* addtovertices: add v to sequence of vertices and return its id */

int addtovertices (VERTICES *vertices, VERTEX v) {
   if (vertices->count == vertices->max) {
      int i;
      VERTEX *vertexNew;
      vertices->max = vertices->count == 0 ? 10 : 2*vertices->count;
      vertexNew = (VERTEX *) mycalloc((unsigned)vertices->max,sizeof(VERTEX));
      for (i = 0; i < vertices->count; i++) vertexNew[i] = vertices->ptr[i];
      if (vertices->ptr != NULL) free((char *) vertices->ptr);
      vertices->ptr = vertexNew;
   }
   vertices->ptr[vertices->count++] = v;
   return (vertices->count-1);
}


/* vnormal: compute unit length surface normal at point */

void vnormal (const R3Pt &in_point, PROCESS *p, R3Vec &out_vec) {
    const double f = p->function(in_point);



    R3Vec vec(0,0,0);

    for ( int i = 0; i < 3; i++ ) {

        vec[i] = p->delta;

        out_vec[i] = p->function( in_point + vec ) - f;

        vec[i] = 0.0;

    }



    out_vec = UnitSafe( out_vec );

}


/* converge: from two points of differing sign, converge to surface */

void converge ( const R3Pt &in_p1, const R3Pt &in_p2, double v,
                double (*function)(const R3Pt &in_pt), 

                R3Pt &out_p)
{
    int i = 0;
    R3Pt pos, neg;
    if (v < 0) {

        pos = in_p2;

        neg = in_p1;

    }
    else {
        pos = in_p1;

        neg = in_p2;

    }
    for (;;) {

        out_p = Lerp( pos, neg, 0.5 );


        if (i++ == RES) return;
        if ((function((out_p))) > 0.0)
             {pos = out_p;}
        else {neg = out_p;}
    }
}

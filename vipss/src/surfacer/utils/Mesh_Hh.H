// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Hh_h
#define Hh_h

#include <WINSystemDefines.H>

#if defined(sgi) && !defined(__GNUG__)
#define SGICC
#endif

/** \defgroup MeshesAsserts Assertions, conversions 
  \ingroup MeshesHelper 

  Stand alone functions for debugging, conversion of pointers to void *
*/
//@{
///
extern int hhassert(int warncode, const char* s, const char* file, int line);
template<class T> inline T assertvaux(T expr, const char* s,
				      const char* file, int line)
{ if (!expr) hhassert(-1,s,file,line); return expr; }
/// Warning: ret: message_printed
#define MESHWarning(s) hhassert(1, s, __FILE__, __LINE__)
/// assertnever: ret: void
#define MESHassertnever(str) \
(void)hhassert(-1, "assertnever(" str ")", __FILE__, __LINE__)
//// assertx: ret: void
#define MESHassertx(expr) \
((void)((expr)?0:hhassert(-1, "assertx(" #expr ")", __FILE__, __LINE__)))
//// assertw: ret: expr_failed
#define MESHassertw(expr) \
((expr)?0:(hhassert(+0, "assertw(" #expr ")", __FILE__, __LINE__),1))
/// assertw1: ret: expr_failed
#define MESHassertw1(expr) \
((expr)?0:(hhassert(+1, "assertw1(" #expr ")", __FILE__, __LINE__),1))
/// assertv: ret: val
#define MESHassertv(e) \
assertvaux(e, (const char*)"assertv(" #e ")", (const char*)__FILE__, __LINE__)

///
extern void FlushStats();
///
extern void FlushWarnings();

class SUniv; 
typedef SUniv* Univ; // bogus class

/*
inline unsigned pointertounsigned(const void* p)    { return *((unsigned*)&p); }
inline const void* unsignedtopointer(unsigned i)    { return (const void*)i; }
inline double pointertodouble(const void* p)        { return *((double*)&p); }
inline const void* doubletopointer(const double& f) { return (const void*)(*((int*)(&f))); }
*/

template<class T>
class Conv {
  public:
	inline static Univ e(T e) { return (Univ)(const SUniv*)e; }
	inline static T d(Univ e) { return (T)(void*)e; }
};

template<>
class Conv<int> {
  public:
	inline static Univ e(int e) { const unsigned long int eLong = e; return Univ((const void*) eLong); }
	inline static int d(Univ e) { return *((unsigned*)&e); }
};

template<>
class Conv<unsigned int> {
  public:
	inline static Univ e(unsigned int e) { const unsigned long int eLong = e; return Univ((const void*) eLong); }
	inline static unsigned int d(Univ e) { return *((unsigned int*)&e); }
};

/// added for GNUG
/*
template<>
class Conv<unsigned> {
  public:
	inline static Univ e(unsigned e) { return Univ((const void*)(e)); }
	inline static unsigned d(Univ e) { return *((unsigned*)&e); }
};
*/

template<>
class Conv<double> {
  public:
	inline static Univ e(double e) { return Univ((const void*)(*((unsigned long int*)(&e)))); }
	inline static double d(Univ e) { return *((double*)&e); }
};

extern void* hmalloc(void* oldp, int size);
extern void hfree(void* e);

#ifndef PI
const double PI=3.14159265358979323846f;
#endif

template<class T>
inline void hqsort(T* ar, int n, int (*cmp)(const T*, const T*))
{ qsort(ar,n,sizeof(T),(int(*)(const void*,const void*))cmp); }

/// GNUG 2.5.0 sometimes requires :: before min() and max()
//template<class T> inline void swap(T* e1, T* e2) { T e=*e1; *e1=*e2; *e2=e; }
template<class T> inline T sign(T e) { return e>0?1:e<0?-1:0; }
template<class T> inline T square(T e) { return e*e; }

//inline double torad(double deg) { return deg/180*PI; }
//inline double todeg(double rad) { return rad/PI*180; }

inline double torad(double deg) { return deg/180*3.14159265358979323846f; }
inline double todeg(double rad) { return rad/3.14159265358979323846f*180; }

// Useful for "NEST { int localvar; ... }" so that indentation is correct
//inline double torad(double deg) { return deg/180*PI; }
//inline double todeg(double rad) { return rad/PI*180; }

/// Useful for "NEST { int localvar; ... }" so that indentation is correct
#define NEST
#define bcase break; case
#define ocase case
#define bdefault break; default
#define ForIndex(i,upb) { int zz=upb; for (int i=0;i<zz;i++) {
#define DummyEndFor }} // for gemacs formatting
#define For {{ for
#define EndFor }}
#define DISABLECOPY(class) class& operator=(class&); class(class&)
#define dummy_use(v) ((void)v)

// SHOWN(x*y);
// SHOW(point+vector);
// SHOWS("Could not open display");
// SHOWF("%s: Argument '%s' ambiguous, '%s' assumed\n",argv[0],arg,assumed);
// SHOWL;
// os << hform(" Endprogram: %dgons %dlines\n",ngons,nlines);

/// Show an expression
#define SHOW(x) (cerr << #x << " = " << (x))
/// Show an expression, and terminate with a newline
#define SHOWN(x) (cerr << #x << " = " << (x) << "\n")
/// Show a string
#define SHOWS(s) (cerr << (s) << "\n")
/// Formatted show, uses small buffer
extern void SHOWF(const char* format, ...);
/// Output a string prefixed with "# " on cerr, and on cout if _is_a_file
extern void SHOWDF(const char* format, ...);
/// Output a string prefixed with "# " on cout only if _is_a_file
extern void SHOWFF(const char* format, ...);
/// Show current line number and file
#define SHOWL (cerr << "Now in " << __FILE__ << " at line " << \
	       __LINE__ << "\n")
/// Format a string.  Uses a few rotating small buffers
/// use { char s[bigv]; sprintf() or ostrstream() } if you need a big buffer
extern const char* hform(const char* format, ...);

/// Allocate a duplicate of string, using new char[], use delete[]!
extern char* newString(const char* s);

extern int GetenvValue(const char* varname);

/// Returns fcntl return value
extern void SpawnOnSockets(int* fdi, int* fdo,
			   FILE** fi, FILE** fo, void (*spawn)());

/// Return absolute time, in seconds
extern double GetPreciseTime();
/// Return absolute time, in seconds
extern int GetTime();
/// Return string of time, without a trailing '\n'
extern const char* CTime();
/// Delay for some number of seconds
extern void Delay(double sec);

/// These could be made inline but would require including some .h files
extern void doubleToStd(void* in, void* ou);
extern void StdTodouble(void* in, void* ou);
extern void IntToStd(void* in, void* ou);
extern void StdToInt(void* in, void* ou);
extern void ShortToStd(void* in, void* ou);
extern void StdToShort(void* in, void* ou);

/// Prevent NaN's from appearing.
inline double myacos(double a)
{
	MESHassertw( ( fabs(a)<1.001 ) );
	if (a<-1) a=-1;
	if (a>1) a=1;
	return acos(a);
}

inline double mysqrt(double a)
{
    ASSERT( a > -1e-5 );
    double b=a<0?0:a;
    return sqrt(b);
}

//Leon found that having extern void unsetenv(...) caused
//problems on his linux box because unsetenv was defined elsewhere
//as having return type int
//#if !defined(__linux__) || !defined(powerpc)
//extern int unsetenv(const char* name);
//#endif
#if !defined(powerpc)
extern int setenv(const char* name, const char* value,
		  int overwrite); // 0=success
#endif

/// Multi-line header, each line begins with "# " and ends with "\n"
extern const char* CreateHeader(int argc, char** argv);

//@}
#endif


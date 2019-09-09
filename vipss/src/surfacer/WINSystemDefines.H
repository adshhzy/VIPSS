#ifndef _PLATFORM_DEPENDENT_DEFS_H
#define _PLATFORM_DEPENDENT_DEFS_H

/** \defgroup SystemDefines Platform dependent nonsense.

   include/WINSystemDefines.H is an amazingly ugly file containing all of the
   definitions/includes/workarounds needed to get my source code to
   compile on the SGI (CC compiler) sun (gcc or native compiler), OS10 and
   PCs (linux gcc or windows). Some things to note: <br>

   I don't use min/max/bool, instead I use WINmin, WINmax, and
   WINbool, which are my own definitions (or default to the
   appropriate definition on the given architecture). <br>

   I include all the system files I might need here
   (string/file/math/stdio/memory). This slows down compile time a
   bit, but it sure beats chasing down all of those include files when
   you try to get it to compile on another platform. <br>

   There's some syntactic differences between the compilers, the worst
   being those compilers which define for(int i ..) with the int i
   defined inside vs. outside the for loop. The solution is that calls
   after the initial for (int i..) should start with for ( FORINT i ).
   
   VERIFY and ASSERT are also a bit different on different platforms. <br>

   And there's the weird list of functions which appear on some
   platforms but not others. <br>
   @{
   */


#if defined(WIN32) || defined(WIN64)
typedef __int64 LONGLONG;
#else
typedef long long LONGLONG;
#endif

#if defined(WIN32) || defined(WIN64)
#ifdef WIN32
#define __WIN32
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#endif
#ifdef WIN64
#define __WIN64
//#define strncpy strncpy_s
#define sscanf sscanf_s
//#define strcat strcat_s
#endif

  #define NOMINMAX
  #define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers

  #ifdef _AFXDLL
    #include <afx.h>
    #include <afxwin.h>
  #else
    #include <assert.h>
    #include <crtdbg.h>
    #define ASSERT _ASSERTE
    #define TRACE printf
    #define VERIFY(x) if ( !(x) ) ASSERT(false);
    #ifndef FALSE
      #define FALSE false
      #define TRUE true
    #endif
  #endif
  #include <string>
  #include <iostream>
  #include <iomanip>
  #include <fstream>
  #include <memory>             // memcpy()
  #include <float.h>
  #include <math.h>
  using namespace std;

#ifndef isnan
	#define isnan _isnan
#endif
#ifndef isinf
	#define isinf !_finite
#endif


  // I have not thought this through.  Look at WSAEWOULDBLOCK...
  #define EWOULDBLOCK EAGAIN

  #define EmptyTemplate template<>
  #define tinline inline
  #define tempbrackets

/** \brief This gets around the difference in scoping between VC6.0 and .net. 
 Use for all loops after the first one. Usage:

 for ( int i; ... blah ) <br>
 ... <br>
 for ( FORINT i ... blah )

 Second (and subsequent) FORINTs default to either int (new compilers) or
 nothing (old compilers) */
#if _MSC_VER <= 1200
  #define FORINT
  #define FORSIZET
#else
  #define FORINT int
  #define FORSIZET size_t
#endif

  #ifndef assert
    #define assert ASSERT
  #endif

  typedef void (__cdecl *HHSIG_PF)(int sig,...);
  // bcopy() not available.
  inline void bcopy(const void* src, void* dst, int nbytes)
  {
      memcpy(dst,src,nbytes);
  }
  // sleep() available as Sleep()
  extern "C" __declspec(dllimport) void __stdcall Sleep(unsigned long millsecs);
  inline void sleep(int sec) { Sleep(sec*1000); }

  /// My version of a boolean.
  typedef bool WINbool;

#include <direct.h>
#include <sys/stat.h>

inline bool FileExists( const string &in_str )
{
    struct _stat status;
    int bExist = (_stat (in_str.c_str(), &status) != -1 );
    if(bExist){
        return true;
    } else {
        return false;
    }
}
#endif


#if defined(__linux__)
#  include <string.h>
#  include <strings.h>
#endif

#if defined(__linux__) || defined(powerpc) || defined(__APPLE__)
#  include <float.h>
#  include <stdio.h>
#  include <stdlib.h>
#  include <unistd.h>
#  include <netinet/in.h>
#  include <cmath>
#  include <climits>
#  include <fstream>
#  include <iostream>
#  include <string>
#  include <sys/stat.h>

using namespace std;

#  define UINT    int

#  define EmptyTemplate template<>
#  define tinline 
#  define tempbrackets <>

// This really should be the otherway 'round in windows, but 
// that doesn't seem to work 
#ifndef _isnan
        #define _isnan isnan
#endif
#ifndef _isinf
        #define _isinf isinf
#endif
#ifndef _finite
#define _finite !isinf
#endif

#  ifndef ASSERT
#    include <assert.h>
#    ifndef NDEBUG
#      define ASSERT(f) assert(f)
#      ifndef DEBUG
#        define DEBUG
#      endif
#    else
#      define ASSERT(f)
#    endif
#  endif



/** \brief This gets around the difference in scoping between VC6.0 and .net. 
 Use for all loops after the first one. Usage:

 for ( int i; ... blah ) <br>
 ... <br>
 for ( FORINT i ... blah )

 Second (and subsequent) FORINTs default to either int (new compilers) or
 nothing (old compilers) */
#  define FORINT int

#  ifdef DEBUG
#    define VERIFY(f) assert(f)
#  else
#    define VERIFY(f) f
#  endif

#  define TRACE printf

  extern const char* CTime();


  typedef void (*HHSIG_PF)(int);
  typedef bool WINbool;

#  ifndef FALSE
#    define FALSE 0
#  endif
#  ifndef TRUE
#    define TRUE 1
#  endif

inline bool FileExists( const string &in_str ) 
{
    struct stat status;
    stat( in_str.c_str(), &status );
    
    if ( S_ISREG( status.st_mode ) ) {
        return true;
    } else {
        return false;
    }
}

inline bool DirectoryExists( const string &in_str ) 
{
    struct stat status;
    stat( in_str.c_str(), &status );
    
    if ( access( in_str.c_str(), 0 ) == 0 ) {
        if ( status.st_mode & S_IFDIR ) {
            return true;
        } else {
            return false;
        }
    }
    return false;
}

inline bool Do_mkdir( const string &in_str, const mode_t mode = S_IRWXU ) 
{
    if ( !DirectoryExists( in_str ) ) {
        /* Directory does not exist */
        if (mkdir( in_str.c_str(), mode) != 0)
            return false;
    }
    
    return true;
}

inline bool Do_mkpath( const string &in_str, const mode_t mode = S_IRWXU ) 
{
    size_t iStart = 1;
    size_t iSplit = in_str.find_first_of( string("/\\"), iStart );
    string strDir;
    while ( iSplit != string::npos ) {
        strDir = in_str;
        strDir.erase( iSplit, string::npos );
        if ( !Do_mkdir( strDir, mode ) ) {
            return false;
        }
        iStart = iSplit+1;
        iSplit = in_str.find_first_of( string("/\\"), iStart );        
    }
    
    return true;
}


#endif

/// A kinda clever way to get around undeclared variable warnings in release mode
  /// (use for, e.g., intput index variables that are not touched/used in release mode)
#if defined(DEBUG) || defined(_DEBUG)
#define DECLAREVAR(v) v
#else
#define DECLAREVAR(v)
#endif


/// Helper functions to read/write booleans in binary
inline WINbool WINBoolRead( istream &in )
{
    char c;
    in >> c;
    if ( c == 't' ) return TRUE;
    else if ( c == 'f' ) return FALSE;
    else ASSERT(FALSE);

    return FALSE;
}

/// Helper functions to read/write booleans in binary
inline void WINBoolWrite( ostream &out, const WINbool &in_b )
{
    if ( in_b == FALSE ) {
        out << " f ";
    } else {
        out << " t ";
    }
}

/** @} */

/** \defgroup SystemDefinesIO Open/Close files 
  \ingroup SystemDefines

 * Open files in binary mode. Checks result.  
 * @{
 */
///
inline void WINWriteBinary(const char *name, ofstream &out)
{ 
    out.clear();

#ifdef sgi
    out.open(name, ios::out); 
#else
    out.open(name, ios::out | ios::binary); 
#endif
    if (!out.good()) {
        cerr << "Bad file " << name << "\n";
        ASSERT(FALSE); 
    }
}
/// Returns TRUE if file was created
inline WINbool WINWriteBinaryTest(const char *name, ofstream &out)
{ 
    out.clear();

#ifdef sgi
    out.open(name, ios::out); 
#else
    out.open(name, ios::out | ios::binary); 
#endif
    if (!out.good()) {
        cerr << "Bad file " << name << "\n";
        return FALSE; 
    }
    return TRUE;
}


///
inline void WINAppendBinary(const char *name, ofstream &out)
{ 
    out.clear();
#ifdef sgi
    out.open(name, ios::app); 
#else
    out.open(name, ios::app | ios::binary); 
#endif

    if (!out.good()) {
        cerr << "Bad file " << name << "\n";
        ASSERT(FALSE); 
    }
}
///
inline void WINReadBinary(const char *name, ifstream &in)
{ 
    in.clear();
#ifdef sgi
    in.open(name, ios::in); 
#else
    in.open(name, ios::in | ios::binary); 
#endif

    if (!in.good()) {
        cerr << "Bad file " << name << "\n";
        ASSERT(FALSE); 
    }
}
/// Returns true if file exists
inline WINbool WINReadBinaryTest(const char *name, ifstream &in)
{ 
    in.clear();
#ifdef sgi
    in.open(name, ios::in); 
#else
    in.open(name, ios::in | ios::binary); 
#endif

    if (!in.good()) {
        cerr << "Bad file " << name << "\n";
        return FALSE; 
    }
    return TRUE;
}
/** @} */

/** \defgroup SystemDefinesSTD Min/max, bool.

  \ingroup SystemDefines

  I'm sick of trying to work around everyone's personal definitions
   of bool, min, and max. So I made my own :).
   WINbool is an enumerated type on suns/sgi, type bool on windows. 
  * @{
  */
///
template<class T> inline T WINmin(T a, T b) { return a<b?a:b; }
///
template<class T> inline T WINmax(T a, T b) { return a>b?a:b; }
///
template<class T> inline T WINminmax(T a, T tMin, T tMax) 
{ 
    return (a > tMax) ? tMax : ( (a < tMin) ? tMin : a ); 
}

/** @} */

#endif





#ifndef _UTILS_IO_DEFS_H
#define _UTILS_IO_DEFS_H

#include <utils/Rn_Defs.H>

/** \class UTILSIo Utils_Io.H utils/Utils_Io.H
  \ingroup UtilitiesFunc
  \brief For reading and writing in binary accross platforms.

   This could be templated, but there are different things to do for floats vs. ints.
   Little indian vs. Big indian. 

  \sa UTILSIoWriteBinaryPtVec, UTILSIoReadBinaryPtVec
<br>   Files: 
- includeutils/Utils_Io.H
- utils/utils/UTILSIo.cpp
*/
class UTILSIo {
private:
    static long    Convert( int            in_i );
    static long    Convert( short          in_i );
    static long    Convert( long           in_l );
    static long    Convert( float          in_f );
    static long    Convert( double         in_d );
    static long    Convert( unsigned short in_us );
    static long    Convert( unsigned long int   in_ui );

    static void ConvertBack( int            &out_i,  long in_l );
    static void ConvertBack( short          &out_i,  long in_l );
    static void ConvertBack( long           &out_l,  long in_l );
    static void ConvertBack( float          &out_f,  long in_l );
    static void ConvertBack( double         &out_d,  long in_l );
    static void ConvertBack( unsigned short &out_ui, long in_l );
    static void ConvertBack( unsigned long int   &out_us, long in_l );

public:
    /** Number of digits to save. Make sure it's the same when converting back.
        @param in_iNDigits 30 works for +-1.0, 29 for +-2.0, etc. Max is 30.*/
    static void    SetNDigits( int in_iNDigits );

    /**@name Single values (double, float, int) */
    //@{

    ///
    static void WriteBinary( int in_i, ofstream &);
    /// 
    static WINbool ReadBinary( int &io_i, ifstream &);
    ///
    static void WriteBinary( float in_f, ofstream &);
    ///
    static void WriteBinary(  short in_is, ofstream &);
    /// 
    static WINbool ReadBinary(  short &io_is, ifstream &);
    /// 
    static WINbool ReadBinary( float &io_f, ifstream &);
    ///
    static void WriteBinary( double in_d, ofstream &);
    /// 
    static WINbool ReadBinary( double &io_d, ifstream &);
    ///
    static void WriteBinary( unsigned short in_is, ofstream &);
    /// 
    static WINbool ReadBinary( unsigned short &io_is, ifstream &);
    //@}

    /**@name Arrays. 
       Find the best value to use then write out the array. ofstream must be opened binary.
        Returns TRUE if numbers were not too big */
    //@{
    ///
    static WINbool WriteBinary(const Array<char> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<int> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<short> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<long> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<float> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<double> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<unsigned long int> &, ofstream &);
    ///
    static WINbool WriteBinary(const Array<unsigned short> &, ofstream &);
    /// 
    static WINbool ReadBinary( Array<char> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<int> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<short> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<long> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<float> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<double> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<unsigned long int> &, ifstream &);
    /// 
    static WINbool ReadBinary( Array<unsigned short> &, ifstream &);
    //@}

    /// Test
    static void    TestConversion();
};

/** \ingroup HelperGeometry

    Only use these if you need to write/read cross platform

   Points and vectors ( double, float ). 
   Find the best value to use then write out the array. 
   ofstream must be opened binary.
   Returns TRUE if numbers were not too big. */
template<class coord, class affine>
inline WINbool UTILSIoWriteBinaryPtVec(const Array<affine> & in_ao, 
                                       ofstream            & out, 
                                       coord                 in_dummy)
{
    int iDim = ( in_ao.num() ) ? in_ao[0].Dim() : 0;
    
    Array<coord> ad( in_ao.num() * iDim );
    
    int iIndx = 0;
    for (int i = 0; i < in_ao.num(); i++) {
        for (int iD = 0; iD < iDim; iD++) {
            ad[iIndx++] = (coord) in_ao[i][iD];
        }
    }
    return UTILSIo::WriteBinary( ad, out );
}

/** \ingroup HelperGeometry

    Only use these if you need to write/read cross platform

   Points and vectors ( double, float ). 
   Find the best value to use then write out the array. 
   ofstream must be opened binary.
   Returns TRUE if numbers were not too big. */
template<class coord, class affine>
inline WINbool UTILSIoReadBinaryPtVec( Array<affine> &io_ao, ifstream &in, coord in_dummy)
{
    Array<coord> ad;

    WINbool bRes = UTILSIo::ReadBinary( ad, in );

    if ( ad.num() == 0 ) {
        io_ao.need(0);
        return bRes;
    }

    io_ao.need(1);
    int iDim = io_ao[0].Dim();
    
    if ( ad.num() % iDim != 0 ) bRes = FALSE;
    
    io_ao.need( ad.num() / iDim );
    
    int iIndx = 0;
    for (int i = 0; i < io_ao.num(); i++) {
        for (int iD = 0; iD < iDim; iD++) {
            io_ao[i][iD] = ad[iIndx++];
        }
    }

    return bRes;
}


#endif

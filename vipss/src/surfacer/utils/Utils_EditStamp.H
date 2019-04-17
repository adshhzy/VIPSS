#ifndef _UTILS_DEFS_EDIT_STAMP_H
#define _UTILS_DEFS_EDIT_STAMP_H


#include <WINSystemDefines.H>

/** \class UTILSEditStamp Utils_EditStamp.H utils/Utils_EditStamp.H
  \ingroup UtilitiesFunc
  \brief A way to do time stamping

  I'm sure there's a boost or std version of this; essentially a time stamp. Everytime
  Edit is called, the time stamp id is incremented by one. Makes lazy/delayed evaluation
  possible in the UniformMflds code.

  Both the class that's changing and the class that needs to know about the first class
  need to have an instance. The changing class calls Edit() whenever something changes. The
  second class can check to see if it's version is the same by keeping an instance around
  and using OutOfDate and Update to check the current status wrt the changing class.
  */
class UTILSEditStamp {
protected:
    int m_iStamp;

public:
    /// 
    void    Edit() {  m_iStamp++; }
    ///
    WINbool OutOfDate( const UTILSEditStamp &in_i ) const { return in_i.m_iStamp == m_iStamp ? FALSE : TRUE; }
    ///
    void Update(  const UTILSEditStamp &in_i ) { m_iStamp = in_i.m_iStamp; }
    /// For reading in new data
    void Reset(   ) { m_iStamp = 0; }

    ///
    UTILSEditStamp() : m_iStamp(0) {}
    ///
    virtual ~UTILSEditStamp() { m_iStamp = -1; }
};

#endif

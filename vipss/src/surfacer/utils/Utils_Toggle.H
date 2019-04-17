#ifndef _UTILS_TOGGLE_DEFS_H
#define _UTILS_TOGGLE_DEFS_H

#include <WINSystemDefines.H>

/** Toggle between true and false.

  \ingroup UtilitiesFunc

  \brief Basically a boolean, but forces initialization so's you don't forget to :)

<br> Files:
- include/utils/Utils_Toggle.H

  */
class UTILSToggle  {
private:
    WINbool m_bBool;

public:
    ///
    WINbool Toggle()
    { 
        m_bBool = ( m_bBool == TRUE ? FALSE : TRUE );
        return m_bBool;
    }

    ///
    WINbool operator=( WINbool in_b ) { m_bBool = in_b; return m_bBool; }
    ///
    WINbool operator==( WINbool in_b ) const
    {
        if ( in_b == m_bBool ) return TRUE;
        return FALSE;
    }

    /// For those functions which reallllly want an int
    int Int() const  { return m_bBool; }

    ///
    WINbool operator()() const { return m_bBool; }

    ///
    UTILSToggle( WINbool in_bSet ) { m_bBool = in_bSet; }
    ///
    ~UTILSToggle() {}
};

#endif

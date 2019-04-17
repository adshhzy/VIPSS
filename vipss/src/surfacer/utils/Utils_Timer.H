// timer.h

#include <sys/types.h>
#include <sys/timeb.h>

#ifdef SGI
#include <time.h>
#endif

/**@name Count seconds and milliseconds
   \ingroup UtilitiesFunc

   May not work under unix system. <br>
   Bracked area to be times with UTILSTickStart and UTILSTickEnd, 
   then call UTILSSecondsFromTick to get time or print out result from UTILSTickEnd. 
<br>   FILES: <br>
- include/utils/Utils_Timer.H
- utils/utils/Timer.cpp
*/
//@{

///
typedef struct {long milliseconds; long seconds;} UTILSTick;
///
UTILSTick UTILSTickStart(void);
/// Call with tick returned from UTILSTickStart.
UTILSTick UTILSTickEnd(UTILSTick tickStart);
/// 
double UTILSSecondsFromTick(UTILSTick tick);

//@}

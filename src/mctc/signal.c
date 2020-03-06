/*
 *  --- taken from the internet ---
 *
 *  It overrides the FORTRAN intrinsic signal routine
 *  and gives the possibility to use custom signal
 *  handlers instead of the standard ones, pretty simple
 *  code but actually impossible to do with FORTRAN only.
 */
#include <signal.h>

typedef void (*sighandler_t)(int);

void signal_( int* signum, sighandler_t handler)
{
   signal(*signum, handler);
}

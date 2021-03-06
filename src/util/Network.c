 /*@@
   @file      Network.c
   @date      Fri 19 Jan 2001
   @author    Thomas Radke
   @desc
              Network related routines
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"
#include "util_Network.h"

#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#ifdef HAVE_NETDB_H
#include <netdb.h>
#elif defined HAVE_WINSOCK2_H
#include <winsock2.h>
#endif /* HAVE_WINSOCK2_H */

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(util_Network_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Util_GetHostName
   @date       Fri Oct 20 12:12:34 2000
   @author     Tom Goodale
   @desc
               Gets the fully qualified name of this host if possible.
   @enddesc
   @calls
   @history
               just copied from thorn HTTPD
   @endhistory
   @var     name
   @vdesc   character buffer to store name in
   @vtype   char *
   @vio     out
   @vcomment

   @endvar
   @var     length
   @vdesc   length of the character buffer
   @vtype   int
   @vio     in
   @vcomment

   @endvar

@@*/
void Util_GetHostName (char *name, int length)
{
  gethostname (name, length);

  /* Does the name include the domain. */
  if (! strchr (name, '.'))
  {
#ifdef HAVE_GETHOSTBYNAME
    struct hostent *thishostent=0;
#ifndef CRAY_XT
    thishostent = gethostbyname (name);
#endif

    if (thishostent)
    {
      strncpy (name, thishostent->h_name, length);
      name[length - 1] = 0;
    }
    else name[0]='\0';
#endif
  }
}

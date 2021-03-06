 /*@@
   @file      SetupCache.c
   @date      Tue Nov 30 10:30:09 1999
   @author    Tom Goodale
   @desc
              Sets up cache stuff for the CCTK
   @enddesc
   @version   $Id: SetupCache.c 3791 2004-06-26 11:16:40Z tradke $
 @@*/

#include "cctk_Config.h"
#include "cctk_Flesh.h"
#include "cctk_Parameters.h"

#include "cctki_Cache.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_SetupCache_c);


 /*@@
   @routine    CCTKi_SetupCache
   @date       Tue Nov 30 10:50:02 1999
   @author     Tom Goodale
   @desc
               Sets the cache information.
   @enddesc
   @calls      Utili_CacheDataSet

   @returntype int
   @returndesc
               0 - success
   @endreturndesc
@@*/
int CCTKi_SetupCache (void)
{
  DECLARE_CCTK_PARAMETERS


  if (manual_cache_setup)
  {
    Utili_CacheDataSet (manual_cacheline_bytes, manual_cache_size);
  }
  else
  {
    Utili_CacheDataSet (CCTK_CACHELINE_BYTES, CCTK_CACHE_SIZE);
  }

  return (0);
}

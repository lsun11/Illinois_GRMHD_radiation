 /*@@
   @file      ShutdownCactus.c
   @date      Mon Sep 28 14:25:48 1998
   @author    Tom Goodale
   @desc
              Contains routines to shutdown cactus.
   @enddesc
   @version   $Id: ShutdownCactus.c 4467 2008-03-01 15:25:31Z schnetter $
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk_Flesh.h"
#include "cctk_Misc.h"
#include "cctk_Parameters.h"
#include "cctk_Schedule.h"

#ifdef CCTK_MPI
#include "mpi.h"
#endif

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_ShutdownCactus_c)

#ifdef CCTK_MPI
extern char MPI_Active;
#endif

#ifdef CCTK_MPI
#define CACTUS_MPI_ERROR(xf)                                                  \
          do                                                                  \
          {                                                                   \
            int errcode;                                                      \
                                                                              \
                                                                              \
            if((errcode = xf) != MPI_SUCCESS)                                 \
            {                                                                 \
              char mpi_error_string[MPI_MAX_ERROR_STRING+1];                  \
              int resultlen;                                                  \
                                                                              \
                                                                              \
              MPI_Error_string(errcode, mpi_error_string, &resultlen);        \
              fprintf(stderr, "MPI Call %s returned error code %d (%s)\n",    \
                              #xf, errcode, mpi_error_string);                \
              fprintf(stderr, "At line %d of file %s\n",                      \
                              __LINE__, __FILE__);                            \
            }                                                                 \
          } while (0)
#endif

 /*@@
   @routine    CCTKi_ShutdownCactus
   @date       Mon Sep 28 14:50:50 1998
   @author     Tom Goodale
   @desc
               Cactus specific shutdown stuff.
   @enddesc
   @calls      CCTK_SchedulePrintTimes

   @var        ConfigData
   @vdesc      Flesh configuration data
   @vtype      tFleshConfig *
   @vio        unused
   @endvar

   @returntype int
   @returndesc
               0  - success
   @endreturndesc
@@*/
int CCTKi_ShutdownCactus(tFleshConfig *ConfigData)
{
  DECLARE_CCTK_PARAMETERS


  /* avoid compiler warning about unused argument */
  ConfigData = ConfigData;

  if (CCTK_Equals (cctk_timer_output, "full"))
  {
    CCTK_SchedulePrintTimes (NULL);
  }

  printf("--------------------------------------------------------------------------------\n"); 
  printf("Done.\n");
  fflush(stdout);
 
#ifdef CCTK_MPI
  if(MPI_Active)
  {
    CACTUS_MPI_ERROR(MPI_Finalize());
  }
#endif

  return 0;
}

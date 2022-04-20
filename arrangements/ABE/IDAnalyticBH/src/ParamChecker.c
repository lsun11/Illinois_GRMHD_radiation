 /*@@
   @file      ParamChecker.F
   @date      March 1999
   @author    Gabrielle Allen
   @desc 
      Check parameters for black hole initial data and give some
      information
   @enddesc 
   @version $Header: /cactus/CactusEinstein/IDAnalyticBH/src/ParamChecker.c,v 1.17 2004/05/17 12:28:57 goodale Exp $
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/IDAnalyticBH/src/ParamChecker.c,v 1.17 2004/05/17 12:28:57 goodale Exp $";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_ParamChecker_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
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
   @routine    IDAnalytic_ParamChecker
   @date       March 1999
   @author     Gabrielle Allen
   @desc 
      Check parameters for black hole initial data and give some
      information
   @enddesc 
   @calls     
   @history 
   @hdate Fri Apr 26 09:48:24 2002 @hauthor Tom Goodale
   @hdesc Modified to new ADMBase and StaticConformal stuff 
   @endhistory 

@@*/

void IDAnalyticBH_ParamChecker(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  char *message;

  /* Do we know how to deal with this type of metric ? */
  if(! CCTK_EQUALS(metric_type, "physical") &&
     ! CCTK_EQUALS(metric_type, "static conformal"))
  {
    CCTK_PARAMWARN("Unknown ADMBase::metric_type - known types are \"physical\" and \"static conformal\"");
  }


  if (CCTK_Equals(initial_data,"schwarzschild") == 1)
  {
    CCTK_INFO("Schwarzschild black hole");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  throat at %f",mass/2.0);
    CCTK_INFO(message);
    free(message);
  }
  if (CCTK_Equals(initial_data,"kerr") == 1)
  {
    CCTK_REAL tmp;
    CCTK_INFO("Kerr black hole");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  m = %f, a = %f",mass,a_Kerr);
    CCTK_INFO(message);
    tmp=pow((mass*mass-a_Kerr*a_Kerr)/4.0,.5);
    sprintf(message,"  throat at %f",tmp);
    CCTK_INFO(message);
    free(message);
  }
  else if (CCTK_Equals(initial_data,"bl_bh") == 1)
  {
    CCTK_INFO("Brill Lindquist black holes");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  %d black holes",bl_nbh);
    CCTK_INFO(message);
    if (bl_nbh > 0)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
              bl_M_1,bl_x0_1,bl_y0_1,bl_z0_1);
      CCTK_INFO(message);
    }
    if (bl_nbh > 1)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
              bl_M_2,bl_x0_2,bl_y0_2,bl_z0_2);
      CCTK_INFO(message);
    }
    if (bl_nbh > 2)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
              bl_M_3,bl_x0_3,bl_y0_3,bl_z0_3);
      CCTK_INFO(message);
    }
    if (bl_nbh > 3)
    {
      sprintf(message,  "  mass %f at (%f,%f,%f)",
              bl_M_4,bl_x0_4,bl_y0_4,bl_z0_4);
      CCTK_INFO(message);
    }
    free(message);
  }
  else if (CCTK_Equals(initial_data,"multiple_misner") == 1)
  {
    CCTK_INFO("Setting up Misner solution for multiple holes");
  }
  else if (CCTK_Equals(initial_data,"misner_bh")==1)
  {
    CCTK_INFO("Two Misner black holes (on z-axis)");
    message = (char *)malloc(200*sizeof(char));
    sprintf(message,"  mu is %f",mu);
    CCTK_INFO(message);
    free(message);
  }

  /*     Remind users about the conformal metric
   *     ---------------------------------------
   */
  
  if (CCTK_EQUALS(metric_type, "static conformal"))
  {
    CCTK_INFO("Black hole initial data uses conformal metric");
    if (CCTK_EQUALS(conformal_storage,"factor+derivs") ||
        CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs"))
    {
      CCTK_INFO("  and conformal derivatives");
    }
    else
    {
      CCTK_INFO("  but no conformal derivatives");
    }
  }
  else
  { 
    CCTK_INFO("Implements non-conformal metric");
    CCTK_INFO("  (Not usually a good idea!)");
  }        

}

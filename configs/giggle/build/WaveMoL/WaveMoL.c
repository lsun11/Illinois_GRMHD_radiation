 /*@@
   @file      WaveMoL.c
   @date      Fri Nov  9 13:33:25 2001
   @author    Ian Hawke
   @desc 
   The equivalent of WaveToy.c
   @enddesc 
 @@*/

  
#include "cctk.h" 
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/*
#ifndef DEBUG_MOL
#define DEBUG_MOL
#endif
*/

static const char *rcsid = "$Header: /cactusdevcvs/CactusExamples/WaveMoL/src/WaveMoL.c,v 1.2 2004/08/19 09:25:48 schnetter Exp $";

CCTK_FILEVERSION(CactusExamples_WaveMoL_WaveMoL_c)

void WaveMoL_Boundaries(CCTK_ARGUMENTS);
void WaveMoL_CalcRHS(CCTK_ARGUMENTS);
void WaveMoL_Energy(CCTK_ARGUMENTS);

 /*@@
   @routine    WaveMoL_Boundaries
   @date       Fri Nov  9 13:34:47 2001
   @author     Ian Hawke
   @desc 
   A slight alteration of the standard routine. Note that I 
   haven't set up the error checking correctly as yet.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void WaveMoL_Boundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int ierr=-1;
  int sw[3];

#ifdef DEBUG_MOL
  printf("We've reached the Boundary enforcement; step must be done!\n");
#endif

  /* Set the stencil width */
  sw[0]=1;
  sw[1]=1;
  sw[2]=1;

  /* Uses all default arguments, so invalid table handle -1 can be passed */
  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "wavemol::scalarevolvemol_scalar", bound);

  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }

  ierr = Boundary_SelectGroupForBC
    (cctkGH, CCTK_ALL_FACES, 1, -1, "wavemol::scalarevolvemol_vector", bound);

  if (ierr < 0) 
  {
    CCTK_WARN(0,"Boundary conditions not applied - giving up!");
  }

  return;
}

 /*@@
   @routine    WaveMoL_CalcRHS
   @date       Fri Nov  9 13:37:27 2001
   @author     Ian Hawke
   @desc 
   The routine that calculates the RHS of the PDE for the MoL step.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void  WaveMoL_CalcRHS(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
      
  int i,j,k;
  int index;
  int istart, jstart, kstart, iend, jend, kend;
  CCTK_REAL dx,dy,dz;
  CCTK_REAL hdxi, hdyi, hdzi;

  /* Set up shorthands */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);

  hdxi = 0.5 / dx;
  hdyi = 0.5 / dy;
  hdzi = 0.5 / dz;
  
  istart = 1;
  jstart = 1;
  kstart = 1;
  
  iend = cctk_lsh[0]-1;
  jend = cctk_lsh[1]-1;
  kend = cctk_lsh[2]-1;

  /* Calculate the right hand sides. */

#ifdef DEBUG_MOL
  printf("About to loop.\n");
#endif

  for (k=0; k<cctk_lsh[2]; k++)
  {
    for (j=0; j<cctk_lsh[1]; j++)
    {
      for (i=0; i<cctk_lsh[0]; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        phirhs[index] = phit[index];
        phitrhs[index] = 0;
        phixrhs[index] = 0;
        phiyrhs[index] = 0;
        phizrhs[index] = 0;
      }
    }
  }  

  for (k=kstart; k<kend; k++)
  {
    for (j=jstart; j<jend; j++)
    {
      for (i=istart; i<iend; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        
        phitrhs[index] = 
            (phix[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] -
             phix[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]) *hdxi
          + (phiy[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] -
             phiy[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]) *hdyi
          + (phiz[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] -
             phiz[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]) *hdzi;
        phixrhs[index] = (phit[CCTK_GFINDEX3D(cctkGH, i+1, j, k)] -
                          phit[CCTK_GFINDEX3D(cctkGH, i-1, j, k)]) *hdxi;
        phiyrhs[index] = (phit[CCTK_GFINDEX3D(cctkGH, i, j+1, k)] -
                          phit[CCTK_GFINDEX3D(cctkGH, i, j-1, k)]) *hdyi;
        phizrhs[index] = (phit[CCTK_GFINDEX3D(cctkGH, i, j, k+1)] -
                          phit[CCTK_GFINDEX3D(cctkGH, i, j, k-1)]) *hdzi;

#ifdef DEBUG_MOL
        printf("ijk %i %i %i, phitrhs %f\n",i,j,k,phitrhs[index]);
#endif  
      } 
    }
  }

#ifdef DEBUG_MOL
  printf("Done with loop\n");
#endif

  return;
}

 /*@@
   @routine    WaveMoL_Energy
   @date       Sat Jun 26 15:39:18 2004
   @author     Ian Hawke
   @desc 
   Computes the energy d_t phi ^2 + d_x^i phi ^2.
   An illustration of using constrained MoL variables rather
   than anything useful.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


void WaveMoL_Energy(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS
      
  int i,j,k;
  int index;

  for (k=0; k<cctk_lsh[2]; k++)
  {
    for (j=0; j<cctk_lsh[1]; j++)
    {
      for (i=0; i<cctk_lsh[0]; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        energy[index] = phit[index]*phit[index] +
          phix[index]*phix[index] +
          phiy[index]*phiy[index] +
          phiz[index]*phiz[index];
      }
    }
  }
  
}

          

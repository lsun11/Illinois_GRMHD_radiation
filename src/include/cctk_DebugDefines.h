/*@@
   @file      cctk_DebugDefines.h
   @date      Sun 28 Dec 2003
   @author    Erik Schnetter
   @desc
              Routines to provide some debugging support for the Cactus code.
   @enddesc
   @version   $Id: cctk_DebugDefines.h 3536 2004-01-19 21:39:52Z schnetter $
 @@*/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

int CCTK_GFIndex1D (const cGH *GH, int i);
int CCTK_GFIndex2D (const cGH *GH, int i, int j);
int CCTK_GFIndex3D (const cGH *GH, int i, int j, int k);
int CCTK_GFIndex4D (const cGH *GH, int i, int j, int k, int l);

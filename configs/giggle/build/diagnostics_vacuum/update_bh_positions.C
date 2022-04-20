//--------------------------------------------------------------------------
// Set Xbh1,Ybh1,Zbh1; and Xbh2,Ybh2,Zbh2 if desired.
//--------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <stddef.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

extern "C" void update_bh_posns(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  
  int horizon_number,foundflag;
  double xx,yy,zz,ahradius;
  
  horizon_number = 1;
  //WARNING: THE FOLLOWING FUNCTION IS TOTALLY BROKEN AFTER CHECKPOINT.  It will say there is a horizon, but only give coordinates 0,0,0!  Thus we now use ah_centroid_x() values.
  foundflag = HorizonCentroid(horizon_number,&xx,&yy,&zz);

  //###########################################################################
  //# HorizonRadiusInDirection() computes the horizon radius in the direction
  //# of each (x,y,z) point, or -1.0 if this horizon wasn't found the most
  //# recent time AHFinderDirect searched for it.  More precisely, for each
  //# (x,y,z), consider the ray from the local coordinate origin through
  //# (x,y,z).  This function computes the Euclidean distance between the
  //# local coordinate origin and this ray's intersection with the horizon,
  //# or -1.0 if this horizon wasn't found the most recent time AHFinderDirect
  //# searched for it.  
  //#

  int N_points = 1;
  double horizdirn_x = 0.0;
  double horizdirn_y = 0.0;
  double horizdirn_z = 1000000.0;
  int dummy = HorizonRadiusInDirection(horizon_number,N_points,&horizdirn_x, &horizdirn_y, &horizdirn_z,&ahradius);
    
  //if(foundflag!=1) xx=yy=zz=0.0;
  if(foundflag!=1) {
    //The following line is useful if we had a checkpoint in which bh_radius_z was not set.  Otherwise, this line is harmless.
    //if(bh_radius_z[0] <= 0.0) bh_radius_z[0] = ahradius;

    printf("WARNING: HORIZON #1 NOT FOUND.  Using last known BH position for diagnostic integrals, which is: %e %e %e; radius = %e\n",bh_posn_x[0],bh_posn_y[0],bh_posn_z[0],bh_radius_z[0]);

    /*
    bh_posn_x[0] = ah_centroid_x[0];
    bh_posn_y[0] = ah_centroid_y[0];
    bh_posn_z[0] = ah_centroid_z[0];
    */
  } else {    
    bh_posn_x[0] = ah_centroid_x[0];
    bh_posn_y[0] = ah_centroid_y[0];
    bh_posn_z[0] = ah_centroid_z[0];
    /*
    bh_posn_x[0] = xx;
    bh_posn_y[0] = yy;
    bh_posn_z[0] = zz;
    */
    bh_radius_z[0] = ahradius;
    printf("HORIZON #1 POSITION: %e %e %e, RADIUS= %e\n",bh_posn_x[0],bh_posn_y[0],bh_posn_z[0],ahradius);
  }

  if(num_BHs>1) {
    horizon_number = 2;
    foundflag = HorizonRadiusInDirection(horizon_number,N_points,&horizdirn_x, &horizdirn_y, &horizdirn_z,&ahradius);
    //WARNING: THE FOLLOWING FUNCTION IS TOTALLY BROKEN AFTER CHECKPOINT.  It will say there is a horizon, but only give coordinates 0,0,0!  Thus we now use ah_centroid_x() values.
    foundflag = HorizonCentroid(horizon_number,&xx,&yy,&zz);
    if(foundflag==1) {
      bh_posn_x[1] = ah_centroid_x[1];
      bh_posn_y[1] = ah_centroid_y[1];
      bh_posn_z[1] = ah_centroid_z[1];
      bh_radius_z[1] = ahradius;
    }
    printf("HORIZON #2 POSITION: %e %e %e, RADIUS= %e\n",xx,yy,zz,ahradius);

    printf("Checking for 3rd (common) horizon!\n");
    for(int i=0;i<3;i++) {
      horizon_number = i+1;
      //WARNING: THE FOLLOWING FUNCTION IS TOTALLY BROKEN AFTER CHECKPOINT.  It will say there isa horizon, but only give coordinates 0,0,0!  Thus we now use ah_centroid_x() values.
      foundflag = HorizonCentroid(horizon_number,&xx,&yy,&zz);
      if(foundflag==1) {
	foundflag = HorizonRadiusInDirection(horizon_number,N_points,&horizdirn_x, &horizdirn_y, &horizdirn_z,&ahradius);
	bh_posn_x[i] = ah_centroid_x[i];
	bh_posn_y[i] = ah_centroid_y[i];
	bh_posn_z[i] = ah_centroid_z[i];
	bh_radius_z[i] = ahradius;
      }
      printf("HORIZON #3 POSITION: %e %e %e, RADIUS= %e\n",xx,yy,zz,ahradius);
    }
  }

}

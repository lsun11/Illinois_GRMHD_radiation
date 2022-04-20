/* $Header: /numrelcvs/AEIThorns/SphericalSurface/src/setup.c,v 1.8 2007/08/23 21:10:04 schnetter Exp $ */

#include <assert.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



void SphericalSurface_Setup (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL const pi = 3.1415926535897932384626433832795028841971693993751;
  
  int group;
  cGroup groupinfo;
  cGroupDynamicData groupdata;
  
  int n;
  int ierr;
  
  
  
  if (nsurfaces == 0) return;
  
  
  
  group = CCTK_GroupIndex ("SphericalSurface::sf_radius");
  assert (group>=0);
  ierr = CCTK_GroupData (group, &groupinfo);
  assert (!ierr);
  ierr = CCTK_GroupDynamicData (cctkGH, group, &groupdata);
  assert (!ierr);
  
  
  
  for (n=0; n<nsurfaces; ++n) {
    
    /* internal consistency checks */
    assert (groupdata.dim == 2);
    assert (groupdata.gsh[0] >= ntheta[n]);
    assert (groupdata.gsh[1] >= nphi[n]);
    
    assert (ntheta[n] >= 3*nghoststheta[n] && ntheta[n] <= maxntheta);
    assert (nphi[n] >= 3*nghostsphi[n] && nphi[n] <= maxnphi);
    
    
    
    /* copy parameters into grid functions */
    sf_ntheta[n] = ntheta[n];
    sf_nphi[n] = nphi[n];
    sf_nghoststheta[n] = nghoststheta[n];
    sf_nghostsphi[n] = nghostsphi[n];
    
    
    
    /* coordinates in the theta direction */
    /* avoid_sf_origin_theta = 1 */
    if (symmetric_z[n]) {
      
      /* upper hemisphere: z>=0, theta in (0, pi/2) */
      sf_delta_theta[n] = pi/2 / (ntheta[n] - 2*nghoststheta[n] - 0.5);
      sf_origin_theta[n] = - (nghoststheta[n] - 0.5) * sf_delta_theta[n];
      
    } else {
      
      /* both hemispheres: theta in (0, pi) */
      sf_delta_theta[n] = pi / (ntheta[n] - 2*nghoststheta[n]);
      sf_origin_theta[n] = - (nghoststheta[n] - 0.5) * sf_delta_theta[n];
      
    }
    
    
    
    /* coordinates in the phi direction */
    /* avoid_sf_origin_phi = 0 */
    if (symmetric_x[n]) {
      if (symmetric_y[n]) {
        
        /* one quadrant: x>=0, y>=0, phi in [0, pi/2] */
        assert (nphi[n] - 2*nghostsphi[n] >= 1);
        sf_delta_phi[n] = pi/2 / (nphi[n] - 2*nghostsphi[n] - 1);
        sf_origin_phi[n] = - nghostsphi[n] * sf_delta_phi[n];
        
      } else {
        
        /* two quadrants: x>=0, phi in [-pi/2, pi/2] */
        assert (nphi[n] - 2*nghostsphi[n] >= 2);
        sf_delta_phi[n] = pi / (nphi[n] - 2*nghostsphi[n] - 1);
        sf_origin_phi[n] = - pi/2 - nghostsphi[n] * sf_delta_phi[n];
        
      }
    } else {
      if (symmetric_y[n]) {
        
        /* two quadrants: y>=0, phi in [0, pi] */
        assert (nphi[n] - 2*nghostsphi[n] >= 2);
        sf_delta_phi[n] = pi / (nphi[n] - 2*nghostsphi[n] - 1);
        sf_origin_phi[n] = - nghostsphi[n] * sf_delta_phi[n];
        
      } else {
        
        /* all quadrants: phi in [0, 2pi) */
        assert (nphi[n] - 2*nghostsphi[n] >= 4);
        sf_delta_phi[n] = 2*pi / (nphi[n] - 2*nghostsphi[n]);
        sf_origin_phi[n] = - nghostsphi[n] * sf_delta_phi[n];
        
      }
    }
    
    
    
    /* mark surface as uninitialised */
    sf_active[n] = 0;
    //printf("RESETTING SURFACE #\%d TO ZERO!\n",n);
    sf_valid[n] = 0;
    
  } /* for n */
}

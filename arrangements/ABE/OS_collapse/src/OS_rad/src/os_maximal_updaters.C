#define NRANSI
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <string.h>
#include <fstream.h>
#include "nr.h"
#include "nrutil.h"
#include "Symmetry.h"
#define N_INT 1000

static char *rcsid = "$Meow...$";
CCTK_FILEVERSION(os_maximal_updaters)
extern "C" void os_maximal_updaters(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  //read in from file
  double *ls_rad_lo,*ls_rad_hi,*lapse_lo,*lapse_hi;
  double *shift_lo,*shift_hi;
  double *rho_lo,*rho_hi;
  //derivs obtained from spline
  double *ddlapse_lo,*ddlapse_hi;
  double *ddshift_lo,*ddshift_hi;
  double *ddrho_lo,*ddrho_hi;
  double lapse_interp_lo,lapse_interp_hi,lapse_interp;
  double shift_interp_lo,shift_interp_hi,shift_interp;
  double rho_interp_lo,rho_interp_hi,rho_interp;
  double ls_time_lo,ls_time_hi;
  double t_frac;
  double riso;
  double shift;
  double sinthcosph,sinthsinph,costh;
  double rad_s, rad_c;
  double junk; // Used to scan through the lapse and shift.
  int mark,i;
  double ONE=1.0,ZERO=0.0;
  long rpts;
    
  //HARDCODED FOR THE TIME BEING
  ifstream data_file;
  data_file.open("data_r10_new");
  data_file.read((char*)&rpts,sizeof(int));

  ls_rad_lo = dvector(1,rpts);
  lapse_lo = dvector(1,rpts);
  shift_lo = dvector(1,rpts);
  rho_lo = dvector(1,rpts);
  ddlapse_lo = dvector(1,rpts);
  ddshift_lo = dvector(1,rpts);
  ddrho_lo = dvector(1,rpts);
  ls_rad_hi = dvector(1,rpts);
  lapse_hi = dvector(1,rpts);
  shift_hi = dvector(1,rpts);
  rho_hi = dvector(1,rpts);
  ddlapse_hi = dvector(1,rpts);
  ddshift_hi = dvector(1,rpts);
  ddrho_hi = dvector(1,rpts);

 
  ls_time_hi=0.0;
  ls_time_lo=0.0;
  
 
 while(ls_time_hi <= cctk_time) {
    ls_time_lo = ls_time_hi;
    data_file.read((char*)&ls_time_hi,sizeof(double));
   
   
 
    for(i=1;i<=rpts;i++) {
      ls_rad_lo[i] = ls_rad_hi[i];
      lapse_lo[i] = lapse_hi[i];
      shift_lo[i] = shift_hi[i];
      rho_lo[i] = rho_hi[i];

      data_file.read((char*)&ls_rad_hi[i],sizeof(double));
      data_file.read((char*)&lapse_hi[i],sizeof(double));
      data_file.read((char*)&shift_hi[i],sizeof(double));
      data_file.read((char*)&rho_hi[i],sizeof(double));
      data_file.read((char*)&junk,sizeof(double));
      data_file.read((char*)&junk,sizeof(double));
      data_file.read((char*)&junk,sizeof(double));    
    }
    
    data_file.read((char*)&mark,sizeof(int));
    if(mark != -62171 || ls_time_lo == ls_time_hi) {
    
    }
    
      // spline
     
      spline(ls_rad_lo,lapse_lo,rpts,0.0,1.e31,ddlapse_lo);
      spline(ls_rad_hi,lapse_hi,rpts,0.0,1.e31,ddlapse_hi);
      spline(ls_rad_lo,shift_lo,rpts,0.0,1.e31,ddshift_lo);
      spline(ls_rad_hi,shift_hi,rpts,0.0,1.e31,ddshift_hi);
      spline(ls_rad_lo,rho_lo,rpts,0.0,1.e31,ddrho_lo);
      spline(ls_rad_hi,rho_hi,rpts,0.0,1.e31,ddrho_hi);
 }
    int istart = 0;
    int jstart = 0;
    int kstart = 0;
    int iend = cctk_lsh[0];
    int jend = cctk_lsh[1];
    int kend = cctk_lsh[2];

    //loop over all gridpoints
    for (int k=kstart; k<kend; k++) for (int j=jstart; j<jend; j++) for (int i=istart; i<iend; i++) {
      int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
      
      //physical isotropic radius
      riso=PhysicalRadius[vindex];
      
      //set radii (not physical radii if fisheye enabled)
      rad_s=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]+z[vindex]*z[vindex]);
      rad_c=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]);
      
      //set angles
      sinthcosph = x[vindex]/rad_s;
      sinthsinph = y[vindex]/rad_s;
      costh = z[vindex]/rad_s;
      
    //do interpolations
     
      splint(ls_rad_lo,lapse_lo,ddlapse_lo,rpts,riso,&lapse_interp_lo);
      splint(ls_rad_hi,lapse_hi,ddlapse_hi,rpts,riso,&lapse_interp_hi);
      splint(ls_rad_lo,shift_lo,ddshift_lo,rpts,riso,&shift_interp_lo);
      splint(ls_rad_hi,shift_hi,ddshift_hi,rpts,riso,&shift_interp_hi);
      splint(ls_rad_lo,rho_lo,ddrho_lo,rpts,riso,&rho_interp_lo);
      splint(ls_rad_hi,rho_hi,ddrho_hi,rpts,riso,&rho_interp_hi);

     
     

      t_frac = (cctk_time - ls_time_lo)/(ls_time_hi - ls_time_lo);
  
      lapm1[vindex] = t_frac*lapse_interp_hi + (1.-t_frac)*lapse_interp_lo - ONE;
      lapm1_p[vindex] = lapm1[vindex];
      lapm1_t[vindex] = lapm1[vindex];
      shift =  t_frac*shift_interp_hi + (1.-t_frac)*shift_interp_lo;
      shiftx[vindex] = shift*sinthcosph;
      shifty[vindex] = shift*sinthsinph;
      shiftz[vindex] = shift*costh;
      shiftx_p[vindex] = shiftx[vindex];
      shifty_p[vindex] = shifty[vindex];
      shiftz_p[vindex] = shiftz[vindex];
      shiftx_t[vindex] = shiftx[vindex];
      shifty_t[vindex] = shifty[vindex];
      shiftz_t[vindex] = shiftz[vindex];
      rho_b_analytic[vindex] = t_frac*rho_interp_hi + (1.-t_frac)*rho_interp_lo;
     }
}


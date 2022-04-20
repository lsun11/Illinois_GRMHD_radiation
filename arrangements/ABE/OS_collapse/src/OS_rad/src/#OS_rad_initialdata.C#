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
CCTK_FILEVERSION(OS_rad_initialdata)
  extern "C" void OS_rad_initialdata(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
 DECLARE_CCTK_PARAMETERS;
 //set EOS stuff
  double k_eos;
  double time; 
  int mark,i,rpts_mat,l,m;
  long rpts;
  //read in from file
  double *rad,*A,*r_areal,*rho_aux,*Sr,*Srr,*lapse,*shift;
  //derivs obtained from spline
  double *ddA,*ddrad,*ddrho_aux,*ddSr,*ddSrr,*ddlapse,*ddshift;
  //there are two radii needed, the spherical radius sqrt(x^2+y^2+z^2), and the cylindrical radius sqrt(x^2+y^2).  We denote these rad_s and rad_c respectively.  Note however, that if we are using fisheye, these are NOT physical radii
  double rad_s,rad_c;
  //we will set riso to be the physical radius, and the above radii will be useful only for finding angles, since the fisheye transformation preserves angles.
  double riso;
  //interpolated values
  double A_interp,r_edge_interp,lapse_interp,shift_interp,rho_aux_interp,Sr_interp,Srr_interp;
  //some angles
  double sinthcosph, sinthsinph, costh;
  //metric quantities in spherical coords
  double **g_pol;
  //transform from spherical to cartesian
  double **trans;
  //some useful constants
  double rounding_func,rho_b_0;
  double HALF=0.5,ONE=1.0,ZERO=0.0,TEN=10.0,f4o3=4.0/3.0,f3o4=3.0/4.0;
 
  trans = dmatrix(1,4,1,4);
  g_pol = dmatrix(1,4,1,4);
  
  printf("filename: %s\n",OS_filename);
 
  rho_b_0=ONE/(f4o3*M_PI*RoM*RoM*RoM); 

  ifstream data_file;
  data_file.open(OS_filename);
  data_file.read((char*)&rpts,sizeof(int));

  //Set up vectors
  rad = dvector(1,rpts);
  A = dvector(1,rpts);
  r_areal = dvector(1,rpts);
  rho_aux = dvector(1,rpts);
  Sr = dvector(1,rpts);
  Srr = dvector(1,rpts);
  lapse = dvector(1,rpts);
  shift = dvector(1,rpts);
  ddA = dvector(1,rpts);
  ddrad = dvector(1,rpts);
  ddrho_aux = dvector(1,rpts);
  ddSr = dvector(1,rpts);
  ddSrr = dvector(1,rpts);
  ddlapse = dvector(1,rpts);
  ddshift = dvector(1,rpts);
  
  
  
  //read in initial data
  data_file.read((char*)&time,sizeof(double));
  for(i=1;i<=rpts;i++) {
    data_file.read((char*)&rad[i],sizeof(double));
    data_file.read((char*)&lapse[i],sizeof(double));
    data_file.read((char*)&shift[i],sizeof(double));
    data_file.read((char*)&rho_aux[i],sizeof(double));
    data_file.read((char*)&Sr[i],sizeof(double));
    data_file.read((char*)&Srr[i],sizeof(double));
    data_file.read((char*)&A[i],sizeof(double));  
    
 
                          
  }
  for(i=1;i<=rpts;i++) {
    r_areal[i] = A[i]*rad[i];
  }
  printf("test1\n");

  data_file.read((char*)&mark,sizeof(int));                                      
if(mark != -62171)
    printf("DANGER! BAD MARK: %d\n",mark);                                              
  //spline
  spline(rad,A,rpts,0.0,1.e31,ddA);
  spline(rad,rad,rpts,0.0,1.e31,ddrad);
  spline(rad,lapse,rpts,0.0,1.e31,ddlapse);
  spline(rad,shift,rpts,0.0,1.e31,ddshift);
  //
  // Now that we've read the matter in, we want to ditch the exterior:
  // it's zero and splining through the sudden drop-off causes severe
  // problems.
  //          
  rpts_mat = -1000;
  for(i=1;i<=rpts;i++) {
    ddrho_aux[i] = 0.0;
    ddSr[i] = 0.0;
    ddSrr[i] = 0.0;
    if (rho_aux[i] == 0.0 && rpts_mat == -1000)
      rpts_mat = i - 1;
  }

  // spline
  if(rpts_mat > 0) {
    spline(rad,rho_aux,rpts_mat,0.0,1.e31,ddrho_aux);
    spline(rad,Sr,rpts_mat,0.0,1.e31,ddSr);
    spline(rad,Srr,rpts_mat,0.0,1.e31,ddSrr);
  }

  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];
  
  splint(r_areal,rad,ddrad,rpts,RoM,&r_edge_interp);
  printf("r_edge_interp: %10.6f\n",r_edge_interp);
  
  //particle tracer stuff
  for (int index=0;index<=narr-1;index++) {
    coord[index + 0*narr] = ZERO;
    splint(r_areal,rad,ddrad,rpts,(index+1)*(0.97*RoM/narr),&coord[index + 1*narr]); 
    coord[index + 2*narr] = ZERO;
    coord[index + 3*narr] = ZERO;
  } 
 
  //loop over all gridpoints
  for (int k=kstart; k<kend; k++) for (int j=jstart; j<jend; j++) for (int i=istart; i<iend; i++) {
    
    int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    //set radii (not physical radii if fisheye enabled)
    rad_s=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]+z[vindex]*z[vindex]);
    rad_c=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]);
    
    //physical isotropic radius
    riso=PhysicalRadius[vindex];
    
    //do interpolations
    splint(rad,A,ddA,rpts,riso,&A_interp);
   
    if (k==kstart && j==jstart && i==istart){
   
    }
    splint(rad,lapse,ddlapse,rpts,riso,&lapse_interp);
    splint(rad,shift,ddshift,rpts,riso,&shift_interp);
    if (riso > rad[rpts_mat]){
      rho_aux_interp = 0.0;
      Sr_interp = 0.0;
      Srr_interp = 0.0;
    }
    else {
      splint(rad,rho_aux,ddrho_aux,rpts_mat,riso,&rho_aux_interp);
      splint(rad,Sr,ddSr,rpts_mat,riso,&Sr_interp);
      splint(rad,Srr,ddSrr,rpts_mat,riso,&Srr_interp);
    }
    
    //    printf("%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",riso,lapse_interp,shift_interp,rho_aux_interp,Sr_interp,Srr_interp,drdt_interp,A_interp);         
    //I set the following because I believe that in the case of OS_rad collapse, rho=rho_b=rho_star at t=0.  But I think I should check this more carefully.

    //rho_b[vindex]=rho_aux_interp;
  
    //  if (rounding_temp > 0.0) {
      rounding_func = ONE/(exp((riso-r_edge_interp)/rounding_temp)+ONE); 
      //} else{
      // rounding_func = ONE;
      // }
    
    rho_b[vindex] = rounding_func*rho_b_0;

    
    //I'm not sure what I'm doing here
    // P[vindex]= .001;
    //vy[vindex]=0.0;
    
    //set angles
    sinthcosph = x[vindex]/rad_s;
    sinthsinph = y[vindex]/rad_s;
    costh = z[vindex]/rad_s;
   
    //lets try setting them to 0
    vx[vindex] = ZERO;
    vy[vindex] = ZERO;
    vz[vindex] = ZERO;
    
    //find various quantities in cartesian coords.  I'm not sure if it actually does anything to set Si and Sij here.  I think they get overwritten later anyway.
    /*
      Sx[vindex] = Sr_interp * sinthcosph;
      Sy[vindex] = Sr_interp * sinthsinph;
      Sz[vindex] = Sr_interp * costh;
      Sxx[vindex] = Srr_interp * sinthcosph * sinthcosph;
      Sxy[vindex] = Srr_interp * sinthcosph * sinthsinph;
      Sxz[vindex] = Srr_interp * sinthcosph * costh;
      Syy[vindex] = Srr_interp * sinthsinph * sinthsinph;
      Syz[vindex] = Srr_interp * sinthsinph * costh;
      Szz[vindex] = Srr_interp * costh * costh;
    */
    
    /*
    shiftx[vindex] = shift_interp * sinthcosph;
    shifty[vindex] = shift_interp * sinthsinph;
    shiftz[vindex] = shift_interp * costh;
    */
    
    shiftx[vindex] = ZERO;
    shifty[vindex] = ZERO;
    shiftz[vindex] = ZERO;
    
    if(CCTK_Equals(slicing_type,"geodesic") && (os_maximal != 1)){
      lapm1[vindex]=ZERO;
    } else {
      lapm1[vindex] = lapse_interp-ONE;
    } 
      //NOTE: this is temporary.  Only for geodesic slicing
    //lapm1[vindex]= ZERO;
       
  
    phi[vindex] = HALF * log(A_interp);

    //this is the tilde metric.  Since the spacetime is conformally flat, it is very simple
    gxx[vindex] = ONE;
    gxy[vindex] = ZERO;                                       
    gxz[vindex] = ZERO;                                                 
    gyy[vindex] = ONE;  
    gyz[vindex] = ZERO;                                                  
    gzz[vindex] = ONE;

    gupxx[vindex] = ONE;
    gupxy[vindex] = ZERO;                                       
    gupxz[vindex] = ZERO;                                                 
    gupyy[vindex] = ONE;  
    gupyz[vindex] = ZERO;                                                  
    gupzz[vindex] = ONE;

    Gammax[vindex] = ZERO;
    Gammay[vindex] = ZERO;
    Gammaz[vindex] = ZERO;

    trK[vindex] = ZERO;
     
    Axx[vindex] = ZERO;
    Axy[vindex] = ZERO;                                       
    Axz[vindex] = ZERO;                                                 
    Ayy[vindex] = ZERO;  
    Ayz[vindex] = ZERO;                                                  
    Azz[vindex] = ZERO;

    
  }

 
  int vindex = CCTK_GFINDEX3D(cctkGH,1,1,1);
  //n
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,0)] = n_index;
  //K_poly
  printf("rho_b_0: %10.6f\n)",rho_b_0);
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,1)] = PoRho*pow(rho_b_0,-ONE/n_index) ;
  //kappaa
  //set it to be opt_depth_a/(rho_b*Radius)  here I assume Mass=1
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,2)] = opt_depth_a/(rho_b_0* RoM);
  //kappas
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,3)] = opt_depth_s/(rho_b_0* RoM);

  //a_R*m_B^4
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,4)] = f3o4/(M_PI*Po4PiB*PoRho*PoRho*PoRho*RoM*RoM*RoM);

  
  //r_edge_interp
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,5)] = r_edge_interp;
  
  //rho_b_0
  Pr[CCTK_GFINDEX3D(cctkGH,0,0,6)] = rho_b_0;

  k_eos= Pr[CCTK_GFINDEX3D(cctkGH,0,0,1)];

  //do EOS stuff
  
  //it seems that I need to set n.  This should be in the .par file, but I'll put it here for now:
  cout << "Can't open " << endl;
  cout << "Assuming polytropic EOS..." << endl;
  cout << "Polytropic index = " << n_index << endl;
  cout << "Polytropic constant K = " << k_eos << endl;
  //rho_tab = new double[1]; P_tab = new double[1]; eps_tab = new double[1];
  //gamma_tab = new double[2]; k_tab = new double[2];
  rho_tab[0] = 0.03;
  P_tab[0] = k_eos * pow(rho_tab[0],1.0+1.0/ n_index);
  eps_tab[0] = n_index*P_tab[0]/rho_tab[0];
  gamma_tab[0] = 1.0 + 1.0/ n_index; k_tab[0] = k_eos;
  gamma_tab[1] = 1.0 + 1.0/ n_index; k_tab[1] = k_eos;

  data_file.close();
	 
 

}





#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sstream>
#include <time.h>
#include <ctype.h>
#include <iostream>
#include <sstream> 
#include <string>
#include <unistd.h>
#include <cmath>
#include <iomanip> 
#include <fstream> 
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
//#include "utilities.h"
using namespace std;
void interpolate_general(double *xco, double *yco, double *zco,double &x1output, double &x2output, double &x3output,double *var,int &nx,int &ny,int *i1, int *j1, int *k1,double &val,
		 double &f11,double &f12,double &f13,double &f21,double &f22, double &f23,double &f31,double &f32,double &f33,int compute_interp_coeffs);

extern "C" void CCTK_FCALL CCTK_FNAME(set_metric_rot_lowlevel)
  (const cGH **cctkGH,double *phi_rot_ang,double *xbh1_initial,double *xbh2_initial,
   double *radmin1,double *radmax1,double *radmin2,double *radmax2,double *radmin3,double *radmax3,
   int *Nradial,int *Ntheta,int *Nphi,
   double *K_rr_rot1,double *K_rth_rot1,double *K_rp_rot1,double *K_thth_rot1,double *K_thp_rot1, double *K_pp_rot1, 
   double *shiftr_rot1,double *shiftth_rot1,double *shiftp_rot1,double *phi_rot1,double *lapm1_rot1,
   double *K_rr_rot2,double *K_rth_rot2,double *K_rp_rot2,double *K_thth_rot2,double *K_thp_rot2, double *K_pp_rot2, 
   double *shiftr_rot2,double *shiftth_rot2,double *shiftp_rot2,double *phi_rot2,double *lapm1_rot2,
   double *K_rr_rot3,double *K_rth_rot3,double *K_rp_rot3,double *K_thth_rot3,double *K_thp_rot3, double *K_pp_rot3, 
   double *shiftr_rot3,double *shiftth_rot3,double *shiftp_rot3,double *phi_rot3,double *lapm1_rot3,
   double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
   double *shiftx,double *shifty,double *shiftz,double *phi,double *lapm1,
   double *X,double *Y,double *Z,
   double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext);

extern "C" void set_metric_rot_lowlevel(const cGH *cctkGH,double phi_rot_ang,double xbh1_initial,double xbh2_initial,
					double radmin1,double radmax1,double radmin2,double radmax2,double radmin3,double radmax3,
					int Nradial,int Ntheta,int Nphi,
					double *K_rr_rot1,double *K_rth_rot1,double *K_rp_rot1,double *K_thth_rot1,double *K_thp_rot1,double *K_pp_rot1,
					double *shiftr_rot1,double *shiftth_rot1,double *shiftp_rot1,double *phi_rot1,double *lapm1_rot1,
					double *K_rr_rot2,double *K_rth_rot2,double *K_rp_rot2,double *K_thth_rot2,double *K_thp_rot2,double *K_pp_rot2,
					double *shiftr_rot2,double *shiftth_rot2,double *shiftp_rot2,double *phi_rot2,double *lapm1_rot2,
					double *K_rr_rot3,double *K_rth_rot3,double *K_rp_rot3,double *K_thth_rot3,double *K_thp_rot3,double *K_pp_rot3,
					double *shiftr_rot3,double *shiftth_rot3,double *shiftp_rot3,double *phi_rot3,double *lapm1_rot3,
					double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
					double *shiftx,double *shifty,double *shiftz,double *phi,double *lapm1,
					double *X,double *Y,double *Z,
					double xmin,double ymin,double zmin,
					double dx,double dy,double dz,int *ext) {
  
  int nx = ext[0];
  int ny = ext[1];
  int nz = ext[2];
  
  double PI=acos(-1.0);

  double thetamin=-PI/1000.;
  double thetamax=PI + PI/1000.;
  double phimin=-PI/1000.;
  double phimax=2.0*PI+PI/1000.;
  
  double logradmin1=log(radmin1);
  double logradmax1=log(radmax1);
  double logradmin2=log(radmin2);
  double logradmax2=log(radmax2);
  double logradmin3=log(radmin3);
  double logradmax3=log(radmax3);
  double radial1_arr[Nradial];
  double radial2_arr[Nradial];
  double radial3_arr[Nradial];
  double theta_arr[Ntheta];
  double phi_arr[Nphi];

  cout << "beginning of set_metric_rot_lowlevel" << endl;

  //  int logarithmic;
  
  //  for (int i=0;i<Nradial*Ntheta*Nphi;i=i+Ntheta*Nphi){
  //   cout << "lapm1_rot3: " << lapm1_rot3[i] << endl;
  // }
  
  //don't hardcode this
#pragma omp parallel for
  for (int i=0;i<Nradial;i++){
    //if (logarithmic==1){
    //radial1_arr[i]=logradmin1+(logradmax1-logradmin1)/(Nradial-1.0)*i;
    //radial2_arr[i]=logradmin2+(logradmax2-logradmin2)/(Nradial-1.0)*i;
      radial3_arr[i]=logradmin3+(logradmax3-logradmin3)/(Nradial-1.0)*i;
      //} else{
      radial1_arr[i]=radmin1+(radmax1-radmin1)/(Nradial-1.0)*i;
      radial2_arr[i]=radmin2+(radmax2-radmin2)/(Nradial-1.0)*i; 
      //radial3_arr[i]=radmin3+(radmax3-radmin3)/(Nradial-1.0)*i;
      //}
      //      cout << "radial3_arr[i]: " << radial3_arr[i] << endl;
  }
#pragma omp parallel for
  for (int i=0;i<Ntheta;i++){
    theta_arr[i]=thetamin+(thetamax-thetamin)/(Ntheta-1.0)*i;
  } 
#pragma omp parallel for
  for (int i=0;i<Nphi;i++){
    phi_arr[i]=phimin+(phimax-phimin)/(Nphi-1.0)*i;
  }
#pragma omp parallel for
  for (int k = 0; k < nz; k++) for (int j = 0; j < ny; j++) for (int i = 0; i < nx; i++) {
    int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    double radL=sqrt(X[index]*X[index]+Y[index]*Y[index]+Z[index]*Z[index])+1e-10;
    double rad_cylL=sqrt(X[index]*X[index]+Y[index]*Y[index])+1e-10;
    
    double Lrx = X[index]/radL;
    double Lry = Y[index]/radL;
    double Lrz = Z[index]/radL;
    double Lpx = -Y[index]/rad_cylL/rad_cylL;
    double Lpy = X[index]/rad_cylL/rad_cylL;
    double Lpz = 0;
    double Lthx = X[index]*Z[index]/(radL*radL*rad_cylL);
    double Lthy = Y[index]*Z[index]/(radL*radL*rad_cylL); 
    double Lthz = -rad_cylL/radL/radL;
   
    double Lxp=-Y[index];
    double Lyp=X[index];
    double Lzp=0.0;
    
    double Lxth=X[index]*Z[index]/rad_cylL;
    double Lyth=Y[index]*Z[index]/rad_cylL;
    double Lzth=-rad_cylL;
    double Lxr=X[index]/radL;
    double Lyr=Y[index]/radL;
    double Lzr=Z[index]/radL;
	
    double x1L=X[index]-(xbh1_initial)*cos(phi_rot_ang);
    double x2L=X[index]-(xbh2_initial)*cos(phi_rot_ang);
    double x3L=X[index];
    double y1L=Y[index]-(xbh1_initial)*sin(phi_rot_ang);
    double y2L=Y[index]-(xbh2_initial)*sin(phi_rot_ang);
    double y3L=Y[index];
    double zL=Z[index];    
    double rad1L=sqrt(x1L*x1L+y1L*y1L+zL*zL)+1e-10;
    double rad2L=sqrt(x2L*x2L+y2L*y2L+zL*zL)+1e-10;
    double rad3L=sqrt(x3L*x3L+y3L*y3L+zL*zL)+1e-10;
    double logradL;
    double radial_coordL;
    double thetaL;
    double phiL;
    int radial_index;
    int theta_index;
    int phi_index;
    double K_rr_rotL,K_rth_rotL,K_rp_rotL,K_thth_rotL,K_thp_rotL,K_pp_rotL,lapm1_rotL,shiftr_rotL,shiftth_rotL,shiftp_rotL,phi_rotL;
    
    //it is a good idea to use the origin centered grid for points near the origin.  You get less problems due to interpolation errors.  We first set points that are a distance less than origin_zone_rad from origin.
    double origin_zone_rad=1.0;
    if (rad3L<=origin_zone_rad){
      // logarithmic=1;
      if (rad3L<=radmin3){
	K_rr_rotL=K_rr_rot3[0];
	K_rth_rotL=K_rth_rot3[0];
	K_rp_rotL=K_rp_rot3[0];
	K_thth_rotL=K_thth_rot3[0];
	K_thp_rotL=K_thp_rot3[0];
	K_pp_rotL=K_pp_rot3[0];
	lapm1_rotL=lapm1_rot3[0];
	shiftr_rotL=shiftr_rot3[0];
	shiftth_rotL=shiftth_rot3[0];
	shiftp_rotL=shiftp_rot3[0];
	phi_rotL=phi_rot3[0];
      } else{
	radL=sqrt(x3L*x3L+y3L*y3L+zL*zL)+1e-10;
	rad_cylL=sqrt(x3L*x3L+y3L*y3L)+1e-10;
	logradL=log(radL);
	thetaL=acos(zL/radL);
	phiL=fmod(acos(x3L/rad_cylL)*((y3L > 0.0)-0.5)*2.0+2.0*PI*(y3L<=0.0)-phi_rot_ang+1000.0*2.0*PI,2.0*PI);   
	
	//	if (logarithmic==1){
	  radial_index=(logradL-logradmin3)/(logradmax3-logradmin3)*(Nradial-1);
	  radial_coordL=logradL;
	  //	} else{
	  //	  radial_index=(radL-radmin3)/(radmax3-radmin3)*(Nradial-1);
	  //	  radial_coordL=radL;
	  //	}
	
	theta_index=(thetaL-thetamin)/(thetamax-thetamin)*(Ntheta-1);
	phi_index=(phiL-phimin)/(phimax-phimin)*(Nphi-1);
	
	int i1[3],j1[3],k1[3];
	if(theta_index==0) theta_index++;
	i1[0] = theta_index - 1; i1[1] = theta_index; i1[2] = theta_index + 1;
	if(phi_index==0) phi_index++;
	j1[0] = phi_index - 1; j1[1] = phi_index; j1[2] = phi_index + 1;
	if(radial_index==0) radial_index++;
	k1[0] = radial_index - 1; k1[1] = radial_index; k1[2] = radial_index + 1;
	double f11,f12,f13,f21,f22,f23,f31,f32,f33;
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_rr_rot3,Ntheta,Nphi,i1,j1,k1,K_rr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,1);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_rth_rot3,Ntheta,Nphi,i1,j1,k1,K_rth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_rp_rot3,Ntheta,Nphi,i1,j1,k1,K_rp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_thth_rot3,Ntheta,Nphi,i1,j1,k1,K_thth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_thp_rot3,Ntheta,Nphi,i1,j1,k1,K_thp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_pp_rot3,Ntheta,Nphi,i1,j1,k1,K_pp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,lapm1_rot3,Ntheta,Nphi,i1,j1,k1,lapm1_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,shiftr_rot3,Ntheta,Nphi,i1,j1,k1,shiftr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,shiftth_rot3,Ntheta,Nphi,i1,j1,k1,shiftth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,shiftp_rot3,Ntheta,Nphi,i1,j1,k1,shiftp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,phi_rot3,Ntheta,Nphi,i1,j1,k1,phi_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
      }
      kxx[index] = Lrx*(Lrx*K_rr_rotL + Lthx*K_rth_rotL + Lpx*K_rp_rotL)
	+ Lthx*(Lrx*K_rth_rotL + Lthx*K_thth_rotL + Lpx*K_thp_rotL)
	+ Lpx*(Lrx*K_rp_rotL + Lthx*K_thp_rotL + Lpx*K_pp_rotL);
      
      kxy[index] = Lrx*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL)
	+ Lthx*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL)
	+ Lpx*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
      
      kxz[index] = Lrx*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthx*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpx*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
      
      kyy[index] = Lry*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL) 
	+ Lthy*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL) 
	+ Lpy*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
      
      kyz[index] = Lry*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthy*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpy*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
      
      kzz[index] = Lrz*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthz*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpz*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
      
      lapm1[index] = lapm1_rotL;
      
      shiftx[index] = Lxr*shiftr_rotL+Lxth*shiftth_rotL+Lxp*shiftp_rotL;
      shifty[index] = Lyr*shiftr_rotL+Lyth*shiftth_rotL+Lyp*shiftp_rotL;
      shiftz[index] = Lzr*shiftr_rotL+Lzth*shiftth_rotL+Lzp*shiftp_rotL;
      
      phi[index] = phi_rotL;
    } else if (rad1L <=exp(logradmax1)){ 
      // logarithmic=0;
      if (rad1L<=radmin1){
	K_rr_rotL=K_rr_rot1[0];
	K_rth_rotL=K_rth_rot1[0];
	K_rp_rotL=K_rp_rot1[0];
	K_thth_rotL=K_thth_rot1[0];
	K_thp_rotL=K_thp_rot1[0];
	K_pp_rotL=K_pp_rot1[0];
	lapm1_rotL=lapm1_rot1[0];
	shiftr_rotL=shiftr_rot1[0];
	shiftth_rotL=shiftth_rot1[0];
	shiftp_rotL=shiftp_rot1[0];
	phi_rotL=phi_rot1[0];
      } else{
	radL=sqrt(x1L*x1L+y1L*y1L+zL*zL)+1e-10;
	rad_cylL=sqrt(x1L*x1L+y1L*y1L)+1e-10;
	logradL=log(radL);
	thetaL=acos(zL/radL);
	phiL=fmod(acos(x1L/rad_cylL)*((y1L > 0.0)-0.5)*2.0+2.0*PI*(y1L<=0.0)-phi_rot_ang+1000.0*2.0*PI,2.0*PI);   
	    
	//	if (logarithmic==1){
	//          radial_index=(logradL-logradmin1)/(logradmax1-logradmin1)*(Nradial-1);
	//	  radial_coordL=logradL;
	//	} else{
          radial_index=(radL-radmin1)/(radmax1-radmin1)*(Nradial-1); 
	  radial_coordL=radL;
	  //	}  
	
	theta_index=(thetaL-thetamin)/(thetamax-thetamin)*(Ntheta-1);
	phi_index=(phiL-phimin)/(phimax-phimin)*(Nphi-1);
		
	int i1[3],j1[3],k1[3];
	if(theta_index==0) theta_index++;
	i1[0] = theta_index - 1; i1[1] = theta_index; i1[2] = theta_index + 1;
	if(phi_index==0) phi_index++;
	j1[0] = phi_index - 1; j1[1] = phi_index; j1[2] = phi_index + 1;
	if(radial_index==0) radial_index++;
	k1[0] = radial_index - 1; k1[1] = radial_index; k1[2] = radial_index + 1;
	double f11,f12,f13,f21,f22,f23,f31,f32,f33;
	
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,K_rr_rot1,Ntheta,Nphi,i1,j1,k1,K_rr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,1);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,K_rth_rot1,Ntheta,Nphi,i1,j1,k1,K_rth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,K_rp_rot1,Ntheta,Nphi,i1,j1,k1,K_rp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,K_thth_rot1,Ntheta,Nphi,i1,j1,k1,K_thth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,K_thp_rot1,Ntheta,Nphi,i1,j1,k1,K_thp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,K_pp_rot1,Ntheta,Nphi,i1,j1,k1,K_pp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,lapm1_rot1,Ntheta,Nphi,i1,j1,k1,lapm1_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,shiftr_rot1,Ntheta,Nphi,i1,j1,k1,shiftr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,shiftth_rot1,Ntheta,Nphi,i1,j1,k1,shiftth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,shiftp_rot1,Ntheta,Nphi,i1,j1,k1,shiftp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial1_arr,thetaL,phiL,radial_coordL,phi_rot1,Ntheta,Nphi,i1,j1,k1,phi_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
      } 
      kxx[index] = Lrx*(Lrx*K_rr_rotL + Lthx*K_rth_rotL + Lpx*K_rp_rotL)
	+ Lthx*(Lrx*K_rth_rotL + Lthx*K_thth_rotL + Lpx*K_thp_rotL)
	+ Lpx*(Lrx*K_rp_rotL + Lthx*K_thp_rotL + Lpx*K_pp_rotL);
	  
      kxy[index] = Lrx*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL)
	+ Lthx*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL)
	+ Lpx*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
	  
      kxz[index] = Lrx*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthx*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpx*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
    
      kyy[index] = Lry*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL) 
	+ Lthy*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL)   
	+ Lpy*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
    
      kyz[index] = Lry*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthy*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpy*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
    
      kzz[index] = Lrz*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthz*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpz*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
      
      lapm1[index] = lapm1_rotL;
    
      shiftx[index] = Lxr*shiftr_rotL+Lxth*shiftth_rotL+Lxp*shiftp_rotL;
      shifty[index] = Lyr*shiftr_rotL+Lyth*shiftth_rotL+Lyp*shiftp_rotL;
      shiftz[index] = Lzr*shiftr_rotL+Lzth*shiftth_rotL+Lzp*shiftp_rotL;
    
      phi[index] = phi_rotL;
      //} else if (rad2L<=fabs(xbh2_initial)*0.5+0.0001){
    } else if (rad2L<=exp(logradmax2)){
      // logarithmic=0;
      if (rad2L<=radmin2){
	K_rr_rotL=K_rr_rot2[0];
	K_rth_rotL=K_rth_rot2[0];
	K_rp_rotL=K_rp_rot2[0];
	K_thth_rotL=K_thth_rot2[0];
	K_thp_rotL=K_thp_rot2[0];
	K_pp_rotL=K_pp_rot2[0];
	lapm1_rotL=lapm1_rot2[0];
	shiftr_rotL=shiftr_rot2[0];
	shiftth_rotL=shiftth_rot2[0];
	shiftp_rotL=shiftp_rot2[0];
	phi_rotL=phi_rot2[0];
      } else{
	radL=sqrt(x2L*x2L+y2L*y2L+zL*zL)+1e-10;
	rad_cylL=sqrt(x2L*x2L+y2L*y2L)+1e-10;
	logradL=log(radL);
	thetaL=acos(zL/radL);
	phiL=fmod(acos(x2L/rad_cylL)*((y2L > 0.0)-0.5)*2.0+2.0*PI*(y2L<=0.0)-phi_rot_ang+1000.0*2.0*PI,2.0*PI);   
	
	//	if (logarithmic==1){
	//          radial_index=(logradL-logradmin2)/(logradmax2-logradmin2)*(Nradial-1); 
	//	  radial_coordL=logradL;
	//	} else{
          radial_index=(radL-radmin2)/(radmax2-radmin2)*(Nradial-1); 
	  radial_coordL=radL;
	  //	}
	theta_index=(thetaL-thetamin)/(thetamax-thetamin)*(Ntheta-1);
	phi_index=(phiL-phimin)/(phimax-phimin)*(Nphi-1);

	int i1[3],j1[3],k1[3];
	if(theta_index==0) theta_index++;
	i1[0] = theta_index - 1; i1[1] = theta_index; i1[2] = theta_index + 1;
	if(phi_index==0) phi_index++;
	j1[0] = phi_index - 1; j1[1] = phi_index; j1[2] = phi_index + 1;
	if(radial_index==0) radial_index++;
	k1[0] = radial_index - 1; k1[1] = radial_index; k1[2] = radial_index + 1;
	double f11,f12,f13,f21,f22,f23,f31,f32,f33;
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,K_rr_rot2,Ntheta,Nphi,i1,j1,k1,K_rr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,1);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,K_rth_rot2,Ntheta,Nphi,i1,j1,k1,K_rth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,K_rp_rot2,Ntheta,Nphi,i1,j1,k1,K_rp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,K_thth_rot2,Ntheta,Nphi,i1,j1,k1,K_thth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,K_thp_rot2,Ntheta,Nphi,i1,j1,k1,K_thp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,K_pp_rot2,Ntheta,Nphi,i1,j1,k1,K_pp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,lapm1_rot2,Ntheta,Nphi,i1,j1,k1,lapm1_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,shiftr_rot2,Ntheta,Nphi,i1,j1,k1,shiftr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,shiftth_rot2,Ntheta,Nphi,i1,j1,k1,shiftth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,shiftp_rot2,Ntheta,Nphi,i1,j1,k1,shiftp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial2_arr,thetaL,phiL,radial_coordL,phi_rot2,Ntheta,Nphi,i1,j1,k1,phi_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
      }
      kxx[index] = Lrx*(Lrx*K_rr_rotL + Lthx*K_rth_rotL + Lpx*K_rp_rotL)
	+ Lthx*(Lrx*K_rth_rotL + Lthx*K_thth_rotL + Lpx*K_thp_rotL)
	+ Lpx*(Lrx*K_rp_rotL + Lthx*K_thp_rotL + Lpx*K_pp_rotL);
    
      kxy[index] = Lrx*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL)
	+ Lthx*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL)
	+ Lpx*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
    
      kxz[index] = Lrx*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthx*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpx*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
    
      kyy[index] = Lry*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL) 
	+ Lthy*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL)  
	+ Lpy*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
    
      kyz[index] = Lry*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthy*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpy*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
    
      kzz[index] = Lrz*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthz*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpz*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
      
      lapm1[index] = lapm1_rotL;
    
      shiftx[index] = Lxr*shiftr_rotL+Lxth*shiftth_rotL+Lxp*shiftp_rotL;
      shifty[index] = Lyr*shiftr_rotL+Lyth*shiftth_rotL+Lyp*shiftp_rotL;
      shiftz[index] = Lzr*shiftr_rotL+Lzth*shiftth_rotL+Lzp*shiftp_rotL;
    
      phi[index] = phi_rotL;
      
    } else {
      // logarithmic=1;
      if (rad3L<=radmin3){
	K_rr_rotL=K_rr_rot3[0];
	K_rth_rotL=K_rth_rot3[0];
	K_rp_rotL=K_rp_rot3[0];
	K_thth_rotL=K_thth_rot3[0];
	K_thp_rotL=K_thp_rot3[0];
	K_pp_rotL=K_pp_rot3[0];
	lapm1_rotL=lapm1_rot3[0];
	shiftr_rotL=shiftr_rot3[0];
	shiftth_rotL=shiftth_rot3[0];
	shiftp_rotL=shiftp_rot3[0];
	phi_rotL=phi_rot3[0];
      } else{
	radL=sqrt(x3L*x3L+y3L*y3L+zL*zL)+1e-10;
	rad_cylL=sqrt(x3L*x3L+y3L*y3L)+1e-10;
	logradL=log(radL);
	thetaL=acos(zL/radL);
	phiL=fmod(acos(x3L/rad_cylL)*((y3L > 0.0)-0.5)*2.0+2.0*PI*(y3L<=0.0)-phi_rot_ang+1000.0*2.0*PI,2.0*PI);   
	
	//	if (logarithmic==1){ 
	  radial_index=(logradL-logradmin3)/(logradmax3-logradmin3)*(Nradial-1); 
	  radial_coordL=logradL;
	  //	} else{
	  //          radial_index=(radL-radmin3)/(radmax3-radmin3)*(Nradial-1); 
	  //	  radial_coordL=radL;
	  //	}   
	theta_index=(thetaL-thetamin)/(thetamax-thetamin)*(Ntheta-1);
	phi_index=(phiL-phimin)/(phimax-phimin)*(Nphi-1);
	
	int i1[3],j1[3],k1[3];
	if(theta_index==0) theta_index++;
	i1[0] = theta_index - 1; i1[1] = theta_index; i1[2] = theta_index + 1;
	if(phi_index==0) phi_index++;
	j1[0] = phi_index - 1; j1[1] = phi_index; j1[2] = phi_index + 1;
	if(radial_index==0) radial_index++;
	k1[0] = radial_index - 1; k1[1] = radial_index; k1[2] = radial_index + 1;
	double f11,f12,f13,f21,f22,f23,f31,f32,f33;
	
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_rr_rot3,Ntheta,Nphi,i1,j1,k1,K_rr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,1);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_rth_rot3,Ntheta,Nphi,i1,j1,k1,K_rth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_rp_rot3,Ntheta,Nphi,i1,j1,k1,K_rp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_thth_rot3,Ntheta,Nphi,i1,j1,k1,K_thth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_thp_rot3,Ntheta,Nphi,i1,j1,k1,K_thp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,K_pp_rot3,Ntheta,Nphi,i1,j1,k1,K_pp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,lapm1_rot3,Ntheta,Nphi,i1,j1,k1,lapm1_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,shiftr_rot3,Ntheta,Nphi,i1,j1,k1,shiftr_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,shiftth_rot3,Ntheta,Nphi,i1,j1,k1,shiftth_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,shiftp_rot3,Ntheta,Nphi,i1,j1,k1,shiftp_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
	interpolate_general(theta_arr,phi_arr,radial3_arr,thetaL,phiL,radial_coordL,phi_rot3,Ntheta,Nphi,i1,j1,k1,phi_rotL,f11,f12,f13,f21,f22,f23,f31,f32,f33,0);
      }
      kxx[index] = Lrx*(Lrx*K_rr_rotL + Lthx*K_rth_rotL + Lpx*K_rp_rotL)
	+ Lthx*(Lrx*K_rth_rotL + Lthx*K_thth_rotL + Lpx*K_thp_rotL)
	+ Lpx*(Lrx*K_rp_rotL + Lthx*K_thp_rotL + Lpx*K_pp_rotL);
    
      kxy[index] = Lrx*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL)
	+ Lthx*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL)
	+ Lpx*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
    
      kxz[index] = Lrx*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthx*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpx*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
    
      kyy[index] = Lry*(Lry*K_rr_rotL + Lthy*K_rth_rotL + Lpy*K_rp_rotL) 
	+ Lthy*(Lry*K_rth_rotL + Lthy*K_thth_rotL + Lpy*K_thp_rotL) 
	+ Lpy*(Lry*K_rp_rotL + Lthy*K_thp_rotL + Lpy*K_pp_rotL);
    
      kyz[index] = Lry*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthy*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpy*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
    
      kzz[index] = Lrz*(Lrz*K_rr_rotL + Lthz*K_rth_rotL + Lpz*K_rp_rotL) 
	+ Lthz*(Lrz*K_rth_rotL + Lthz*K_thth_rotL + Lpz*K_thp_rotL) 
	+ Lpz*(Lrz*K_rp_rotL + Lthz*K_thp_rotL + Lpz*K_pp_rotL);
      
      lapm1[index] = lapm1_rotL;
    
      shiftx[index] = Lxr*shiftr_rotL+Lxth*shiftth_rotL+Lxp*shiftp_rotL;
      shifty[index] = Lyr*shiftr_rotL+Lyth*shiftth_rotL+Lyp*shiftp_rotL;
      shiftz[index] = Lzr*shiftr_rotL+Lzth*shiftth_rotL+Lzp*shiftp_rotL;
      
      phi[index] = phi_rotL;
      
      
    }
  }
  cout << "end of set_metric_rot_lowlevel" << endl;
}

extern "C" void CCTK_FCALL CCTK_FNAME(set_metric_rot_lowlevel)
  (const cGH **cctkGH,double *phi_rot_ang,double *xbh1_initial,double *xbh2_initial,
   double *radmin1,double *radmax1,double *radmin2,double *radmax2,double *radmin3,double *radmax3,
   int *Nradial,int *Ntheta,int *Nphi,
   double *K_rr_rot1,double *K_rth_rot1,double *K_rp_rot1,double *K_thth_rot1,double *K_thp_rot1,double *K_pp_rot1,
   double *shiftr_rot1,double *shiftth_rot1,double *shiftp_rot1,double *phi_rot1,double *lapm1_rot1,
   double *K_rr_rot2,double *K_rth_rot2,double *K_rp_rot2,double *K_thth_rot2,double *K_thp_rot2,double *K_pp_rot2,
   double *shiftr_rot2,double *shiftth_rot2,double *shiftp_rot2,double *phi_rot2,double *lapm1_rot2,
   double *K_rr_rot3,double *K_rth_rot3,double *K_rp_rot3,double *K_thth_rot3,double *K_thp_rot3,double *K_pp_rot3,
   double *shiftr_rot3,double *shiftth_rot3,double *shiftp_rot3,double *phi_rot3,double *lapm1_rot3,
   double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
   double *shiftx,double *shifty,double *shiftz,double *phi,double *lapm1,
   double *X,double *Y,double *Z,
   double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext)
{  
  set_metric_rot_lowlevel(*cctkGH,*phi_rot_ang,*xbh1_initial,*xbh2_initial,
			  *radmin1,*radmax1,*radmin2,*radmax2,*radmin3,*radmax3,
			  *Nradial,*Ntheta,*Nphi,
			  K_rr_rot1,K_rth_rot1,K_rp_rot1,K_thth_rot1,K_thp_rot1,K_pp_rot1,
			  shiftr_rot1,shiftth_rot1,shiftp_rot1,phi_rot1,lapm1_rot1,
			  K_rr_rot2,K_rth_rot2,K_rp_rot2,K_thth_rot2,K_thp_rot2,K_pp_rot2,
			  shiftr_rot2,shiftth_rot2,shiftp_rot2,phi_rot2,lapm1_rot2,
			  K_rr_rot3,K_rth_rot3,K_rp_rot3,K_thth_rot3,K_thp_rot3,K_pp_rot3,
			  shiftr_rot3,shiftth_rot3,shiftp_rot3,phi_rot3,lapm1_rot3,
			  kxx,kxy,kxz,kyy,kyz,kzz,shiftx,shifty,shiftz,phi,lapm1,
			  X,Y,Z,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext);
}

// -------------------------------------------------------------------------------

void interpolate_general(double *xco, double *yco, double *zco,
		 double &x1output, double &x2output, double &x3output,
		 double *var,int &nx,int &ny,int *i1, int *j1, int *k1,
		 double &val,
		 double &f11,double &f12,double &f13,double &f21,double &f22, double &f23,double &f31,double &f32,double &f33,
		 int compute_interp_coeffs)
{
  
  int nxny=nx*ny;
  double stencil[3][3][3];
  if(compute_interp_coeffs==1) {
  
    // Here we set coefficients for interpolation
    //
    // Interpolation coefficients in each dimension are set by solving the following 
    // equation for (a,b,c)[= e.g., (f11,f12,f13)]:
    // a f(x1) + b f(x2) + c f(x3) = f(x5), where (x1,x2,x3) are the coordinate points in
    //    the interpolation stencil, and x5 is the interpolated destination point.
    // If f(x) = k = constant, we have
    // a + b + c = 1
    // If f(x) = k x, we have
    // a x1 + b x2 + c x3 = x5
    // Similarly for f(x) = k x^2, we have:
    // a x1^2 + b x2^2 + c x3^2 = x5^2

    // Using Mathematica, we can solve these equations easily:
    //   CForm[FullSimplify[Solve[{a + b + c == 1, a*x1 + b*x2 + c*x3 == x5, a*x1^2 + b*x2^2 + c*x3^2 == x5^2}, {a, b, c}]]]
    // This has solution: 
    // a = ((x2 - x5)*( x3 - x5))/((x1 - x2)*(x1 - x3))
    // b = ((x1 - x5)*(-x3 + x5))/((x1 - x2)*(x2 - x3))
    // c = ((x1 - x5)*( x2 - x5))/((x1 - x3)*(x2 - x3))

    double x1a[3],x2a[3],x3a[3];
  
    x1a[0] = xco[i1[0]]; x1a[1] = xco[i1[1]]; x1a[2] = xco[i1[2]]; 
    x2a[0] = yco[j1[0]]; x2a[1] = yco[j1[1]]; x2a[2] = yco[j1[2]]; 
    x3a[0] = zco[k1[0]]; x3a[1] = zco[k1[1]]; x3a[2] = zco[k1[2]]; 

    double x1,x2,x3,x5;
    
    x1 = x1a[0];
    x2 = x1a[1];
    x3 = x1a[2];
    x5 = x1output;

    //fij indicates direction i, point number j in the interpolation stencil, as in j=1->a, j=2->b, j=3->c
    f11 = ((x2 - x5)*( x3 - x5))/((x1 - x2)*(x1 - x3));
    f12 = ((x1 - x5)*(-x3 + x5))/((x1 - x2)*(x2 - x3));
    f13 = ((x1 - x5)*( x2 - x5))/((x1 - x3)*(x2 - x3));
  
    x1 = x2a[0];
    x2 = x2a[1];
    x3 = x2a[2];
    x5 = x2output;
    
    f21 = ((x2 - x5)*( x3 - x5))/((x1 - x2)*(x1 - x3));
    f22 = ((x1 - x5)*(-x3 + x5))/((x1 - x2)*(x2 - x3));
    f23 = ((x1 - x5)*( x2 - x5))/((x1 - x3)*(x2 - x3));
    
    x1 = x3a[0];
    x2 = x3a[1];
    x3 = x3a[2];
    x5 = x3output;
    
    f31 = ((x2 - x5)*( x3 - x5))/((x1 - x2)*(x1 - x3));
    f32 = ((x1 - x5)*(-x3 + x5))/((x1 - x2)*(x2 - x3));
    f33 = ((x1 - x5)*( x2 - x5))/((x1 - x3)*(x2 - x3));
  }

  //Second-order Lagrange interpolation <-> 3x3x3 stencil
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) {   
    int ind = i1[i] + nx*j1[j] + nxny*k1[k]; 
    stencil[i][j][k] = var[ind];
  }

 
  val =  
    + f11*f21*f31*stencil[0][0][0]
    + f11*f21*f32*stencil[0][0][1]
    + f11*f21*f33*stencil[0][0][2]
    + f11*f22*f31*stencil[0][1][0]
    + f11*f22*f32*stencil[0][1][1]
    + f11*f22*f33*stencil[0][1][2]
    + f11*f23*f31*stencil[0][2][0]
    + f11*f23*f32*stencil[0][2][1]
    + f11*f23*f33*stencil[0][2][2]
    + f12*f21*f31*stencil[1][0][0]
    + f12*f21*f32*stencil[1][0][1]
    + f12*f21*f33*stencil[1][0][2]
    + f12*f22*f31*stencil[1][1][0]
    + f12*f22*f32*stencil[1][1][1]
    + f12*f22*f33*stencil[1][1][2]
    + f12*f23*f31*stencil[1][2][0]
    + f12*f23*f32*stencil[1][2][1]
    + f12*f23*f33*stencil[1][2][2]
    + f13*f21*f31*stencil[2][0][0]
    + f13*f21*f32*stencil[2][0][1]
    + f13*f21*f33*stencil[2][0][2]
    + f13*f22*f31*stencil[2][1][0]
    + f13*f22*f32*stencil[2][1][1]
    + f13*f22*f33*stencil[2][1][2]
    + f13*f23*f31*stencil[2][2][0]
    + f13*f23*f32*stencil[2][2][1]
    + f13*f23*f33*stencil[2][2][2];

}
//------------------------------------------------------------------------------------

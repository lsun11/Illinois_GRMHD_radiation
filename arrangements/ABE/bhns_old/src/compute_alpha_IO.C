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

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_compute_alpha_IO)
  (const cGH **cctkGH, double &time,
   double &varpi_min,double &dvarpi,double &dz,
   int &N_varpi,int &N_phi,int &N_Z,
   double *rho_b_int,double *P_int,double *h_int,
   double *vx_int,double *vy_int,double *vz_int,double *u0_int,
   //double *smallb2_int,
   double *sbt_int,double *sbx_int,double *sby_int,double *sbz_int,
   double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
   double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
   double *gyy_int,double *gyz_int, double *gzz_int,
   double *gupxx_int,double *gupxy_int,double *gupxz_int,
   double *gupyy_int,double *gupyz_int, double *gupzz_int, double &bh_posn_x, double &bh_posn_y, double &bh_posn_z, double &alpha_rho_cut_off);
//   double *sigma_xx_int,double *sigma_xy_int,double *sigma_xz_int,
//   double *sigma_yy_int,double *sigma_yz_int, double *sigma_zz_int,
//   double *eta_vis_int);

extern "C" void bhns_compute_alpha_IO(const cGH *cctkGH, double &time,
				      double &varpi_min,double &dvarpi,double &dz,
				      int &N_varpi,int &N_phi,int &N_Z,
				      double *rho_b_int,double *P_int,double *h_int,
				      double *vx_int,double *vy_int,double *vz_int,double *u0_int,
				      //double *smallb2_int,
				      double *sbt_int,double *sbx_int,double *sby_int,double *sbz_int,
				      double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
				      double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
				      double *gyy_int,double *gyz_int, double *gzz_int,
				      double *gupxx_int,double *gupxy_int,double *gupxz_int,
				      double *gupyy_int,double *gupyz_int, double *gupzz_int,
				      double &bh_posn_x, double &bh_posn_y, double &bh_posn_z,double &alpha_rho_cut_off){
  //				 double *sigma_xx_int,double *sigma_xy_int,double *sigma_xz_int,
  //				 double *sigma_yy_int,double *sigma_yz_int, double *sigma_zz_int,
  //				 double *eta_vis_int){
	

  printf("BH position %e %e %e\n",bh_posn_x,bh_posn_y,bh_posn_z);

  cout << "N_varpi: " << N_varpi << endl;
  cout << "N_phi: " << N_phi << endl;
  cout << "N_Z: " << N_Z << endl;
  cout << "varpi_min: " << varpi_min << endl;
  cout << "dvarpi: "<< dvarpi <<endl;
  cout << "dz: " << dz << endl;



  ofstream outfile1;
  char filenamestring[100];
  sprintf(filenamestring,"alpha_1d_t%10.6e.dat",time);
  std::ostringstream filename1;
  filename1 << filenamestring << ends;
  outfile1.open(filename1.str().c_str(),ios::out | ios::ate);

  outfile1 << endl;
  outfile1 << "#1d alpha_visc" << endl;
   
  outfile1 << endl;
  outfile1 << endl;
  outfile1 << "#Time =" << time << endl;

  //outfile1 << "#varpi   alpha_hydro   alpha_em" << endl;
  outfile1 << "#varpi \t alpha_hydro \t alpha_em \t <T_visc^r^phi/P> " << endl;


  double PI = 3.14159265358979323846;
  double dphi = 2.0 * PI / N_phi;
  int n=0;
  
  cout << "*************" << endl;


  //NOTE: Makeing this openMP seems to break things.  Is this unavoidable?
  // ROMAN: Well it should work, but we have to be careful with private variables that are declared before the pragma statements
  //  #pragma omp parallel for
  for (int i = 1; i <= N_varpi; i++){
    double varpi = varpi_min+(i-0.5)*dvarpi;
    double hydro_T_rph_avgL = 0.0;
    double em_T_rph_avgL = 0.0;
    double visc_T_rph_avgL = 0.0;
    double P_avgL = 0.0;

    for (int j = 1; j <= N_phi; j++) {
      double phiangle = (j - 0.5)*dphi;
      double xL = varpi*cos(phiangle)  + bh_posn_x;
      double yL = varpi*sin(phiangle)  + bh_posn_y;
      for (int k = 1; k <= N_Z; k++) {
	double zL = dz*(k-0.5) +  bh_posn_z;
	double rL = sqrt(xL*xL+yL*yL+zL*zL);
	
	double Lx_r = xL/rL;
	double Ly_r = yL/rL;
	double Lz_r = zL/rL;

	double Lx_th = zL*xL/sqrt(xL*xL+yL*yL);	// ==  z cos(phi)?
	double Ly_th = zL*yL/sqrt(xL*xL+yL*yL);	// ==  z sin(phi)?
	double Lz_th = -sqrt(xL*xL+yL*yL);	// == -r sin(theta) ?

	double Lx_ph = -yL;//ok
	double Ly_ph = xL;//ok
	double Lz_ph = 0.0;//ok

	double Lr_x = xL/rL;//ok
	double Lr_y = yL/rL;//ok
	double Lr_z = zL/rL;//ok

	double Lth_x = zL*xL/sqrt(xL*xL+yL*yL)/(rL*rL);
	double Lth_y = zL*yL/sqrt(xL*xL+yL*yL)/(rL*rL);
	double Lth_z = -sqrt(xL*xL+yL*yL)/(rL*rL);

	double Lph_x = -yL/(xL*xL + yL*yL);//ok
	double Lph_y =  xL/(xL*xL + yL*yL);//ok
	double Lph_z = 0.0;//ok
	
	double rho_b_intL = rho_b_int[n];    
	double P_intL = P_int[n];
	double h_intL = h_int[n];
	double vx_intL = vx_int[n];
	double vy_intL = vy_int[n];
	double vz_intL = vz_int[n];
	double u0_intL = u0_int[n];
	
	double lapse_intL = lapse_int[n]+1.0;
        double shiftx_intL = shiftx_int[n];
        double shifty_intL = shifty_int[n];
        double shiftz_intL = shiftz_int[n];

	double alphaLm1 = 1.0/lapse_intL;
	double sqrtmfourpi = 1.0/sqrt(4.0*M_PI);

	// ROMAN:
	// All components of small b are nan (uninitialized on the first time steps)
	//	printf("sbt_int[n]=%e,sbx_int[n]=%e,sby_int[n]=%e,sbz_int[n]=%e\n",sbt_int[n],sbx_int[n],sby_int[n],sbz_int[n]);

	double sbt_intL = sbt_int[n]*alphaLm1*sqrtmfourpi;
	double sbx_intL = sbx_int[n]*alphaLm1*sqrtmfourpi;
	double sby_intL = sby_int[n]*alphaLm1*sqrtmfourpi;
	double sbz_intL = sbz_int[n]*alphaLm1*sqrtmfourpi;

	double phi_intL = phi_int[n];
	double Psi4 = exp(4.0*phi_intL);
	double Psi6 = exp(6.0*phi_intL);

	double gxx_intL = gxx_int[n];
	double gxy_intL = gxy_int[n];
	double gxz_intL = gxz_int[n];
	double gyy_intL = gyy_int[n];
	double gyz_intL = gyz_int[n];
	double gzz_intL = gzz_int[n];
	
	double u_xL = u0_intL*Psi4*(gxx_intL*(shiftx_intL+vx_intL) + 
			    gxy_intL*(shifty_intL+vy_intL) + 
			    gxz_intL*(shiftz_intL+vz_intL)); 
	double u_yL = u0_intL*Psi4*(gxy_intL*(shiftx_intL+vx_intL) + 
			    gyy_intL*(shifty_intL+vy_intL) + 
			    gyz_intL*(shiftz_intL+vz_intL));
	double u_zL = u0_intL*Psi4*(gxz_intL*(shiftx_intL+vx_intL) +
				    gyz_intL*(shifty_intL+vy_intL) +
				    gzz_intL*(shiftz_intL+vz_intL));
	
	double u_0L = -(vx_intL*u_xL+vy_intL*u_yL+vz_intL*u_zL+1.0/u0_intL);
     
	double u_rL = Lx_r*u_xL + Ly_r*u_yL + Lz_r*u_zL;
	double u_thL = Lx_th*u_xL + Ly_th*u_yL + Lz_th*u_zL;
	double u_phL = Lx_ph*u_xL + Ly_ph*u_yL + Lz_ph*u_zL;
	
	double gupxx_intL = gupxx_int[n];
	double gupxy_intL = gupxy_int[n];
	double gupxz_intL = gupxz_int[n];
	double gupyy_intL = gupyy_int[n];
	double gupyz_intL = gupyz_int[n];
	double gupzz_intL = gupzz_int[n];

	//note that this is not the real lowered shift since it is lowered by the conformal 3 metric
	double shift_xL = gxx_intL*shiftx_intL+gxy_intL*shifty_intL+gxz_intL*shiftz_intL;
	double shift_yL = gxy_intL*shiftx_intL+gyy_intL*shifty_intL+gyz_intL*shiftz_intL;
	double shift_zL = gxz_intL*shiftx_intL+gyz_intL*shifty_intL+gzz_intL*shiftz_intL;
    
	double g4tt = -lapse_intL*lapse_intL + Psi4*(shift_xL*shiftx_intL+shift_yL*shifty_intL+shift_zL*shiftz_intL);
	double g4tx = Psi4*shift_xL;
	double g4ty = Psi4*shift_yL;
	double g4tz = Psi4*shift_zL;
	double g4xx = Psi4*gxx_intL;
	double g4xy = Psi4*gxy_intL;
	double g4xz = Psi4*gxz_intL;
	double g4yy = Psi4*gyy_intL;
	double g4yz = Psi4*gyz_intL;
	double g4zz = Psi4*gzz_intL;

	// Compute b^2
	double b2 = -(lapse_intL*sbt_intL)*(lapse_intL*sbt_intL) + 
	  Psi4*( gxx_intL*(sbx_intL+shiftx_intL*sbt_intL)*(sbx_intL+shiftx_intL*sbt_intL) +
		 2.0*gxy_intL*(sbx_intL+shiftx_intL*sbt_intL)*(sby_intL+shifty_intL*sbt_intL) +
		 2.0*gxz_intL*(sbx_intL+shiftx_intL*sbt_intL)*(sbz_intL+shiftz_intL*sbt_intL) +
		 gyy_intL*(sby_intL+shifty_intL*sbt_intL)*(sby_intL+shifty_intL*sbt_intL) +
		 2.0*gyz_intL*(sby_intL+shifty_intL*sbt_intL)*(sbz_intL+shiftz_intL*sbt_intL) +
		 gzz_intL*(sbz_intL+shiftz_intL*sbt_intL)*(sbz_intL+shiftz_intL*sbt_intL) );

	double g4tr = Lx_r*g4tx + Ly_r*g4ty + Lz_r*g4tz;
	double g4tth = Lx_th*g4tx+ Ly_th*g4ty + Lz_th*g4tz;
	double g4tph = Lx_ph*g4tx+ Ly_ph*g4ty + Lz_ph*g4tz;
	double g4rr = Lx_r*(Lx_r*g4xx + Ly_r*g4xy + Lz_r*g4xz)
	  + Ly_r*(Lx_r*g4xy + Ly_r*g4yy + Lz_r*g4yz)
	  + Lz_r*(Lx_r*g4xz + Ly_r*g4yz + Lz_r*g4zz);
	double g4rth = Lx_r*(Lx_th*g4xx + Ly_th*g4xy + Lz_th*g4xz)
          + Ly_r*(Lx_th*g4xy + Ly_th*g4yy + Lz_th*g4yz)
          + Lz_r*(Lx_th*g4xz + Ly_th*g4yz + Lz_th*g4zz);
	double g4rph = Lx_r*(Lx_ph*g4xx + Ly_ph*g4xy + Lz_ph*g4xz)
          + Ly_r*(Lx_ph*g4xy + Ly_ph*g4yy + Lz_ph*g4yz)
          + Lz_r*(Lx_ph*g4xz + Ly_ph*g4yz + Lz_ph*g4zz);
	double g4thth = Lx_th*(Lx_th*g4xx + Ly_th*g4xy + Lz_th*g4xz)
          + Ly_th*(Lx_th*g4xy + Ly_th*g4yy + Lz_th*g4yz)
          + Lz_th*(Lx_th*g4xz + Ly_th*g4yz + Lz_th*g4zz);
	double g4thph = Lx_ph*(Lx_th*g4xx + Ly_th*g4xy + Lz_th*g4xz)
          + Ly_ph*(Lx_th*g4xy + Ly_th*g4yy + Lz_th*g4yz)
          + Lz_ph*(Lx_th*g4xz + Ly_th*g4yz + Lz_th*g4zz);
	double g4phph = Lx_ph*(Lx_ph*g4xx + Ly_ph*g4xy + Lz_ph*g4xz)
          + Ly_ph*(Lx_ph*g4xy + Ly_ph*g4yy + Lz_ph*g4yz)
          + Lz_ph*(Lx_ph*g4xz + Ly_ph*g4yz + Lz_ph*g4zz);

	//better check that gupij is actually set correctly
	double g4uptt = -1.0/(lapse_intL*lapse_intL);
	double g4uptx = 1.0/(lapse_intL*lapse_intL)*shiftx_intL;
	double g4upty = 1.0/(lapse_intL*lapse_intL)*shifty_intL;
	double g4uptz = 1.0/(lapse_intL*lapse_intL)*shiftz_intL;
	double g4upxx = gupxx_intL/Psi4-1.0/(lapse_intL*lapse_intL)*shiftx_intL*shiftx_intL;
	double g4upxy = gupxy_intL/Psi4-1.0/(lapse_intL*lapse_intL)*shiftx_intL*shifty_intL;
	double g4upxz = gupxz_intL/Psi4-1.0/(lapse_intL*lapse_intL)*shiftx_intL*shiftz_intL;
	double g4upyy = gupyy_intL/Psi4-1.0/(lapse_intL*lapse_intL)*shifty_intL*shifty_intL;
	double g4upyz = gupyz_intL/Psi4-1.0/(lapse_intL*lapse_intL)*shifty_intL*shiftz_intL;
	double g4upzz = gupzz_intL/Psi4-1.0/(lapse_intL*lapse_intL)*shiftz_intL*shiftz_intL;

	double rho_starL = lapse_intL*Psi6*u0_intL*rho_b_intL;
    
	double ett  = u0_intL;
	double etr  = u0_intL*( Lr_x*vx_intL + Lr_y *vy_intL + Lr_z *vz_intL); 
	double etth = u0_intL*(Lth_x*vx_intL + Lth_y*vy_intL + Lth_z*vz_intL);
	double etph = u0_intL*(Lph_x*vx_intL + Lph_y*vy_intL + Lph_z*vz_intL);
	  
	double etx = Lx_r*etr + Lx_th*etth + Lx_ph*etph;
	double ety = Ly_r*etr + Ly_th*etth + Ly_ph*etph;;
	double etz = Lz_r*etr + Lz_th*etth + Lz_ph*etph;;
	
	double u_phiL = Lx_ph*u_xL + Ly_ph*u_yL + Lz_ph*u_zL;
	
	double epht = u_phiL;	// eq.(7)
	double ephr = 0.0;	// eq.(7)
	double ephth = 0.0;	// eq.(7)
	double ephph = -u_0L;	// eq.(7)

	double NORM_phi = sqrt(epht*epht*g4tt+ephr*ephr*g4rr+ephth*ephth*g4thth+ephph*ephph*g4phph
			       + 2.0*(epht*ephr*g4tr + epht*ephth*g4tth + epht*ephph*g4tph
				      + ephr*ephth*g4rth + ephr*ephph*g4rph + ephth*ephph*g4thph));

	if(isnan(NORM_phi)){
	  //printf("epht: %e ephph: %e g4tt: %e  g4rr: %e  g4thth: %e  g4phph: %e g4tr: %e  g4tth: %e  g4tph: %e  g4rth: %e  g4rph: %e  g4thph: %e\n",epht,ephph,g4tt, g4rr, g4thth, g4phph,g4tr, g4tth, g4tph, g4rth, g4rph, g4thph);
	  // printf("NORM_phi = nan, xyz, %e %e %e rho_b_intL: %e\n",xL,yL,zL,rho_b_intL);
	}

	epht = epht/NORM_phi;
	ephr = ephr/NORM_phi;
	ephth = ephth/NORM_phi;
	ephph = ephph/NORM_phi;

	double ephx = Lx_r*ephr + Lx_th*ephth + Lx_ph*ephph;
        double ephy = Ly_r*ephr + Ly_th*ephth + Ly_ph*ephph;
        double ephz = Lz_r*ephr + Lz_th*ephth + Lz_ph*ephph;

	double eph_t = g4tt*epht + g4tx*ephx + g4ty*ephy + g4tz*ephz;
	double eph_x = g4tx*epht + g4xx*ephx + g4xy*ephy + g4xz*ephz;
	double eph_y = g4ty*epht + g4xy*ephx + g4yy*ephy + g4yz*ephz;
	double eph_z = g4tz*epht + g4xz*ephx + g4yz*ephy + g4zz*ephz;

	double A = epht*g4tt + ephph*g4tph;   // eq.(8)
	double B = epht*g4tr + ephph*g4rph;   // eq.(9)
	double C = epht*g4tph + ephph*g4phph; // eq.(10)
	
	double ert = u_phL*B-u_rL*C;	// eq. (12)
	double err = u_0L*C-u_phL*A;	// eq. (12)
	double erth = 0.0;		// eq. (12)
	double erph = u_rL*A-u_0L*B;	// eq. (12)

	double NORM_r = sqrt(ert*ert*g4tt+err*err*g4rr+erth*erth*g4thth+erph*erph*g4phph
			     + 2.0*(ert*err*g4tr + ert*erth*g4tth + ert*erph*g4tph
				    + err*erth*g4rth + err*erph*g4rph + erth*erph*g4thph));

	ert = ert/NORM_r;
        err = err/NORM_r;
        erth = erth/NORM_r;
        erph = erph/NORM_r;

	double erx = Lx_r*err + Lx_th*erth + Lx_ph*erph;
        double ery = Ly_r*err + Ly_th*erth + Ly_ph*erph;
        double erz = Lz_r*err + Lz_th*erth + Lz_ph*erph;

	double er_t = g4tt*ert + g4tx*erx + g4ty*ery + g4tz*erz;
	double er_x = g4tx*ert + g4xx*erx + g4xy*ery + g4xz*erz;
	double er_y = g4ty*ert + g4xy*erx + g4yy*ery + g4yz*erz;
        double er_z = g4tz*ert + g4xz*erx + g4yz*ery + g4zz*erz;

	// REF astro-ph/0503420v2 eq (33)
	double hydro_Ttt = rho_b_intL*h_intL*u0_intL*u0_intL+P_intL*g4uptt;//ok
	double hydro_Ttx = rho_b_intL*h_intL*u0_intL*u0_intL*vx_intL+P_intL*g4uptx;
	double hydro_Tty = rho_b_intL*h_intL*u0_intL*u0_intL*vy_intL+P_intL*g4upty;
	double hydro_Ttz = rho_b_intL*h_intL*u0_intL*u0_intL*vz_intL+P_intL*g4uptz;
	double hydro_Txx = rho_b_intL*h_intL*u0_intL*u0_intL*vx_intL*vx_intL+P_intL*g4upxx;
	double hydro_Txy = rho_b_intL*h_intL*u0_intL*u0_intL*vx_intL*vy_intL+P_intL*g4upxy;
	double hydro_Txz = rho_b_intL*h_intL*u0_intL*u0_intL*vx_intL*vz_intL+P_intL*g4upxz;
	double hydro_Tyy = rho_b_intL*h_intL*u0_intL*u0_intL*vy_intL*vy_intL+P_intL*g4upyy;
	double hydro_Tyz = rho_b_intL*h_intL*u0_intL*u0_intL*vy_intL*vz_intL+P_intL*g4upyz;
	double hydro_Tzz = rho_b_intL*h_intL*u0_intL*u0_intL*vz_intL*vz_intL+P_intL*g4upzz;

	// REF astro-ph/0503420v2 eq (32) or eq (33)
	double em_Ttt = b2*(u0_intL*u0_intL+0.5*g4uptt)-sbt_intL*sbt_intL;//ok
	double em_Ttx = b2*(u0_intL*u0_intL*vx_intL+0.5*g4uptx)-sbt_intL*sbx_intL;//ok
	double em_Tty = b2*(u0_intL*u0_intL*vy_intL+0.5*g4upty)-sbt_intL*sby_intL;//ok
	double em_Ttz = b2*(u0_intL*u0_intL*vz_intL+0.5*g4uptz)-sbt_intL*sbz_intL;//ok
	double em_Txx = b2*(u0_intL*u0_intL*vx_intL*vx_intL+0.5*g4upxx)-sbx_intL*sbx_intL;//ok
	double em_Txy = b2*(u0_intL*u0_intL*vx_intL*vy_intL+0.5*g4upxy)-sbx_intL*sby_intL;//ok
	double em_Txz = b2*(u0_intL*u0_intL*vx_intL*vz_intL+0.5*g4upxz)-sbx_intL*sbz_intL;//ok
	double em_Tyy = b2*(u0_intL*u0_intL*vy_intL*vy_intL+0.5*g4upyy)-sby_intL*sby_intL;//ok
	double em_Tyz = b2*(u0_intL*u0_intL*vy_intL*vz_intL+0.5*g4upyz)-sby_intL*sbz_intL;//ok
	double em_Tzz = b2*(u0_intL*u0_intL*vz_intL*vz_intL+0.5*g4upzz)-sbz_intL*sbz_intL;//ok

	int orthonormal = 1;// use orthonormal tetrad
	double hydro_T_orthonormal_rph;
	double em_T_orthonormal_rph;
	double hydro_T_rph;
	double em_T_rph;
	double visc_T_rph;

	if (orthonormal == 1){
       	hydro_T_rph = er_t*(eph_t*hydro_Ttt + eph_x*hydro_Ttx + eph_y*hydro_Tty + eph_z*hydro_Ttz)
	  + er_x*(eph_t*hydro_Ttx + eph_x*hydro_Txx + eph_y*hydro_Txy + eph_z*hydro_Txz)
	  + er_y*(eph_t*hydro_Tty + eph_x*hydro_Txy + eph_y*hydro_Tyy + eph_z*hydro_Tyz)
	  + er_z*(eph_t*hydro_Ttz + eph_x*hydro_Txz + eph_y*hydro_Tyz + eph_z*hydro_Tzz);
	
	if (rho_b_intL>alpha_rho_cut_off){
	em_T_rph = er_t*(eph_t*em_Ttt + eph_x*em_Ttx + eph_y*em_Tty + eph_z*em_Ttz)
          + er_x*(eph_t*em_Ttx + eph_x*em_Txx + eph_y*em_Txy + eph_z*em_Txz)
          + er_y*(eph_t*em_Tty + eph_x*em_Txy + eph_y*em_Tyy + eph_z*em_Tyz)
          + er_z*(eph_t*em_Ttz + eph_x*em_Txz + eph_y*em_Tyz + eph_z*em_Tzz);
	} else{
	  em_T_rph = 0.0;
	}
	  //printf("em_Txx: %e em_Txy: %e em_Txz: %e em_Tyy: %e em_Tyz: %e em_Tzz: %e \n",em_Txx,em_Txy,em_Txz,em_Tyy,em_Tyz,em_Tzz);
	//printf("er_t: %e er_x: %e er_y: %e er_z: %e eph_t: %e eph_x: %e eph_y: %e eph_z: %e\n",er_t,er_x,er_y,er_z,eph_t,eph_x,eph_y,eph_z);
	
	visc_T_rph = 0.0;

        hydro_T_rph_avgL = hydro_T_rph_avgL + hydro_T_rph*dz*dphi*u0_intL*lapse_intL*Psi6;
	//printf("em_T_rph_avgL: %e em_T_rph: %e dz: %e dphi: %e u0_intL: %e lapse_intL: %e Psi6: %e\n",em_T_rph_avgL,em_T_rph,dz,dphi,u0_intL,lapse_intL,Psi6);
	if (isnan(em_T_rph)){
	  //printf("em_Ttt: %e em_Ttx: %e em_Tty: %e em_Ttz: %e em_Txx: %e em_Txy: %e em_xz: %e em_Tyy: %e em_Tyz: %e em_Tzz: %e \n",em_Ttt,em_Ttx,em_Tty,em_Ttz,em_Txx,em_Txy,em_Txz,em_Tyy,em_Tyz,em_Tzz); 
	  //printf("em_T_rph_avgL: %e em_T_rph: %e er_t: %e er_x: %e er_y: %e er_z: %e eph_t: %e eph_x: %e eph_y: %e eph_z: %e epht: %e ephph: %e\n",em_T_rph_avgL,em_T_rph,er_t,er_x,er_y,er_z,eph_t,eph_x,eph_y,eph_z,epht,ephph); 
	  //printf("ert: %e err: %e  erth: %e  erph: %e u_rL: %e u_0L: %e u_phL: %e A: %e B: %e C: %e\n",ert,err,erth,erph,u_rL,u_0L,u_phL,A,B,C);
	}
	em_T_rph_avgL = em_T_rph_avgL + em_T_rph*dz*dphi*u0_intL*lapse_intL*Psi6;
        visc_T_rph_avgL = visc_T_rph_avgL + visc_T_rph*dz*dphi*u0_intL*lapse_intL*Psi6;
	}
	else {// non-orthonormal
     	hydro_T_rph = 
	  + Lr_x*(Lph_x*hydro_Txx + Lph_y*hydro_Txy + Lph_z*hydro_Txz)
	  + Lr_y*(Lph_x*hydro_Txy + Lph_y*hydro_Tyy + Lph_z*hydro_Tyz)
	  + Lr_z*(Lph_x*hydro_Txz + Lph_y*hydro_Tyz + Lph_z*hydro_Tzz);
	
	em_T_rph = 
          + Lr_x*(Lph_x*em_Txx + Lph_y*em_Txy + Lph_z*em_Txz)
          + Lr_y*(Lph_x*em_Txy + Lph_y*em_Tyy + Lph_z*em_Tyz)
          + Lr_z*(Lph_x*em_Txz + Lph_y*em_Tyz + Lph_z*em_Tzz);

	// viscous stress
	visc_T_rph = 0.0;

        hydro_T_rph_avgL = hydro_T_rph_avgL + hydro_T_rph*dz*dphi*u0_intL*lapse_intL*Psi6;
        //printf("em_T_rph_avgL: %e em_T_rph: %e dz: %e dphi: %e u0_intL: %e lapse_intL: %e Psi6: %e\n",em_T_rph_avgL,em_T_rph,dz,dphi,u0_intL,lapse_intL,Psi6);
	em_T_rph_avgL = em_T_rph_avgL + em_T_rph*dz*dphi*u0_intL*lapse_intL*Psi6;

	visc_T_rph_avgL = visc_T_rph_avgL + visc_T_rph*dz*dphi*u0_intL*lapse_intL*Psi6;
	}

        P_avgL = P_avgL + P_intL*dz*dphi*u0_intL*lapse_intL*Psi6;

	n=n+1;
      }
    }
    //    printf("P_avgL: %e em_T_rph_avgL %e\n",P_avgL,em_T_rph_avgL);

      outfile1 << varpi << "\t" << hydro_T_rph_avgL/P_avgL << "\t" << em_T_rph_avgL/P_avgL << "\t"  << "\t" << visc_T_rph_avgL/P_avgL << endl;
  }
  
  //close file
  outfile1.close();
}
extern "C" void CCTK_FCALL CCTK_FNAME(bhns_compute_alpha_IO)
  (const cGH **cctkGH, double &time,
   double &varpi_min,double &dvarpi,double &dz,
   int &N_varpi,int &N_phi,int &N_Z,
   double *rho_b_int,double *P_int,double *h_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
   double *sbt_int,double *sbx_int,double *sby_int,double *sbz_int,
   double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
   double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
   double *gyy_int,double *gyz_int, double *gzz_int,
   double *gupxx_int,double *gupxy_int,double *gupxz_int,
   double *gupyy_int,double *gupyz_int, double *gupzz_int,
   double &bh_posn_x, double &bh_posn_y, double &bh_posn_z, double &alpha_rho_cut_off)
   //   double *sigma_xx_int,double *sigma_xy_int,double *sigma_xz_int,
   //   double *sigma_yy_int,double *sigma_yz_int,double *sigma_zz_int,
   //   double *eta_vis_int)
{  
  bhns_compute_alpha_IO(*cctkGH, time,
		   varpi_min,dvarpi,dz,
		   N_varpi,N_phi,N_Z,
		   rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,
		   sbt_int,sbx_int,sby_int,sbz_int,
		   lapse_int,shiftx_int,shifty_int,shiftz_int,
		   phi_int,gxx_int,gxy_int,gxz_int,
		   gyy_int,gyz_int, gzz_int,
		   gupxx_int,gupxy_int,gupxz_int,
		   gupyy_int,gupyz_int, gupzz_int,
		   bh_posn_x, bh_posn_y, bh_posn_z, alpha_rho_cut_off);
  //		   sigma_xx_int,sigma_xy_int,sigma_xz_int,
  //		   sigma_yy_int,sigma_yz_int,sigma_zz_int,
  //		   eta_vis_int);
}


//  LocalWords:  intL

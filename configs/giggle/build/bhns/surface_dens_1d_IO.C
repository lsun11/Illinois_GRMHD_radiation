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

extern "C" void CCTK_FCALL bhns_surf_dens_1d_io_
  (const cGH **cctkGH, double &time,
   double &varpi_min,double &dvarpi,double &dz,
   int &N_varpi,int &N_phi,int &N_Z,
   double *rho_b_int,double *P_int,double *h_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
   double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
   double *lapsex_int,double *shiftxx_int,double *shiftyx_int,double *shiftzx_int,
   double *lapsey_int,double *shiftxy_int,double *shiftyy_int,double *shiftzy_int,
   double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
   double *gyy_int,double *gyz_int, double *gzz_int,
   double *gupxx_int,double *gupxy_int,double *gupxz_int,
   double *gupyy_int,double *gupyz_int, double *gupzz_int,
   double *phix_int,double *gxxx_int,double *gxyx_int,double *gxzx_int,
   double *gyyx_int,double *gyzx_int, double *gzzx_int,
   double *phiy_int,double *gxxy_int,double *gxyy_int,double *gxzy_int,
   double *gyyy_int,double *gyzy_int, double *gzzy_int, double &bh_posn_x, double &bh_posn_y, double &bh_posn_z);

extern "C" void bhns_surf_dens_1d_IO(const cGH *cctkGH, double &time,
				     double &varpi_min,double &dvarpi,double &dz,
				     int &N_varpi,int &N_phi,int &N_Z,
				     double *rho_b_int,double *P_int,double *h_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
				     double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
				     double *lapsex_int,double *shiftxx_int,double *shiftyx_int,double *shiftzx_int,
				     double *lapsey_int,double *shiftxy_int,double *shiftyy_int,double *shiftzy_int,
				     double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
				     double *gyy_int,double *gyz_int, double *gzz_int,
				     double *gupxx_int,double *gupxy_int,double *gupxz_int,
				     double *gupyy_int,double *gupyz_int, double *gupzz_int,
				     double *phix_int,double *gxxx_int,double *gxyx_int,double *gxzx_int,
				     double *gyyx_int,double *gyzx_int, double *gzzx_int,
				     double *phiy_int,double *gxxy_int,double *gxyy_int,double *gxzy_int,
				     double *gyyy_int,double *gyzy_int, double *gzzy_int, 
				     double &bh_posn_x, double &bh_posn_y, double &bh_posn_z){
 

  printf("syrface- BH position %e %e %e\n",bh_posn_x,bh_posn_y,bh_posn_z);

  ofstream outfile1;
  char filenamestring[100];
  sprintf(filenamestring,"surf_dens_1d_t%10.6e.dat",time);
  std::ostringstream filename1;
  filename1 << filenamestring << ends;
  outfile1.open(filename1.str().c_str(),ios::out | ios::ate);

  ofstream outfile2;
  sprintf(filenamestring,"surf_dens_polar_t%10.6e.dat",time);
  std::ostringstream filename2;
  filename2 << filenamestring << ends;
  outfile2.open(filename2.str().c_str(),ios::out | ios::ate);

  ofstream outfile3;
  sprintf(filenamestring,"surf_dens_polar_avg_t%10.6e.dat",time);
  std::ostringstream filename3;
  filename3 << filenamestring << ends;
  outfile3.open(filename3.str().c_str(),ios::out | ios::ate);

  outfile1 << endl;
  outfile1 << "#1d surface density" << endl;
  outfile2 << endl;
  outfile2 << "#polar surface density" << endl;
  outfile3 << endl;
  outfile3 << "#polar surface avg density" << endl;
   
  outfile1 << endl;
  outfile1 << endl;
  outfile1 << "#Time =" << time << endl;
  outfile2 << endl;
  outfile2 << endl;
  outfile2 << "#Time =" << time << endl;
  outfile3 << endl;
  outfile3 << endl;
  outfile3 << "#Time =" << time << endl;

  outfile1 << "#varpi   Sigma   torque_density" << endl;
  outfile2 << "#varpi   phiangle    Sigma   torque_density" << endl;
  outfile3 << "#varpi   phiangle    Sigma   torque_density" << endl;

  double PI = 3.14159265358979323844;
  double dphi = 2.0 * PI / N_phi;
  int n=0; // ROMAN: FIXME: if openMP is used in below loop this should be private
  
  cout << "*************" << endl;

  //NOTE: Makeing this openMP seems to break things.  Is this unavoidable?
  //  #pragma omp parallel for 
  // ROMAN: the below line should fix it: //but I think it does not
  //#pragma omp parallel for private(n,Sigma_avgL,torque_densL,varpi)
  for (int i = 1; i <= N_varpi; i++){
    double varpi = varpi_min+(i-0.5)*dvarpi;
    double Sigma_avgL = 0.0;
    double torque_densL = 0.0;
    /*double torque_dens1L = 0.0;
    double torque_dens2L = 0.0;
    double torque_dens3L = 0.0;
    double torque_dens4L = 0.0;
    double torque_dens_oldL =0.0;
    */

    for (int j = 1; j <= N_phi; j++) {
      double phiangle = (j - 0.5)*dphi;
      double SigmaL = 0.0;
      double xL = varpi*cos(phiangle);
      double yL = varpi*sin(phiangle);
      for (int k = 1; k <= N_Z; k++) {
	double ZL = dz*(k-0.5) +  bh_posn_z;
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
  
	double phi_intL = phi_int[n];
	double Psi4 = exp(4.0*phi_intL);
	double Psi6 = exp(6.0*phi_intL);


	SigmaL = SigmaL + rho_b_intL*lapse_intL*Psi6*u0_intL*dz;
	Sigma_avgL = Sigma_avgL + rho_b_intL*lapse_intL*Psi6*u0_intL*dz*dphi/(2.0*PI);

	double phix_intL = phix_int[n];
	double phiy_intL = phiy_int[n];
	
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
     
	double gupxx_intL = gupxx_int[n];
	double gupxy_intL = gupxy_int[n];
	double gupxz_intL = gupxz_int[n];
	double gupyy_intL = gupyy_int[n];
	double gupyz_intL = gupyz_int[n];
	double gupzz_intL = gupzz_int[n];
	
	double gxxx_intL = gxxx_int[n];
	double gxyx_intL = gxyx_int[n];
	double gxzx_intL = gxzx_int[n];
	double gyyx_intL = gyyx_int[n];
	double gyzx_intL = gyzx_int[n];
	double gzzx_intL = gzzx_int[n];
	
	double gxxy_intL = gxxy_int[n];
	double gxyy_intL = gxyy_int[n];
	double gxzy_intL = gxzy_int[n];
	double gyyy_intL = gyyy_int[n];
	double gyzy_intL = gyzy_int[n];
	double gzzy_intL = gzzy_int[n];

	double lapsex_intL = lapsex_int[n];
	double shiftxx_intL = shiftxx_int[n];
	double shiftyx_intL = shiftyx_int[n];
	double shiftzx_intL = shiftzx_int[n];

	double lapsey_intL = lapsey_int[n];
	double shiftxy_intL = shiftxy_int[n];
	double shiftyy_intL = shiftyy_int[n];
	double shiftzy_intL = shiftzy_int[n];

	double g4xxx = Psi4*(4.0*phix_intL*gxx_intL+gxxx_intL);
	double g4xxy = Psi4*(4.0*phiy_intL*gxx_intL+gxxy_intL);

	double g4xyx = Psi4*(4.0*phix_intL*gxy_intL+gxyx_intL);
	double g4xyy = Psi4*(4.0*phiy_intL*gxy_intL+gxyy_intL);

	double g4xzx = Psi4*(4.0*phix_intL*gxz_intL+gxzx_intL);
	double g4xzy = Psi4*(4.0*phiy_intL*gxz_intL+gxzy_intL);

	double g4yyx = Psi4*(4.0*phix_intL*gyy_intL+gyyx_intL);
	double g4yyy = Psi4*(4.0*phiy_intL*gyy_intL+gyyy_intL);

	double g4yzx = Psi4*(4.0*phix_intL*gyz_intL+gyzx_intL);
	double g4yzy = Psi4*(4.0*phiy_intL*gyz_intL+gyzy_intL);

	double g4zzx = Psi4*(4.0*phix_intL*gzz_intL+gzzx_intL);
	double g4zzy = Psi4*(4.0*phiy_intL*gzz_intL+gzzy_intL);
    
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

	double dx_gxj_shiftj = gxxx_intL*shiftx_intL+gxyx_intL*shifty_intL*gxzx_intL*shiftz_intL;
	double dx_gyj_shiftj = gxyx_intL*shiftx_intL+gyyx_intL*shifty_intL*gyzx_intL*shiftz_intL;
	double dx_gzj_shiftj = gxzx_intL*shiftx_intL+gyzx_intL*shifty_intL*gzzx_intL*shiftz_intL;
    
	double dy_gxj_shiftj = gxxy_intL*shiftx_intL+gxyy_intL*shifty_intL*gxzy_intL*shiftz_intL;
	double dy_gyj_shiftj = gxyy_intL*shiftx_intL+gyyy_intL*shifty_intL*gyzy_intL*shiftz_intL;
	double dy_gzj_shiftj = gxzy_intL*shiftx_intL+gyzy_intL*shifty_intL*gzzy_intL*shiftz_intL;
    
	double dx_shiftj_gxj = gxx_intL*shiftxx_intL+gxy_intL*shiftyx_intL+gxz_intL*shiftzx_intL;
	double dx_shiftj_gyj = gxy_intL*shiftxx_intL+gyy_intL*shiftyx_intL+gyz_intL*shiftzx_intL;
	double dx_shiftj_gzj = gxz_intL*shiftxx_intL+gyz_intL*shiftyx_intL+gzz_intL*shiftzx_intL;
    
	double dy_shiftj_gxj = gxx_intL*shiftxy_intL+gxy_intL*shiftyy_intL+gxz_intL*shiftzy_intL;
	double dy_shiftj_gyj = gxy_intL*shiftxy_intL+gyy_intL*shiftyy_intL+gyz_intL*shiftzy_intL;
	double dy_shiftj_gzj = gxz_intL*shiftxy_intL+gyz_intL*shiftyy_intL+gzz_intL*shiftzy_intL;
    
	double g4txx = g4xxx*shiftx_intL+g4xyx*shifty_intL+g4xzx*shiftz_intL
	  +g4xx*shiftxx_intL+g4xy*shiftyx_intL+g4xz*shiftzx_intL;
	double g4tyx = g4xyx*shiftx_intL+g4yyx*shifty_intL+g4yzx*shiftz_intL
	  +g4xy*shiftxx_intL+g4yy*shiftyx_intL+g4yz*shiftzx_intL;
	double g4tzx = g4xzx*shiftx_intL+g4yzx*shifty_intL+g4zzx*shiftz_intL
	  +g4xz*shiftxx_intL+g4yz*shiftyx_intL+g4zz*shiftzx_intL;
    
	double g4txy = g4xxy*shiftx_intL+g4xyy*shifty_intL+g4xzy*shiftz_intL
	  +g4xx*shiftxy_intL+g4xy*shiftyy_intL+g4xz*shiftzy_intL;
	double g4tyy = g4xyy*shiftx_intL+g4yyy*shifty_intL+g4yzy*shiftz_intL
	  +g4xy*shiftxy_intL+g4yy*shiftyy_intL+g4yz*shiftzy_intL;
	double g4tzy = g4xzy*shiftx_intL+g4yzy*shifty_intL+g4zzy*shiftz_intL
	  +g4xz*shiftxy_intL+g4yz*shiftyy_intL+g4zz*shiftzy_intL;
  
	double g4ttx = -2.0*lapse_intL*lapsex_intL +  
	  g4xxx*shiftx_intL*shiftx_intL+g4yyx*shifty_intL*shifty_intL+g4zzx*shiftz_intL*shiftz_intL+
	  2.0*(g4xyx*shiftx_intL*shifty_intL+g4xzx*shiftx_intL*shiftz_intL+g4yzx*shifty_intL*shiftz_intL)+
	  2.0*Psi4*(shift_xL*shiftxx_intL+shift_yL*shiftyx_intL+shift_zL*shiftzx_intL);
	double g4tty = -2.0*lapse_intL*lapsey_intL +  
	  g4xxy*shiftx_intL*shiftx_intL+g4yyy*shifty_intL*shifty_intL+g4zzy*shiftz_intL*shiftz_intL+
	  2.0*(g4xyy*shiftx_intL*shifty_intL+g4xzy*shiftx_intL*shiftz_intL+g4yzy*shifty_intL*shiftz_intL)+
	  2.0*Psi4*(shift_xL*shiftxy_intL+shift_yL*shiftyy_intL+shift_zL*shiftzy_intL);
    
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
    
	double integrand1 = rho_starL*h_intL*(-vy_intL*u_xL+vx_intL*u_yL);
	double integrandtt = (-yL*g4ttx+xL*g4tty)*(rho_starL*h_intL*u0_intL+lapse_intL*Psi6*P_intL*g4uptt);
	double integrandtx = (-yL*g4txx+xL*g4txy)*(rho_starL*h_intL*vx_intL*u0_intL+lapse_intL*Psi6*P_intL*shiftx_intL);
	double integrandty = (-yL*g4tyx+xL*g4tyy)*(rho_starL*h_intL*vy_intL*u0_intL+lapse_intL*Psi6*P_intL*shifty_intL);
	double integrandtz = (-yL*g4tzx+xL*g4tzy)*(rho_starL*h_intL*vz_intL*u0_intL+lapse_intL*Psi6*P_intL*shiftz_intL);
	double integrandxx = (-yL*g4xxx+xL*g4xxy)*(rho_starL*h_intL*vx_intL*vx_intL*u0_intL+lapse_intL*Psi6*P_intL*g4upxx);
	double integrandyy = (-yL*g4yyx+xL*g4yyy)*(rho_starL*h_intL*vy_intL*vy_intL*u0_intL+lapse_intL*Psi6*P_intL*g4upyy);
	double integrandzz = (-yL*g4zzx+xL*g4zzy)*(rho_starL*h_intL*vz_intL*vz_intL*u0_intL+lapse_intL*Psi6*P_intL*g4upzz);
	double integrandxy = (-yL*g4xyx+xL*g4xyy)*(rho_starL*h_intL*vx_intL*vy_intL*u0_intL+lapse_intL*Psi6*P_intL*g4upxy);
	double integrandxz = (-yL*g4xzx+xL*g4xzy)*(rho_starL*h_intL*vx_intL*vz_intL*u0_intL+lapse_intL*Psi6*P_intL*g4upxz);
	double integrandyz = (-yL*g4yzx+xL*g4yzy)*(rho_starL*h_intL*vy_intL*vz_intL*u0_intL+lapse_intL*Psi6*P_intL*g4upyz);

	double total_integrand = integrand1+0.5*(integrandtt+2.0*(integrandtx+integrandty+integrandtz)+
						       integrandxx+integrandyy+integrandzz+2.0*(integrandxy+integrandxz+integrandyz));
	torque_densL = torque_densL+total_integrand*varpi*dz*dphi; 
	n=n+1; 
      }
      //insert output here
      outfile2 << varpi << " " << phiangle << " " << SigmaL << " " << torque_densL << endl;
    }
    outfile2 << endl;
    outfile1 << varpi << " " << Sigma_avgL << " " <<torque_densL << endl; 
    
    for (int j1 = 1; j1 <= N_phi; j1++) {
      double phiangle = (j1 - 0.5)*dphi;
      outfile3 << varpi << " " << phiangle << " " << Sigma_avgL << " " << torque_densL << endl;
      //insert output here
    }
    outfile3 << endl;
  }
  

  //close files
  outfile1.close();
  outfile2.close();
  outfile3.close();
}
extern "C" void CCTK_FCALL bhns_surf_dens_1d_io_
  (const cGH **cctkGH, double &time,
   double &varpi_min,double &dvarpi,double &dz,
   int &N_varpi,int &N_phi,int &N_Z,
   double *rho_b_int,double *P_int,double *h_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
   double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
   double *lapsex_int,double *shiftxx_int,double *shiftyx_int,double *shiftzx_int,
   double *lapsey_int,double *shiftxy_int,double *shiftyy_int,double *shiftzy_int,
   double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
   double *gyy_int,double *gyz_int, double *gzz_int,
   double *gupxx_int,double *gupxy_int,double *gupxz_int,
   double *gupyy_int,double *gupyz_int, double *gupzz_int,
   double *phix_int,double *gxxx_int,double *gxyx_int,double *gxzx_int,
   double *gyyx_int,double *gyzx_int, double *gzzx_int,
   double *phiy_int,double *gxxy_int,double *gxyy_int,double *gxzy_int,
   double *gyyy_int,double *gyzy_int, double *gzzy_int,
   double &bh_posn_x, double &bh_posn_y, double &bh_posn_z)
{  
  bhns_surf_dens_1d_IO(*cctkGH, time,
		       varpi_min,dvarpi,dz,
		       N_varpi,N_phi,N_Z,
		       rho_b_int,P_int,h_int,vx_int,vy_int,vz_int,u0_int,
		       lapse_int,shiftx_int,shifty_int,shiftz_int,
		       lapsex_int,shiftxx_int,shiftyx_int,shiftzx_int,
		       lapsey_int,shiftxy_int,shiftyy_int,shiftzy_int,
		       phi_int,gxx_int,gxy_int,gxz_int,
		       gyy_int,gyz_int, gzz_int,
		       gupxx_int,gupxy_int,gupxz_int,
		       gupyy_int,gupyz_int, gupzz_int,
		       phix_int,gxxx_int,gxyx_int,gxzx_int,
		       gyyx_int,gyzx_int, gzzx_int,
		       phiy_int,gxxy_int,gxyy_int,gxzy_int,
		       gyyy_int,gyzy_int, gzzy_int,
		       bh_posn_x, bh_posn_y, bh_posn_z);
}


//  LocalWords:  intL

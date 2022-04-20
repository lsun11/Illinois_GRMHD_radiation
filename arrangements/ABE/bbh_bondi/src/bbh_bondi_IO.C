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

extern "C" void CCTK_FCALL CCTK_FNAME(bbh_bondi_IO)
  (const cGH **cctkGH, double &time,int &iteration, int &gridnum,
   double &output_radmin,double &output_radmax,
   int &output_Nlograd,int &output_Nphi,int &output_Ntheta,
   double *rho_b_int,double *P_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
   double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
   double *lapsex_int,double *shiftxx_int,double *shiftyx_int,double *shiftzx_int,
   double *lapsey_int,double *shiftxy_int,double *shiftyy_int,double *shiftzy_int,
   double *lapsez_int,double *shiftxz_int,double *shiftyz_int,double *shiftzz_int,
   double *lapset_int,double *shiftxt_int,double *shiftyt_int,double *shiftzt_int,
   double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
   double *gyy_int,double *gyz_int, double *gzz_int,
   double *gupxx_int,double *gupxy_int,double *gupxz_int,
   double *gupyy_int,double *gupyz_int, double *gupzz_int,
   double *phix_int,double *gxxx_int,double *gxyx_int,double *gxzx_int,
   double *gyyx_int,double *gyzx_int, double *gzzx_int,
   double *phiy_int,double *gxxy_int,double *gxyy_int,double *gxzy_int,
   double *gyyy_int,double *gyzy_int, double *gzzy_int,
   double *phiz_int,double *gxxz_int,double *gxyz_int,double *gxzz_int,
   double *gyyz_int,double *gyzz_int, double *gzzz_int,
   double *phit_int,double *gxxt_int,double *gxyt_int,double *gxzt_int,
   double *gyyt_int,double *gyzt_int,double *gzzt_int,
   double &BH_posx,double &BH_posy);

extern "C" void bbh_bondi_IO(const cGH *cctkGH, double &time,int &iteration, int &gridnum,
			     double &output_radmin,double &output_radmax,
			     int &output_Nlograd,int &output_Nphi,int &output_Ntheta,
			     double *rho_b_int,double *P_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
			     double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
			     double *lapsex_int,double *shiftxx_int,double *shiftyx_int,double *shiftzx_int,
			     double *lapsey_int,double *shiftxy_int,double *shiftyy_int,double *shiftzy_int,
			     double *lapsez_int,double *shiftxz_int,double *shiftyz_int,double *shiftzz_int,
			     double *lapset_int,double *shiftxt_int,double *shiftyt_int,double *shiftzt_int,
			     double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
			     double *gyy_int,double *gyz_int, double *gzz_int,
			     double *gupxx_int,double *gupxy_int,double *gupxz_int,
			     double *gupyy_int,double *gupyz_int, double *gupzz_int,
			     double *phix_int,double *gxxx_int,double *gxyx_int,double *gxzx_int,
			     double *gyyx_int,double *gyzx_int, double *gzzx_int,
			     double *phiy_int,double *gxxy_int,double *gxyy_int,double *gxzy_int,
			     double *gyyy_int,double *gyzy_int, double *gzzy_int,
			     double *phiz_int,double *gxxz_int,double *gxyz_int,double *gxzz_int,
			     double *gyyz_int,double *gyzz_int, double *gzzz_int,
			     double *phit_int,double *gxxt_int,double *gxyt_int,double *gxzt_int,
			     double *gyyt_int,double *gyzt_int,double *gzzt_int,
			     double &BH_posx,double &BH_posy){
 

  cout << "output_Nlograd: " << output_Nlograd << endl;
  cout << "output_Nphi: " << output_Nphi << endl;
  cout << "output_Ntheta: " << output_Ntheta << endl;
  
  //  printf ("hello 1\n");
  printf("iteration: %d\n",iteration);
  ofstream outfile;
  char filenamestring[100];
  sprintf(filenamestring,"outfile%d.iter%d.dat",gridnum,iteration);
  std::ostringstream filename;filename << filenamestring << ends;
  outfile.open(filename.str().c_str(),ios::out | ios::ate);

  double PI=acos(-1.0);
  double thetamin=-PI/100000.;
  double thetamax=PI/2.0 + PI/100000.;
  double phimin=-PI/100000.;
  double phimax=2.0*PI+PI/100000.;

  double output_logradmin = log(output_radmin);
  //  double output_logradmin2 = log(output_radmin2);
  // double output_logradmin3 = log(output_radmin3);
  double output_logradmax = log(output_radmax);
  // double output_logradmax2 = log(output_radmax2);
  // double output_logradmax3 = log(output_radmax3);


  outfile.write((char *) &time,sizeof(double));
  outfile.write((char *) &BH_posx,sizeof(double));
  outfile.write((char *) &BH_posy,sizeof(double));
  //  outfile.write((char *) &BH2_posx,sizeof(double));
  // outfile.write((char *) &BH2_posy,sizeof(double));
  outfile.write((char *) &output_Nlograd,sizeof(int));
  outfile.write((char *) &output_Nphi,sizeof(int));
  outfile.write((char *) &output_Ntheta,sizeof(int));
  outfile.write((char *) &output_logradmin,sizeof(double));
  // outfile.write((char *) &output_logradmin2,sizeof(double));
  //  outfile.write((char *) &output_logradmin3,sizeof(double));
  outfile.write((char *) &output_logradmax,sizeof(double));
  outfile.write((char *) &thetamin,sizeof(double));
  outfile.write((char *) &thetamax,sizeof(double));
  outfile.write((char *) &phimin,sizeof(double));
  outfile.write((char *) &phimax,sizeof(double));
  //  outfile.write((char *) &output_logradmax2,sizeof(double));
  //  outfile.write((char *) &output_logradmax3,sizeof(double));
  
 //  cout << "output_logradmin1: "<< output_logradmin1 << endl;
//   cout << "output_logradmin2: "<< output_logradmin2 << endl;
//   cout << "output_logradmin3: "<< output_logradmin3 << endl;
//   cout << "output_logradmax1: "<< output_logradmax1 << endl;
//   cout << "output_logradmax2: "<< output_logradmax2 << endl;
//   cout << "output_logradmax3: "<< output_logradmax3 << endl;
  int n=0;
  //  printf ("hello 2\n");
  
  //NOTE: Makeing this openMP seems to break things.  Is this unavoidable?
  //  #pragma omp parallel for
  for (int k = 0; k < output_Nlograd; k++) for (int j = 0; j < output_Nphi; j++) for (int i = 0; i < output_Ntheta; i++) {
    // printf ("test1 i: %d, j: %d, k: %d\n",i,j,k);
    double output_logradL=output_logradmin+(output_logradmax-output_logradmin)/(output_Nlograd-1.0)*(k);
    //    double output_lograd2L=output_logradmin2+(output_logradmax2-output_logradmin2)/(output_Nlograd-1.0)*k;
    //   double output_lograd3L=output_logradmin3+(output_logradmax3-output_logradmin3)/(output_Nlograd-1.0)*k;
    double phiL=phimin+(phimax-phimin)/(output_Nphi-1.0)*(j);
    double thetaL = thetamin+(thetamax-thetamin)/(output_Ntheta-1.0)*(i);
    double output_radL=exp(output_logradL);
  
    double phi_intL = phi_int[n];
    double Psi4 = exp(4.0*phi_intL);

    double phix_intL = phix_int[n];
    double phiy_intL = phiy_int[n];
    double phiz_intL = phiz_int[n];
    double phit_intL = phit_int[n];

    double gxx_intL = gxx_int[n];
    double gxy_intL = gxy_int[n];
    double gxz_intL = gxz_int[n];
    double gyy_intL = gyy_int[n];
    double gyz_intL = gyz_int[n];
    double gzz_intL = gzz_int[n];

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

    double gxxz_intL = gxxz_int[n];
    double gxyz_intL = gxyz_int[n];
    double gxzz_intL = gxzz_int[n];
    double gyyz_intL = gyyz_int[n];
    double gyzz_intL = gyzz_int[n];
    double gzzz_intL = gzzz_int[n];

    double gxxt_intL = gxxt_int[n];
    double gxyt_intL = gxyt_int[n];
    double gxzt_intL = gxzt_int[n];
    double gyyt_intL = gyyt_int[n];
    double gyzt_intL = gyzt_int[n];
    double gzzt_intL = gzzt_int[n];

    double lapse_intL = lapse_int[n]+1.0;
    double shiftx_intL = shiftx_int[n];
    double shifty_intL = shifty_int[n];
    double shiftz_intL = shiftz_int[n];
    
    

    double lapsex_intL = lapsex_int[n];
    double shiftxx_intL = shiftxx_int[n];
    double shiftyx_intL = shiftyx_int[n];
    double shiftzx_intL = shiftzx_int[n];

    double lapsey_intL = lapsey_int[n];
    double shiftxy_intL = shiftxy_int[n];
    double shiftyy_intL = shiftyy_int[n];
    double shiftzy_intL = shiftzy_int[n];

    double lapsez_intL = lapsez_int[n];
    double shiftxz_intL = shiftxz_int[n];
    double shiftyz_intL = shiftyz_int[n];
    double shiftzz_intL = shiftzz_int[n];
 
    double lapset_intL = lapset_int[n];
    double shiftxt_intL = shiftxt_int[n];
    double shiftyt_intL = shiftyt_int[n];
    double shiftzt_intL = shiftzt_int[n];

    double g4xxx = Psi4*(4.0*phix_intL*gxx_intL+gxxx_intL);
    double g4xxy = Psi4*(4.0*phiy_intL*gxx_intL+gxxy_intL);
    double g4xxz = Psi4*(4.0*phiz_intL*gxx_intL+gxxz_intL);
    double g4xxt = Psi4*(4.0*phit_intL*gxx_intL+gxxt_intL);

    // printf("Psi4: %10.6e\n",Psi4);
//     printf("gxx_intL: %10.6e\n",gxx_intL);
//     printf("phit_intL: %10.6e\n",phit_intL);
//     printf("gxxt_intL: %10.6e\n",gxxt_intL);

    double g4xyx = Psi4*(4.0*phix_intL*gxy_intL+gxyx_intL);
    double g4xyy = Psi4*(4.0*phiy_intL*gxy_intL+gxyy_intL);
    double g4xyz = Psi4*(4.0*phiz_intL*gxy_intL+gxyz_intL);
    double g4xyt = Psi4*(4.0*phit_intL*gxy_intL+gxyt_intL);

    double g4xzx = Psi4*(4.0*phix_intL*gxz_intL+gxzx_intL);
    double g4xzy = Psi4*(4.0*phiy_intL*gxz_intL+gxzy_intL);
    double g4xzz = Psi4*(4.0*phiz_intL*gxz_intL+gxzz_intL);
    double g4xzt = Psi4*(4.0*phit_intL*gxz_intL+gxzt_intL);

    double g4yyx = Psi4*(4.0*phix_intL*gyy_intL+gyyx_intL);
    double g4yyy = Psi4*(4.0*phiy_intL*gyy_intL+gyyy_intL);
    double g4yyz = Psi4*(4.0*phiz_intL*gyy_intL+gyyz_intL);
    double g4yyt = Psi4*(4.0*phit_intL*gyy_intL+gyyt_intL);

    double g4yzx = Psi4*(4.0*phix_intL*gyz_intL+gyzx_intL);
    double g4yzy = Psi4*(4.0*phiy_intL*gyz_intL+gyzy_intL);
    double g4yzz = Psi4*(4.0*phiz_intL*gyz_intL+gyzz_intL);
    double g4yzt = Psi4*(4.0*phit_intL*gyz_intL+gyzt_intL);

    double g4zzx = Psi4*(4.0*phix_intL*gzz_intL+gzzx_intL);
    double g4zzy = Psi4*(4.0*phiy_intL*gzz_intL+gzzy_intL);
    double g4zzz = Psi4*(4.0*phiz_intL*gzz_intL+gzzz_intL);
    double g4zzt = Psi4*(4.0*phit_intL*gzz_intL+gzzt_intL);
    
//     printf("g4xxt: %10.6e\n",g4xxt);
//     printf("g4xyt: %10.6e\n",g4xyt);
//     printf("g4xzt: %10.6e\n",g4xzt);
//     printf("g4yyt: %10.6e\n",g4yyt);
//     printf("g4yzt: %10.6e\n",g4yzt);
//     printf("g4zzt: %10.6e\n",g4zzt);

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
    
    double dz_gxj_shiftj = gxxz_intL*shiftx_intL+gxyz_intL*shifty_intL*gxzz_intL*shiftz_intL;
    double dz_gyj_shiftj = gxyz_intL*shiftx_intL+gyyz_intL*shifty_intL*gyzz_intL*shiftz_intL;
    double dz_gzj_shiftj = gxzz_intL*shiftx_intL+gyzz_intL*shifty_intL*gzzz_intL*shiftz_intL;

    double dt_gxj_shiftj = gxxt_intL*shiftx_intL+gxyt_intL*shifty_intL*gxzt_intL*shiftz_intL;
    double dt_gyj_shiftj = gxyt_intL*shiftx_intL+gyyt_intL*shifty_intL*gyzt_intL*shiftz_intL;
    double dt_gzj_shiftj = gxzt_intL*shiftx_intL+gyzt_intL*shifty_intL*gzzt_intL*shiftz_intL;
    
    double dx_shiftj_gxj = gxx_intL*shiftxx_intL+gxy_intL*shiftyx_intL+gxz_intL*shiftzx_intL;
    double dx_shiftj_gyj = gxy_intL*shiftxx_intL+gyy_intL*shiftyx_intL+gyz_intL*shiftzx_intL;
    double dx_shiftj_gzj = gxz_intL*shiftxx_intL+gyz_intL*shiftyx_intL+gzz_intL*shiftzx_intL;
    
//     printf ("test2 gxx_intL: %10.6e\n",gxx_intL);
//     printf ("test2 gxy_intL: %10.6e\n",gxy_intL);
//     printf ("test2 gxz_intL: %10.6e\n",gxz_intL);
//     printf ("test2 shiftxy_intL: %10.6e\n",shiftxy_intL);
//     printf ("test2 shiftyy_intL: %10.6e\n",shiftyy_intL);
//     printf ("test2 shiftzy_intL: %10.6e\n",shiftzy_intL);
    double dy_shiftj_gxj = gxx_intL*shiftxy_intL+gxy_intL*shiftyy_intL+gxz_intL*shiftzy_intL;
//     printf ("test2 dy_shiftj_gxj: %10.6e\n",dy_shiftj_gxj);
    double dy_shiftj_gyj = gxy_intL*shiftxy_intL+gyy_intL*shiftyy_intL+gyz_intL*shiftzy_intL;
    double dy_shiftj_gzj = gxz_intL*shiftxy_intL+gyz_intL*shiftyy_intL+gzz_intL*shiftzy_intL;
    
    double dz_shiftj_gxj = gxx_intL*shiftxz_intL+gxy_intL*shiftyz_intL+gxz_intL*shiftzz_intL;
    double dz_shiftj_gyj = gxy_intL*shiftxz_intL+gyy_intL*shiftyz_intL+gyz_intL*shiftzz_intL;
    double dz_shiftj_gzj = gxz_intL*shiftxz_intL+gyz_intL*shiftyz_intL+gzz_intL*shiftzz_intL;
 
    double dt_shiftj_gxj = gxx_intL*shiftxt_intL+gxy_intL*shiftyt_intL+gxz_intL*shiftzt_intL;
    double dt_shiftj_gyj = gxy_intL*shiftxt_intL+gyy_intL*shiftyt_intL+gyz_intL*shiftzt_intL;
    double dt_shiftj_gzj = gxz_intL*shiftxt_intL+gyz_intL*shiftyt_intL+gzz_intL*shiftzt_intL;
    
    //double g4txx = Psi4*(4.0*shift_xL*phix_intL + dx_gxj_shiftj + dx_shiftj_gxj);
    //double g4tyx = Psi4*(4.0*shift_yL*phix_intL + dx_gyj_shiftj + dx_shiftj_gyj);
    //double g4tzx = Psi4*(4.0*shift_zL*phix_intL + dx_gzj_shiftj + dx_shiftj_gzj);
    double g4txx = g4xxx*shiftx_intL+g4xyx*shifty_intL+g4xzx*shiftz_intL
      +g4xx*shiftxx_intL+g4xy*shiftyx_intL+g4xz*shiftzx_intL;
    double g4tyx = g4xyx*shiftx_intL+g4yyx*shifty_intL+g4yzx*shiftz_intL
      +g4xy*shiftxx_intL+g4yy*shiftyx_intL+g4yz*shiftzx_intL;
    double g4tzx = g4xzx*shiftx_intL+g4yzx*shifty_intL+g4zzx*shiftz_intL
      +g4xz*shiftxx_intL+g4yz*shiftyx_intL+g4zz*shiftzx_intL;
    
    //double g4txy = Psi4*(4.0*shift_xL*phiy_intL + dy_gxj_shiftj + dy_shiftj_gxj);
    //double g4tyy = Psi4*(4.0*shift_yL*phiy_intL + dy_gyj_shiftj + dy_shiftj_gyj);
    //double g4tzy = Psi4*(4.0*shift_zL*phiy_intL + dy_gzj_shiftj + dy_shiftj_gzj);
    double g4txy = g4xxy*shiftx_intL+g4xyy*shifty_intL+g4xzy*shiftz_intL
      +g4xx*shiftxy_intL+g4xy*shiftyy_intL+g4xz*shiftzy_intL;
    double g4tyy = g4xyy*shiftx_intL+g4yyy*shifty_intL+g4yzy*shiftz_intL
      +g4xy*shiftxy_intL+g4yy*shiftyy_intL+g4yz*shiftzy_intL;
    double g4tzy = g4xzy*shiftx_intL+g4yzy*shifty_intL+g4zzy*shiftz_intL
      +g4xz*shiftxy_intL+g4yz*shiftyy_intL+g4zz*shiftzy_intL;
  
    //double g4txz = Psi4*(4.0*shift_xL*phiz_intL + dz_gxj_shiftj + dz_shiftj_gxj);
    //double g4tyz = Psi4*(4.0*shift_yL*phiz_intL + dz_gyj_shiftj + dz_shiftj_gyj);
    //double g4tzz = Psi4*(4.0*shift_zL*phiz_intL + dz_gzj_shiftj + dz_shiftj_gzj);
    double g4txz = g4xxz*shiftx_intL+g4xyz*shifty_intL+g4xzz*shiftz_intL
      +g4xx*shiftxz_intL+g4xy*shiftyz_intL+g4xz*shiftzz_intL;
    double g4tyz = g4xyz*shiftx_intL+g4yyz*shifty_intL+g4yzz*shiftz_intL
      +g4xy*shiftxz_intL+g4yy*shiftyz_intL+g4yz*shiftzz_intL;
    double g4tzz = g4xzz*shiftx_intL+g4yzz*shifty_intL+g4zzz*shiftz_intL
      +g4xz*shiftxz_intL+g4yz*shiftyz_intL+g4zz*shiftzz_intL;
 
    double g4txt = Psi4*(4.0*shift_xL*phit_intL + dt_gxj_shiftj + dt_shiftj_gxj);
    double g4tyt = Psi4*(4.0*shift_yL*phit_intL + dt_gyj_shiftj + dt_shiftj_gyj);
    double g4tzt = Psi4*(4.0*shift_zL*phit_intL + dt_gzj_shiftj + dt_shiftj_gzj);
   
    double g4ttx = -2.0*lapse_intL*lapsex_intL +  
      g4xxx*shiftx_intL*shiftx_intL+g4yyx*shifty_intL*shifty_intL+g4zzx*shiftz_intL*shiftz_intL+
      2.0*(g4xyx*shiftx_intL*shifty_intL+g4xzx*shiftx_intL*shiftz_intL+g4yzx*shifty_intL*shiftz_intL)+
      2.0*Psi4*(shift_xL*shiftxx_intL+shift_yL*shiftyx_intL+shift_zL*shiftzx_intL);
    double g4tty = -2.0*lapse_intL*lapsey_intL +  
      g4xxy*shiftx_intL*shiftx_intL+g4yyy*shifty_intL*shifty_intL+g4zzy*shiftz_intL*shiftz_intL+
      2.0*(g4xyy*shiftx_intL*shifty_intL+g4xzy*shiftx_intL*shiftz_intL+g4yzy*shifty_intL*shiftz_intL)+
      2.0*Psi4*(shift_xL*shiftxy_intL+shift_yL*shiftyy_intL+shift_zL*shiftzy_intL);
    double g4ttz = -2.0*lapse_intL*lapsez_intL +  
      g4xxz*shiftx_intL*shiftx_intL+g4yyz*shifty_intL*shifty_intL+g4zzz*shiftz_intL*shiftz_intL+
      2.0*(g4xyz*shiftx_intL*shifty_intL+g4xzz*shiftx_intL*shiftz_intL+g4yzz*shifty_intL*shiftz_intL)+
      2.0*Psi4*(shift_xL*shiftxz_intL+shift_yL*shiftyz_intL+shift_zL*shiftzz_intL);
    double g4ttt = -2.0*lapse_intL*lapset_intL +  
      g4xxt*shiftx_intL*shiftx_intL+g4yyt*shifty_intL*shifty_intL+g4zzt*shiftz_intL*shiftz_intL+
      2.0*(g4xyt*shiftx_intL*shifty_intL+g4xzt*shiftx_intL*shiftz_intL+g4yzt*shifty_intL*shiftz_intL)+
      2.0*Psi4*(shift_xL*shiftxt_intL+shift_yL*shiftyt_intL+shift_zL*shiftzt_intL);
    
//   double g4ttx = -2.0*lapse_intL*lapsex_intL +
//     Psi4*(4.0*(shift_xL*shiftx_intL+shift_yL*shifty_intL+shift_zL*shiftz_intL)*phix_intL +
//     dx_gxj_shiftj*shiftx_intL+dx_gyj_shiftj*shifty_intL+dx_gzj_shiftj*shiftz_intL +
//     2.0 * (dx_shiftj_gxj*shiftx_intL+dx_shiftj_gyj*shifty_intL+dx_shiftj_gzj*shiftz_intL));
    
//     double g4tty = -2.0*lapse_intL*lapsey_intL +
//       Psi4*(4.0*(shift_xL*shiftx_intL+shift_yL*shifty_intL+shift_zL*shiftz_intL)*phiy_intL +
// 	    dy_gxj_shiftj*shiftx_intL+dy_gyj_shiftj*shifty_intL+dy_gzj_shiftj*shiftz_intL +
// 	    2.0 * (dy_shiftj_gxj*shiftx_intL+dy_shiftj_gyj*shifty_intL+dy_shiftj_gzj*shiftz_intL));
//     double g4ttz = -2.0*lapse_intL*lapsez_intL +
//       Psi4*(4.0*(shift_xL*shiftx_intL+shift_yL*shifty_intL+shift_zL*shiftz_intL)*phiz_intL +
// 	    dz_gxj_shiftj*shiftx_intL+dz_gyj_shiftj*shifty_intL+dz_gzj_shiftj*shiftz_intL +
// 	    2.0 * (dz_shiftj_gxj*shiftx_intL+dz_shiftj_gyj*shifty_intL+dz_shiftj_gzj*shiftz_intL));
//     double g4ttt = -2.0*lapse_intL*lapset_intL +
//       Psi4*(4.0*(shift_xL*shiftx_intL+shift_yL*shifty_intL+shift_zL*shiftz_intL)*phit_intL +
// 	    dt_gxj_shiftj*shiftx_intL+dt_gyj_shiftj*shifty_intL+dt_gzj_shiftj*shiftz_intL +
// 	    2.0 * (dt_shiftj_gxj*shiftx_intL+dt_shiftj_gyj*shifty_intL+dt_shiftj_gzj*shiftz_intL));
    
    double Gam_ttt = 0.5*g4ttt;
    double Gam_ttx = 0.5*g4ttx;
    double Gam_tty = 0.5*g4tty;
    double Gam_ttz = 0.5*g4ttz;
    double Gam_txx = 0.5*(g4txx+g4txx-g4xxt);
    double Gam_txy = 0.5*(g4txy+g4tyx-g4xyt);
    double Gam_txz = 0.5*(g4txz+g4tzx-g4xzt);
    double Gam_tyy = 0.5*(g4tyy+g4tyy-g4yyt);
    double Gam_tyz = 0.5*(g4tyz+g4tzy-g4yzt);
    double Gam_tzz = 0.5*(g4tzz+g4tzz-g4zzt);

    double Gam_xtt = 0.5*(g4txt+g4txt-g4ttx);
    double Gam_xtx = 0.5*g4xxt;
    double Gam_xty = 0.5*(g4txy+g4xyt-g4tyx);
    double Gam_xtz = 0.5*(g4txz+g4xzt-g4tzx);
    double Gam_xxx = 0.5*g4xxx;
    double Gam_xxy = 0.5*g4xxy;
    double Gam_xxz = 0.5*g4xxz;
    double Gam_xyy = 0.5*(g4xyy+g4xyy-g4yyx);
    double Gam_xyz = 0.5*(g4xyz+g4xzy-g4yzx);
    double Gam_xzz = 0.5*(g4xzz+g4xzz-g4zzx);

    double Gam_ytt = 0.5*(g4tyt+g4tyt-g4tty);
    double Gam_ytx = 0.5*(g4tyx+g4xyt-g4txy);
    double Gam_yty = 0.5*g4yyt;
    double Gam_ytz = 0.5*(g4tyz+g4yzt-g4tzy);
    double Gam_yxx = 0.5*(g4xyx+g4xyx-g4xxy);
    double Gam_yxy = 0.5*g4yyx;
    double Gam_yxz = 0.5*(g4xyz+g4yzx-g4xzy);
    double Gam_yyy = 0.5*g4yyy;
    double Gam_yyz = 0.5*g4yyz;
    double Gam_yzz = 0.5*(g4yzz+g4yzz-g4zzy);

    double Gam_ztt = 0.5*(g4tzt+g4tzt-g4ttz);
    double Gam_ztx = 0.5*(g4tzx+g4xzt-g4txz);
    double Gam_zty = 0.5*(g4tzy+g4yzt-g4tyz);
    double Gam_ztz = 0.5*g4zzt;
    double Gam_zxx = 0.5*(g4xzx+g4xzx-g4xxz);
    double Gam_zxy = 0.5*(g4xzy+g4yzx-g4xyz);
    double Gam_zxz = 0.5*g4zzx;
    double Gam_zyy = 0.5*(g4yzy+g4yzy-g4yyz);
    double Gam_zyz = 0.5*g4zzy;
    double Gam_zzz = 0.5*g4zzz;

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

    double Gamt_tt = g4uptt*Gam_ttt+g4uptx*Gam_xtt+g4upty*Gam_ytt+g4uptz*Gam_ztt;
    double Gamt_tx = g4uptt*Gam_ttx+g4uptx*Gam_xtx+g4upty*Gam_ytx+g4uptz*Gam_ztx;
    double Gamt_ty = g4uptt*Gam_tty+g4uptx*Gam_xty+g4upty*Gam_yty+g4uptz*Gam_zty;
    double Gamt_tz = g4uptt*Gam_ttz+g4uptx*Gam_xtz+g4uptz*Gam_ytz+g4uptz*Gam_ztz;
    double Gamt_xx = g4uptt*Gam_txx+g4uptx*Gam_xxx+g4upty*Gam_yxx+g4uptz*Gam_zxx;
    double Gamt_xy = g4uptt*Gam_txy+g4uptx*Gam_xxy+g4upty*Gam_yxy+g4uptz*Gam_zxy;
    double Gamt_xz = g4uptt*Gam_txz+g4uptx*Gam_xxz+g4upty*Gam_yxz+g4uptz*Gam_zxz;
    double Gamt_yy = g4uptt*Gam_tyy+g4uptx*Gam_xyy+g4upty*Gam_yyy+g4uptz*Gam_zyy;
    double Gamt_yz = g4uptt*Gam_tyz+g4uptx*Gam_xyz+g4upty*Gam_yyz+g4uptz*Gam_zyz;
    double Gamt_zz = g4uptt*Gam_tzz+g4uptx*Gam_xzz+g4upty*Gam_yzz+g4uptz*Gam_zzz;

    double Gamx_tt = g4uptx*Gam_ttt+g4upxx*Gam_xtt+g4upxy*Gam_ytt+g4upxz*Gam_ztt;
    double Gamx_tx = g4uptx*Gam_ttx+g4upxx*Gam_xtx+g4upxy*Gam_ytx+g4upxz*Gam_ztx;
    double Gamx_ty = g4uptx*Gam_tty+g4upxx*Gam_xty+g4upxy*Gam_yty+g4upxz*Gam_zty;
    double Gamx_tz = g4uptx*Gam_ttz+g4upxx*Gam_xtz+g4upxy*Gam_ytz+g4upxz*Gam_ztz;
    double Gamx_xx = g4uptx*Gam_txx+g4upxx*Gam_xxx+g4upxy*Gam_yxx+g4upxz*Gam_zxx;
    double Gamx_xy = g4uptx*Gam_txy+g4upxx*Gam_xxy+g4upxy*Gam_yxy+g4upxz*Gam_zxy;
    double Gamx_xz = g4uptx*Gam_txz+g4upxx*Gam_xxz+g4upxy*Gam_yxz+g4upxz*Gam_zxz;
    double Gamx_yy = g4uptx*Gam_tyy+g4upxx*Gam_xyy+g4upxy*Gam_yyy+g4upxz*Gam_zyy;
    double Gamx_yz = g4uptx*Gam_tyz+g4upxx*Gam_xyz+g4upxy*Gam_yyz+g4upxz*Gam_zyz;
    double Gamx_zz = g4uptx*Gam_tzz+g4upxx*Gam_xzz+g4upxy*Gam_yzz+g4upxz*Gam_zzz;
  
    double Gamy_tt = g4upty*Gam_ttt+g4upxy*Gam_xtt+g4upyy*Gam_ytt+g4upyz*Gam_ztt;
    double Gamy_tx = g4upty*Gam_ttx+g4upxy*Gam_xtx+g4upyy*Gam_ytx+g4upyz*Gam_ztx;
    double Gamy_ty = g4upty*Gam_tty+g4upxy*Gam_xty+g4upyy*Gam_yty+g4upyz*Gam_zty;
    double Gamy_tz = g4upty*Gam_ttz+g4upxy*Gam_xtz+g4upyz*Gam_ytz+g4upyz*Gam_ztz;
    double Gamy_xx = g4upty*Gam_txx+g4upxy*Gam_xxx+g4upyy*Gam_yxx+g4upyz*Gam_zxx;
    double Gamy_xy = g4upty*Gam_txy+g4upxy*Gam_xxy+g4upyy*Gam_yxy+g4upyz*Gam_zxy;
    double Gamy_xz = g4upty*Gam_txz+g4upxy*Gam_xxz+g4upyy*Gam_yxz+g4upyz*Gam_zxz;
    double Gamy_yy = g4upty*Gam_tyy+g4upxy*Gam_xyy+g4upyy*Gam_yyy+g4upyz*Gam_zyy;
    double Gamy_yz = g4upty*Gam_tyz+g4upxy*Gam_xyz+g4upyy*Gam_yyz+g4upyz*Gam_zyz;
    double Gamy_zz = g4upty*Gam_tzz+g4upxy*Gam_xzz+g4upyy*Gam_yzz+g4upyz*Gam_zzz;

    double Gamz_tt = g4uptz*Gam_ttt+g4upxz*Gam_xtt+g4upyz*Gam_ytt+g4upzz*Gam_ztt;
    double Gamz_tx = g4uptz*Gam_ttx+g4upxz*Gam_xtx+g4upyz*Gam_ytx+g4upzz*Gam_ztx;
    double Gamz_ty = g4uptz*Gam_tty+g4upxz*Gam_xty+g4upyz*Gam_yty+g4upzz*Gam_zty;
    double Gamz_tz = g4uptz*Gam_ttz+g4upxz*Gam_xtz+g4upyz*Gam_ytz+g4upzz*Gam_ztz;
    double Gamz_xx = g4uptz*Gam_txx+g4upxz*Gam_xxx+g4upyz*Gam_yxx+g4upzz*Gam_zxx;
    double Gamz_xy = g4uptz*Gam_txy+g4upxz*Gam_xxy+g4upyz*Gam_yxy+g4upzz*Gam_zxy;
    double Gamz_xz = g4uptz*Gam_txz+g4upxz*Gam_xxz+g4upyz*Gam_yxz+g4upzz*Gam_zxz;
    double Gamz_yy = g4uptz*Gam_tyy+g4upxz*Gam_xyy+g4upyz*Gam_yyy+g4upzz*Gam_zyy;
    double Gamz_yz = g4uptz*Gam_tyz+g4upxz*Gam_xyz+g4upyz*Gam_yyz+g4upzz*Gam_zyz;
    double Gamz_zz = g4uptz*Gam_tzz+g4upxz*Gam_xzz+g4upyz*Gam_yzz+g4upzz*Gam_zzz;
        
    //    printf("output_radL: %10.6f\n",output_radL);

    if ((i==23)&&(j==23)&&(k==10)){
      double XL=output_radL*sin(thetaL)*cos(phiL)+BH_posx;
      double YL=output_radL*sin(thetaL)*sin(phiL)+BH_posy; 
      double ZL=output_radL*cos(thetaL);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tt: %10.6e\n",XL,YL,ZL,g4tt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tx: %10.6e\n",XL,YL,ZL,g4tx); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4ty: %10.6e\n",XL,YL,ZL,g4ty); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tz: %10.6e\n",XL,YL,ZL,g4tz); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xx: %10.6e\n",XL,YL,ZL,g4xx); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xy: %10.6e\n",XL,YL,ZL,g4xy); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xz: %10.6e\n",XL,YL,ZL,g4xz); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yy: %10.6e\n",XL,YL,ZL,g4yy); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yz: %10.6e\n",XL,YL,ZL,g4yz); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4zz: %10.6e\n",XL,YL,ZL,g4zz); 
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4ttt: %10.6e\n",XL,YL,ZL,g4ttt);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4txt: %10.6e\n",XL,YL,ZL,g4txt);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tyt: %10.6e\n",XL,YL,ZL,g4tyt);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tzt: %10.6e\n",XL,YL,ZL,g4tzt);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4ttx: %10.6e\n",XL,YL,ZL,g4ttx);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4txx: %10.6e\n",XL,YL,ZL,g4txx);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tyx: %10.6e\n",XL,YL,ZL,g4tyx);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tzx: %10.6e\n",XL,YL,ZL,g4tzx);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tty: %10.6e\n",XL,YL,ZL,g4tty);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4txy: %10.6e\n",XL,YL,ZL,g4txy);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tyy: %10.6e\n",XL,YL,ZL,g4tyy);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tzy: %10.6e\n",XL,YL,ZL,g4tzy);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4ttz: %10.6e\n",XL,YL,ZL,g4ttz);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4txz: %10.6e\n",XL,YL,ZL,g4txz);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tyz: %10.6e\n",XL,YL,ZL,g4tyz);  
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4tzz: %10.6e\n",XL,YL,ZL,g4tzz);  
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, lapsex_intL: %10.6e\n",XL,YL,ZL,lapsex_intL);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, shiftxx_intL: %10.6e\n",XL,YL,ZL,shiftxx_intL);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, shiftyx_intL: %10.6e\n",XL,YL,ZL,shiftyx_intL);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, shiftzx_intL: %10.6e\n",XL,YL,ZL,shiftzx_intL);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xxx: %10.6e\n",XL,YL,ZL,g4xxx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xyx: %10.6e\n",XL,YL,ZL,g4xyx); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xzx: %10.6e\n",XL,YL,ZL,g4xzx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yyx: %10.6e\n",XL,YL,ZL,g4yyx); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yzx: %10.6e\n",XL,YL,ZL,g4yzx); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4zzx: %10.6e\n",XL,YL,ZL,g4zzx); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xxy: %10.6e\n",XL,YL,ZL,g4xxy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xyy: %10.6e\n",XL,YL,ZL,g4xyy); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xzy: %10.6e\n",XL,YL,ZL,g4xzy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yyy: %10.6e\n",XL,YL,ZL,g4yyy); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yzy: %10.6e\n",XL,YL,ZL,g4yzy); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4zzy: %10.6e\n",XL,YL,ZL,g4zzy); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xxz: %10.6e\n",XL,YL,ZL,g4xxz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xyz: %10.6e\n",XL,YL,ZL,g4xyz); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4xzz: %10.6e\n",XL,YL,ZL,g4xzz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yyz: %10.6e\n",XL,YL,ZL,g4yyz); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4yzz: %10.6e\n",XL,YL,ZL,g4yzz); 
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, g4zzz: %10.6e\n",XL,YL,ZL,g4zzz); 
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ttt: %10.6e\n",XL,YL,ZL,Gam_ttt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ttx: %10.6e\n",XL,YL,ZL,Gam_ttx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_tty: %10.6e\n",XL,YL,ZL,Gam_tty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ttz: %10.6e\n",XL,YL,ZL,Gam_ttz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_txx: %10.6e\n",XL,YL,ZL,Gam_txx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_txy: %10.6e\n",XL,YL,ZL,Gam_txy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_txz: %10.6e\n",XL,YL,ZL,Gam_txz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_tyy: %10.6e\n",XL,YL,ZL,Gam_tyy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_tyz: %10.6e\n",XL,YL,ZL,Gam_tyz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_tzz: %10.6e\n",XL,YL,ZL,Gam_tzz);
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xtt: %10.6e\n",XL,YL,ZL,Gam_xtt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xtx: %10.6e\n",XL,YL,ZL,Gam_xtx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xty: %10.6e\n",XL,YL,ZL,Gam_xty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xtz: %10.6e\n",XL,YL,ZL,Gam_xtz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xxx: %10.6e\n",XL,YL,ZL,Gam_xxx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xxy: %10.6e\n",XL,YL,ZL,Gam_xxy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xxz: %10.6e\n",XL,YL,ZL,Gam_xxz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xyy: %10.6e\n",XL,YL,ZL,Gam_xyy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xyz: %10.6e\n",XL,YL,ZL,Gam_xyz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_xzz: %10.6e\n",XL,YL,ZL,Gam_xzz);
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ytt: %10.6e\n",XL,YL,ZL,Gam_ytt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ytx: %10.6e\n",XL,YL,ZL,Gam_ytx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yty: %10.6e\n",XL,YL,ZL,Gam_yty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ytz: %10.6e\n",XL,YL,ZL,Gam_ytz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yxx: %10.6e\n",XL,YL,ZL,Gam_yxx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yxy: %10.6e\n",XL,YL,ZL,Gam_yxy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yxz: %10.6e\n",XL,YL,ZL,Gam_yxz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yyy: %10.6e\n",XL,YL,ZL,Gam_yyy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yyz: %10.6e\n",XL,YL,ZL,Gam_yyz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_yzz: %10.6e\n",XL,YL,ZL,Gam_yzz);
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ztt: %10.6e\n",XL,YL,ZL,Gam_ztt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ztx: %10.6e\n",XL,YL,ZL,Gam_ztx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zty: %10.6e\n",XL,YL,ZL,Gam_zty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_ztz: %10.6e\n",XL,YL,ZL,Gam_ztz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zxx: %10.6e\n",XL,YL,ZL,Gam_zxx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zxy: %10.6e\n",XL,YL,ZL,Gam_zxy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zxz: %10.6e\n",XL,YL,ZL,Gam_zxz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zyy: %10.6e\n",XL,YL,ZL,Gam_zyy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zyz: %10.6e\n",XL,YL,ZL,Gam_zyz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gam_zzz: %10.6e\n",XL,YL,ZL,Gam_zzz);  
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_tt: %10.6e\n",XL,YL,ZL,Gamt_tt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_tx: %10.6e\n",XL,YL,ZL,Gamt_tx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_ty: %10.6e\n",XL,YL,ZL,Gamt_ty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_tz: %10.6e\n",XL,YL,ZL,Gamt_tz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_xx: %10.6e\n",XL,YL,ZL,Gamt_xx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_xy: %10.6e\n",XL,YL,ZL,Gamt_xy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_xz: %10.6e\n",XL,YL,ZL,Gamt_xz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_yy: %10.6e\n",XL,YL,ZL,Gamt_yy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_yz: %10.6e\n",XL,YL,ZL,Gamt_yz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamt_zz: %10.6e\n",XL,YL,ZL,Gamt_zz);
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_tt: %10.6e\n",XL,YL,ZL,Gamx_tt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_tx: %10.6e\n",XL,YL,ZL,Gamx_tx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_ty: %10.6e\n",XL,YL,ZL,Gamx_ty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_tz: %10.6e\n",XL,YL,ZL,Gamx_tz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_xx: %10.6e\n",XL,YL,ZL,Gamx_xx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_xy: %10.6e\n",XL,YL,ZL,Gamx_xy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_xz: %10.6e\n",XL,YL,ZL,Gamx_xz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_yy: %10.6e\n",XL,YL,ZL,Gamx_yy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_yz: %10.6e\n",XL,YL,ZL,Gamx_yz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamx_zz: %10.6e\n",XL,YL,ZL,Gamx_zz);
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_tt: %10.6e\n",XL,YL,ZL,Gamy_tt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_tx: %10.6e\n",XL,YL,ZL,Gamy_tx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_ty: %10.6e\n",XL,YL,ZL,Gamy_ty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_tz: %10.6e\n",XL,YL,ZL,Gamy_tz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_xx: %10.6e\n",XL,YL,ZL,Gamy_xx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_xy: %10.6e\n",XL,YL,ZL,Gamy_xy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_xz: %10.6e\n",XL,YL,ZL,Gamy_xz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_yy: %10.6e\n",XL,YL,ZL,Gamy_yy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_yz: %10.6e\n",XL,YL,ZL,Gamy_yz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamy_zz: %10.6e\n",XL,YL,ZL,Gamy_zz);
      printf("\n");
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_tt: %10.6e\n",XL,YL,ZL,Gamz_tt);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_tx: %10.6e\n",XL,YL,ZL,Gamz_tx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_ty: %10.6e\n",XL,YL,ZL,Gamz_ty);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_tz: %10.6e\n",XL,YL,ZL,Gamz_tz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_xx: %10.6e\n",XL,YL,ZL,Gamz_xx);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_xy: %10.6e\n",XL,YL,ZL,Gamz_xy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_xz: %10.6e\n",XL,YL,ZL,Gamz_xz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_yy: %10.6e\n",XL,YL,ZL,Gamz_yy);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_yz: %10.6e\n",XL,YL,ZL,Gamz_yz);
      printf("X: %10.6e, Y: %10.6e, Z: %10.6e, Gamz_zz: %10.6e\n",XL,YL,ZL,Gamz_zz);  
    }
    //   double output_rad2L=exp(output_lograd2L);
    //   double output_rad3L=exp(output_lograd3L);
    
    
    //try getting rid of rad,phi,theta stuff.  It can be reconstructed later
    //outfile.write((char *) &output_rad1L,sizeof(double));
    //outfile.write((char *) &output_rad2L,sizeof(double));
    //outfile.write((char *) &output_rad3L,sizeof(double));
    //outfile.write((char *) &phiL,sizeof(double));
    //outfile.write((char *) &thetaL,sizeof(double));
    
//     printf ("test2 output_radL: %10.6e\n",output_radL);
//     printf ("test2 g4uptx: %10.6e\n",g4uptx);
//     printf ("test2 g4upxx: %10.6e\n",g4upxx);
//     printf ("test2 g4upxy: %10.6e\n",g4upxy);
//     printf ("test2 g4upxz: %10.6e\n",g4upxz);
//     printf ("test2 Psi4: %10.6e\n",Psi4);
//     printf ("test2 shiftx_intL: %10.6e\n",shiftx_intL);
//     printf ("test2 shifty_intL: %10.6e\n",shifty_intL);
//     printf ("test2 shiftz_intL: %10.6e\n",shiftz_intL);
//     printf ("test2 shift_xL: %10.6e\n",shift_xL); 
//     printf ("test2 shift_yL: %10.6e\n",shift_yL); 
//     printf ("test2 shift_zL: %10.6e\n",shift_zL); 
//     printf ("test2 phiy_intL: %10.6e\n",phiy_intL);
//     printf ("test2 gxxy_intL: %10.6e\n",gxxy_intL);
//     printf ("test2 gxyy_intL: %10.6e\n",gxyy_intL);
//     printf ("test2 gxzy_intL: %10.6e\n",gxzy_intL);
//     printf ("test2 dy_gxj_shiftj: %10.6e\n",dy_gxj_shiftj);
//     printf ("test2 dy_gyj_shiftj: %10.6e\n",dy_gyj_shiftj);
//     printf ("test2 dy_gzj_shiftj: %10.6e\n",dy_gzj_shiftj);
//     printf ("test2 dy_shiftj_gxj: %10.6e\n",dy_shiftj_gxj);
//     printf ("test2 dy_shiftj_gyj: %10.6e\n",dy_shiftj_gyj); 
//     printf ("test2 dy_shiftj_gzj: %10.6e\n",dy_shiftj_gzj); 

//     printf ("test2 g4tty: %10.6e\n",g4tty);
//     printf ("test2 Gam_tty: %10.6e\n",Gam_tty);
//     printf ("test2 Gam_xty: %10.6e\n",Gam_xty);
//     printf ("test2 Gam_yty: %10.6e\n",Gam_yty);
//     printf ("test2 Gam_zty: %10.6e\n",Gam_zty);
//     printf ("test2 Gamx_ty: %10.6e\n",Gamx_ty); 
    
    double uxL = vx_int[n]*u0_int[n];
    double uyL = vy_int[n]*u0_int[n];
    double uzL = vz_int[n]*u0_int[n];

    outfile.write((char *) &rho_b_int[n],sizeof(double));
    outfile.write((char *) &P_int[n],sizeof(double));
    outfile.write((char *) &uxL,sizeof(double));
    outfile.write((char *) &uyL,sizeof(double));
    outfile.write((char *) &uzL,sizeof(double));
    outfile.write((char *) &u0_int[n],sizeof(double));
    outfile.write((char *) &g4tt,sizeof(double));
    outfile.write((char *) &g4tx,sizeof(double));
    outfile.write((char *) &g4ty,sizeof(double));
    outfile.write((char *) &g4tz,sizeof(double));
    outfile.write((char *) &g4xx,sizeof(double));
    outfile.write((char *) &g4xy,sizeof(double));
    outfile.write((char *) &g4xz,sizeof(double));
    outfile.write((char *) &g4yy,sizeof(double));
    outfile.write((char *) &g4yz,sizeof(double));
    outfile.write((char *) &g4zz,sizeof(double));
    outfile.write((char *) &Gamt_tt,sizeof(double));
    outfile.write((char *) &Gamt_tx,sizeof(double));
    outfile.write((char *) &Gamt_ty,sizeof(double));
    outfile.write((char *) &Gamt_tz,sizeof(double));
    outfile.write((char *) &Gamt_xx,sizeof(double));
    outfile.write((char *) &Gamt_xy,sizeof(double));
    outfile.write((char *) &Gamt_xz,sizeof(double));
    outfile.write((char *) &Gamt_yy,sizeof(double));
    outfile.write((char *) &Gamt_yz,sizeof(double));
    outfile.write((char *) &Gamt_zz,sizeof(double));
    outfile.write((char *) &Gamx_tt,sizeof(double));
    outfile.write((char *) &Gamx_tx,sizeof(double));
    outfile.write((char *) &Gamx_ty,sizeof(double));
    outfile.write((char *) &Gamx_tz,sizeof(double));
    outfile.write((char *) &Gamx_xx,sizeof(double));
    outfile.write((char *) &Gamx_xy,sizeof(double));
    outfile.write((char *) &Gamx_xz,sizeof(double));
    outfile.write((char *) &Gamx_yy,sizeof(double));
    outfile.write((char *) &Gamx_yz,sizeof(double));
    outfile.write((char *) &Gamx_zz,sizeof(double));
    outfile.write((char *) &Gamy_tt,sizeof(double));
    outfile.write((char *) &Gamy_tx,sizeof(double));
    outfile.write((char *) &Gamy_ty,sizeof(double));
    outfile.write((char *) &Gamy_tz,sizeof(double));
    outfile.write((char *) &Gamy_xx,sizeof(double));
    outfile.write((char *) &Gamy_xy,sizeof(double));
    outfile.write((char *) &Gamy_xz,sizeof(double));
    outfile.write((char *) &Gamy_yy,sizeof(double));
    outfile.write((char *) &Gamy_yz,sizeof(double));
    outfile.write((char *) &Gamy_zz,sizeof(double));
    outfile.write((char *) &Gamz_tt,sizeof(double));
    outfile.write((char *) &Gamz_tx,sizeof(double));
    outfile.write((char *) &Gamz_ty,sizeof(double));
    outfile.write((char *) &Gamz_tz,sizeof(double));
    outfile.write((char *) &Gamz_xx,sizeof(double));
    outfile.write((char *) &Gamz_xy,sizeof(double));
    outfile.write((char *) &Gamz_xz,sizeof(double));
    outfile.write((char *) &Gamz_yy,sizeof(double));
    outfile.write((char *) &Gamz_yz,sizeof(double));
    outfile.write((char *) &Gamz_zz,sizeof(double));

    n = n + 1;
    // printf ("test2 i: %d, j: %d, k: %d\n",i,j,k);
  }
  //  printf ("hello 3\n");

  //close files
  outfile.close();
}
extern "C" void CCTK_FCALL CCTK_FNAME(bbh_bondi_IO)
  (const cGH **cctkGH, double &time,int &iteration,int &gridnum,
   double &output_radmin,double &output_radmax,
   int &output_Nlograd,int &output_Nphi,int &output_Ntheta,
   double *rho_b_int,double *P_int,double *vx_int,double *vy_int,double *vz_int,double *u0_int,
   double *lapse_int,double *shiftx_int,double *shifty_int,double *shiftz_int,
   double *lapsex_int,double *shiftxx_int,double *shiftyx_int,double *shiftzx_int,
   double *lapsey_int,double *shiftxy_int,double *shiftyy_int,double *shiftzy_int,
   double *lapsez_int,double *shiftxz_int,double *shiftyz_int,double *shiftzz_int,
   double *lapset_int,double *shiftxt_int,double *shiftyt_int,double *shiftzt_int,
   double *phi_int,double *gxx_int,double *gxy_int,double *gxz_int,
   double *gyy_int,double *gyz_int, double *gzz_int,
   double *gupxx_int,double *gupxy_int,double *gupxz_int,
   double *gupyy_int,double *gupyz_int, double *gupzz_int,
   double *phix_int,double *gxxx_int,double *gxyx_int,double *gxzx_int,
   double *gyyx_int,double *gyzx_int, double *gzzx_int,
   double *phiy_int,double *gxxy_int,double *gxyy_int,double *gxzy_int,
   double *gyyy_int,double *gyzy_int, double *gzzy_int,
   double *phiz_int,double *gxxz_int,double *gxyz_int,double *gxzz_int,
   double *gyyz_int,double *gyzz_int, double *gzzz_int,
   double *phit_int,double *gxxt_int,double *gxyt_int,double *gxzt_int,
   double *gyyt_int,double *gyzt_int,double *gzzt_int,
   double &BH_posx,double &BH_posy)
{  
  bbh_bondi_IO(*cctkGH, time, iteration,gridnum,
	       output_radmin,output_radmax,
	       output_Nlograd,output_Nphi,output_Ntheta,
	       rho_b_int,P_int,vx_int,vy_int,vz_int,u0_int,
	       lapse_int,shiftx_int,shifty_int,shiftz_int,
	       lapsex_int,shiftxx_int,shiftyx_int,shiftzx_int,
	       lapsey_int,shiftxy_int,shiftyy_int,shiftzy_int,
	       lapsez_int,shiftxz_int,shiftyz_int,shiftzz_int,
	       lapset_int,shiftxt_int,shiftyt_int,shiftzt_int,
	       phi_int,gxx_int,gxy_int,gxz_int,
	       gyy_int,gyz_int, gzz_int,
	       gupxx_int,gupxy_int,gupxz_int,
	       gupyy_int,gupyz_int, gupzz_int,
	       phix_int,gxxx_int,gxyx_int,gxzx_int,
	       gyyx_int,gyzx_int, gzzx_int,
	       phiy_int,gxxy_int,gxyy_int,gxzy_int,
	       gyyy_int,gyzy_int, gzzy_int,
	       phiz_int,gxxz_int,gxyz_int,gxzz_int,
	       gyyz_int,gyzz_int, gzzz_int,
	       phit_int,gxxt_int,gxyt_int,gxzt_int,
	       gyyt_int,gyzt_int,gzzt_int,
	       BH_posx,BH_posy);
}


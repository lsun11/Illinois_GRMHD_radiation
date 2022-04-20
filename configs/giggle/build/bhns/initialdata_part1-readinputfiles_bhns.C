//----------------------------------------------------------------------------
//
// $Id: $
//
//----------------------------------------------------------------------------
//
// Contains methods for integrating the OV equations
//
//----------------------------------------------------------------------------

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
//#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "nrutil.h"

#include "Symmetry.h"

//#define PI acos(-1.0) // <-- BAD.  Use M_PI.
#define N_INT 1000

void Setup_Hybrid_EOS_magnetar_bhns(double K_poly, double n, char * eos_file,int neos,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, int ergo_star);
void polar_derivs_bhns(double **f, double **f1, double **f2,
		  double **f12,
		  int NS, int NU, double RADEQUAT);
double bicub_int_bhns(double **f, int NS, int NU,
		 double RADEQUAT,
		 double riso, double theta,
		 double **f1, double **f2, double **f12);
void bcuint_bhns(double *y, double *y1, double *y2, double *y12, double x1l,
	    double x1u, double x2l, double x2u, double x1, double x2,
	    double & ansy, double & ansy1, double & ansy2);
void bcucof_bhns(double *y, double *y1, double *y2, double *y12,
	    double d1, double d2, double **c);
double findmin_bhns(double *array, int &size);

static char *rcsid = "$Meow...$";
CCTK_FILEVERSION(Magnetar_initialdata_readfiles_bhns)


extern "C" void Magnetar_initialdata_readfiles_bhns(CCTK_ARGUMENTS) 
{
  using namespace std;
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  printf("Magnetar_initialdata_readfiles_bhns");

  if(genID_cmdline_output_enable==1) {

    printf("USE Magnetar_initialdata_readfiles_bhns to generate regrid information!!!\n");
    double dx = CCTK_DELTA_SPACE(0);
    double dy = CCTK_DELTA_SPACE(1);
    double dz = CCTK_DELTA_SPACE(2);
    //    double xmin = findmin_bhns(x,cctk_lsh[0]); //x[CCTK_GFINDEX3D(cctkGH,0,0,0)];
    //    double ymin = findmin_bhns(y,cctk_lsh[1]); //y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
    //    double zmin = findmin_bhns(z,cctk_lsh[2]); //z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
    double xmin = x[CCTK_GFINDEX3D(cctkGH,0,0,0)];
    double ymin = y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
    double zmin = z[CCTK_GFINDEX3D(cctkGH,0,0,0)];
    double xmax= xmin + cctk_lsh[0]*dx;
    double ymax= ymin + cctk_lsh[1]*dy;
    double zmax= zmin + cctk_lsh[2]*dz;    
    static int nStatic;
    nStatic += 1;

    int filename_checksum = cctk_lsh[0]*301+cctk_lsh[1]*10+cctk_lsh[2]*1000000  + 1241*CCTK_MyProc(cctkGH) + (int)(999999.*dx) + (int)(100100100.*dy) + (int)(800000.*dz) + 
      (int)(fabs(xmin)*10.)     + (int)(fabs(ymin)*101.)      + (int)(fabs(zmin)*1001.) +
      (int)(fabs(xmax)*999999.) + (int)(fabs(ymax)*30000.)    + (int)(fabs(zmax)*10100.)+
      (int)(xmax*xmax*191919. ) + (int)(ymax*ymax*30220. )    + (int)(zmax*zmax*20202.) +nStatic;    
    
    printf("./rns_initialdata %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,cctk_lsh[0],cctk_lsh[1],cctk_lsh[2],filename_checksum);

  } else if (use_new_bhns_initial_data != -1){
           
    int NS,NU;
    double COMEGA, ROTFUNC, RADEQUAT_TEMP,n;
    
    double *ri, *mui, **pol_p0, **pol_alphi, **pol_rhoi, **pol_gami, **pol_omgi;
    double **pol_omgv, **pol_vOm;
    double **dr_p0, **dt_p0, **drt_p0;
    double **dr_alphi, **dt_alphi, **drt_alphi;
    double **dr_rhoi, **dt_rhoi, **drt_rhoi;
    double **dr_gami, **dt_gami, **drt_gami;
    double **dr_omgi, **dt_omgi, **drt_omgi;
    double **dr_omgv, **dt_omgv, **drt_omgv;
    ifstream pdat("rns.dat");
    if( !pdat) {
      cerr << "Error opening RNS data file" << endl;
      exit (1);
    }
    double Mass = 0.2;
    char buf[100],c;
    double RestMass, AngM;
    pdat.get(buf,100,'='); pdat.get(c); pdat >> RestMass;
    pdat.get(buf,100,'='); pdat.get(c); pdat >> AngM;
    pdat.get(buf,100,'='); pdat.get(c); pdat >> Mass;
    
    pdat >> NS >> NU;
    pdat >> RADEQUAT_TEMP >> n >> COMEGA >> ROTFUNC;
    cout << "RestMass = " << RestMass << endl;
    cout << "Ang Mom = " << AngM << endl;
    cout << "ADM mass = " << Mass << endl;
    cout << "NS, NU = " << NS << "   " << NU << endl;
    cout << "COMEGA = " << COMEGA
	 << ", ROTFUNC = " << ROTFUNC
	 << ", n = " << n << endl;
    /*
     *n_poly = n;
     printf("N = %e\n",*n_poly);
     *RADEQUAT = RADEQUAT_TEMP;
     printf("N = %e\n",*n_poly);
    */
    Pr[CCTK_GFINDEX3D(cctkGH,0,0,0)] = n;
    Pr[CCTK_GFINDEX3D(cctkGH,0,0,1)] = RADEQUAT_TEMP;
    ROTFUNC = ROTFUNC/RADEQUAT_TEMP;
    // Setup hybrid EOS_magnetar_bhns table
    Setup_Hybrid_EOS_magnetar_bhns(K_poly, n, "eos.tab",neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,ergo_star);
    //---------------------------


    //---------------------------
    
    ri = new double[NS];
    mui = new double[NU+1];
    pol_p0 = new double*[NS];
    pol_alphi = new double*[NS];
    pol_rhoi = new double*[NS];
    pol_gami = new double*[NS];
    pol_omgi = new double*[NS];
    pol_omgv = new double*[NS];
    
    dr_p0 = new double*[NS];
    dr_alphi = new double*[NS];
    dr_rhoi = new double*[NS];
    dr_gami = new double*[NS];
    dr_omgi = new double*[NS];
    dr_omgv = new double*[NS];
    
    dt_p0 = new double*[NS];
    dt_alphi = new double*[NS];
    dt_rhoi = new double*[NS];
    dt_gami = new double*[NS];
    dt_omgi = new double*[NS];
    dt_omgv = new double*[NS];
    
    drt_p0 = new double*[NS];
    drt_alphi = new double*[NS];
    drt_rhoi = new double*[NS];
    drt_gami = new double*[NS];
    drt_omgi = new double*[NS];
    drt_omgv = new double*[NS];
    
    for(int i=0; i<NS; i++) {
      pol_p0[i] = new double[NU+1];
      pol_alphi[i] = new double[NU+1];
      pol_rhoi[i] = new double[NU+1];
      pol_gami[i] = new double[NU+1];
      pol_omgi[i] = new double[NU+1];
      pol_omgv[i] = new double[NU+1];
      
      dr_p0[i] = new double[NU+1];
      dr_alphi[i] = new double[NU+1];
      dr_rhoi[i] = new double[NU+1];
      dr_gami[i] = new double[NU+1];
      dr_omgi[i] = new double[NU+1];
      dr_omgv[i] = new double[NU+1];
      
      dt_p0[i] = new double[NU+1];
      dt_alphi[i] = new double[NU+1];
      dt_rhoi[i] = new double[NU+1];
      dt_gami[i] = new double[NU+1];
      dt_omgi[i] = new double[NU+1];
      dt_omgv[i] = new double[NU+1];
      
      drt_p0[i] = new double[NU+1];
      drt_alphi[i] = new double[NU+1];
      drt_rhoi[i] = new double[NU+1];
      drt_gami[i] = new double[NU+1];
      drt_omgi[i] = new double[NU+1];
      drt_omgv[i] = new double[NU+1];
    }
    for(int i=0; i<NS; i++) {
      for(int j=0; j<=NU; j++) {
	pdat >> ri[i] >> mui[j] >> pol_p0[i][j] >> pol_alphi[i][j]
	     >> pol_rhoi[i][j] >> pol_gami[i][j] >> pol_omgi[i][j]
	     >> pol_omgv[i][j];
      }
      if(i%10==0) printf("Reading equilibrium star input file: %d%% finished...\n",(int)((double)i/NS*100.0));
    }
    for(int i=0; i<NS; i++)
      for(int j=0; j<=NU; j++) {
	pol_omgi[i][j] = pol_omgi[i][j]/RADEQUAT_TEMP;
	pol_omgv[i][j] = pol_omgv[i][j]/RADEQUAT_TEMP;
      }
    polar_derivs_bhns(pol_p0,dr_p0,dt_p0,drt_p0,NS,NU,RADEQUAT_TEMP);
    polar_derivs_bhns(pol_alphi,dr_alphi,dt_alphi,drt_alphi,NS,NU,RADEQUAT_TEMP);
    polar_derivs_bhns(pol_rhoi,dr_rhoi,dt_rhoi,drt_rhoi,NS,NU,RADEQUAT_TEMP);
    polar_derivs_bhns(pol_gami,dr_gami,dt_gami,drt_gami,NS,NU,RADEQUAT_TEMP);
    polar_derivs_bhns(pol_omgi,dr_omgi,dt_omgi,drt_omgi,NS,NU,RADEQUAT_TEMP);
    polar_derivs_bhns(pol_omgv,dr_omgv,dt_omgv,drt_omgv,NS,NU,RADEQUAT_TEMP);
    
    double L[3][3];
    double R,Rho,R1,R2,phi1,phi2,t,u;
    double en_den,alp,rho_m,gam,omg_m,omg_v;
    int l,m;
    double grr,gtt,gpp;
    double riso;
    
    int istart = 0;
    int jstart = 0;
    int kstart = 0;
    int iend = cctk_lsh[0];
    int jend = cctk_lsh[1];
    int kend = cctk_lsh[2];
    
    /*  
	printf("87, etc: HELLO!?\t%e\t%e\t%e\t%e\n",x[87],x[88],x[89],x[90]);
	printf("87, etc: HELLO!?\t%d\t%d\t%d\t%d\t%d\n",CCTK_GFINDEX3D(cctkGH,15,1,1),CCTK_GFINDEX3D(cctkGH,16,1,1),CCTK_GFINDEX3D(cctkGH,17,1,1),CCTK_GFINDEX3D(cctkGH,18,1,1),cctk_lsh[0]);
	printf("139, etc: HELLO!?\t%e\t%e\t%e\t%e\n",x[139],x[140],x[141],x[142]);
    */
    for (int k=kstart; k<kend; k++) for (int j=jstart; j<jend; j++) for (int i=istart; i<iend; i++) {
	  int vindex = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  R=sqrt(x[vindex]*x[vindex]+y[vindex]*y[vindex]+z[vindex]*z[vindex]);
	  riso = PhysicalRadius[vindex];
	  Rho = sqrt( x[vindex]*x[vindex]+y[vindex]*y[vindex] );
	  
	  //Here's a fix in case the grid includes the origin point: (0,0,0).  FIXME: THIS FIX DOES NOT WORK YET!
	  //if(fabs(R)<CCTK_DELTA_SPACE(0)*1e-8 || fabs(riso)<CCTK_DELTA_SPACE(0)*1e-8) { R=1e-14; riso=1e-15; }
	  
	  double sint = Rho/R;
	  
	  //Here's another fix in case the grid includes the origin point: (0,0,0).  FIXME: THIS FIX DOES NOT WORK YET!
	  //if(fabs(Rho)<CCTK_DELTA_SPACE(0)*1e-8) { Rho=1e-16; }
	  
	  double cost = sqrt(1.0 - sint*sint);
	  double theta_ang = acos(cost);
	  double phi_ang = atan2(y[vindex],x[vindex]);
	  //
	  // Interpolate metric functions to gridpoint
	  //
	  en_den = bicub_int_bhns(pol_p0,NS,NU,RADEQUAT_TEMP,riso,theta_ang,
				  dr_p0,dt_p0,drt_p0);
	  alp    = bicub_int_bhns(pol_alphi,NS,NU,RADEQUAT_TEMP,riso,theta_ang,
				  dr_alphi,dt_alphi,drt_alphi);
	  rho_m  = bicub_int_bhns(pol_rhoi,NS,NU,RADEQUAT_TEMP,riso,theta_ang,
				  dr_rhoi,dt_rhoi,drt_rhoi);
	  gam    = bicub_int_bhns(pol_gami,NS,NU,RADEQUAT_TEMP,riso,theta_ang,
		       dr_gami,dt_gami,drt_gami);
	  omg_m  = bicub_int_bhns(pol_omgi,NS,NU,RADEQUAT_TEMP,riso,theta_ang,
				  dr_omgi,dt_omgi,drt_omgi);
	  omg_v  = bicub_int_bhns(pol_omgv,NS,NU,RADEQUAT_TEMP,riso,theta_ang,
				  dr_omgv,dt_omgv,drt_omgv);
	  
	  rho_b[vindex] = en_den;
	  P[vindex]     = rho_m;
	  vy[vindex] = omg_v;	    
	  

	  //	    
	  // Get physical metric in Cartesian coordinates
    //
	  L[0][0] = x[vindex]/R;
	  L[0][1] = y[vindex]/R;
	  L[0][2] = z[vindex]/R;
	  L[1][0] = x[vindex]*z[vindex]/(R*R * Rho);
	  L[1][1] = y[vindex]*z[vindex]/(R*R * Rho);
	  L[1][2] = - Rho/(R*R);
	  L[2][0] = - y[vindex]/(Rho*Rho);
	  L[2][1] = x[vindex]/(Rho*Rho);
	  L[2][2] = 0.0;
	  
	  //WARNING: DO NOT USE thetaangle as defined below in the bicub_int statements!
	  //         This screws up equatorial symmetries.
	  double thetaangle=atan2(Rho,z[vindex]);
	  double phiangle=atan2(y[vindex],x[vindex]);
	  
	  double sintheta = Rho/R; //sin(thetaangle);
	  thetaangle=asin(sintheta);
	  double costheta = cos(thetaangle);
	  double sinphi = sin(phiangle);
	  double cosphi = cos(phiangle);
	  
	  L[0][0] = sintheta*cosphi;
	  L[0][1] = sintheta*sinphi;
	  L[0][2] = costheta;
	  L[1][0] = costheta*cosphi/R;
	  L[1][1] = costheta*sinphi/R;
	  L[1][2] = - sintheta/R;
	  L[2][0] = - sinphi/(R*sintheta);
	  L[2][1] = cosphi/(R*sintheta);
	  L[2][2] = 0.0;
	  
	  grr = exp(2.0*alp);
	  grr = grr*RadiusDerivative[vindex]*RadiusDerivative[vindex];
	  gtt = riso*riso*exp(2.0*alp);
	  gpp = riso*riso*sintheta*sintheta*exp(gam-rho_m);
	  gxx[vindex] = pow(L[0][0],2)*grr 
	    + pow(L[1][0],2)*gtt + pow(L[2][0],2)*gpp;
	  //if(x[vindex]>0 && y[vindex]>0 && z[vindex]>0) printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",x[vindex],y[vindex],z[vindex],gxx[vindex],pow(L[0][0],2),grr,gtt,gpp);
	  gxy[vindex] = L[0][0]*L[0][1]*grr 
	    + L[1][0]*L[1][1]*gtt + L[2][0]*L[2][1]*gpp;
	  gxz[vindex] = L[0][0]*L[0][2]*grr 
	    + L[1][0]*L[1][2]*gtt + L[2][0]*L[2][2]*gpp; 
	  if(i==-17 && j==17 && k==4) {
	    for(int ww=0;ww<3;ww++) for(int xx=0;xx<3;xx++) {
		printf("%d %d: %.17e\n",ww,xx,L[ww][xx]);
	      }
	    printf("gxz[vindex]=%.17e... %.17e\n",gxz[vindex],L[0][0]*L[0][2]*grr 
		   + L[1][0]*L[1][2]*gtt + L[2][0]*L[2][2]*gpp);
	    printf("1 = %.17e\n",L[0][0]*L[0][2]*grr);
	    printf("2 = %.17e\n",L[1][0]*L[1][2]*gtt);
	    printf("3 = %.17e\n",L[2][0]*L[2][2]*gpp);
	    printf("grr,gtt,gpp = %e,%e,%e\n",grr,gtt,gpp);
	    exit(1);
	  }
	  gyy[vindex] = pow(L[0][1],2)*grr 
	    + pow(L[1][1],2)*gtt + pow(L[2][1],2)*gpp;
	  gyz[vindex] = L[0][1]*L[0][2]*grr 
	    + L[1][1]*L[1][2]*gtt + L[2][1]*L[2][2]*gpp;
	  gzz[vindex] = pow(L[0][2],2)*grr 
	    + pow(L[1][2],2)*gtt + pow(L[2][2],2)*gpp;
	  if(fabs(sintheta)<1e-8) { printf("HI. %e %e %e %e %e %e\n",gxx[vindex],gxy[vindex],gxz[vindex],gyy[vindex],gyz[vindex],gzz[vindex]); }
	  //
	  // Set gauge functions
	  //
	  lapm1[vindex] = exp(0.5*(gam+rho_m)) - 1.0;
	  shiftx[vindex] = y[vindex]*omg_m;
	  shifty[vindex] = -x[vindex]*omg_m;
	  shiftz[vindex] = 0.0;
	  
	  if(gxx[vindex]!=gxx[vindex]) { printf("FOUND A NAN!\n"); exit(1); }
	  //printf("HI. %e %e %e %e %e %e %e %e %e\n",gxx[vindex],gxy[vindex],gxz[vindex],gyy[vindex],gyz[vindex],gzz[vindex],lapm1[vindex],shiftx[vindex],shifty[vindex]);
	  
	  //if(fabs(Rho)<1e-8) { printf("HI. %e %e %e %e %e %e %e %e %e\n",gxx[vindex],gxy[vindex],gxz[vindex],gyy[vindex],gyz[vindex],gzz[vindex],lapm1[vindex],shiftx[vindex],shifty[vindex]); }
	}
    printf("HEWWO!\t%e\n",gxz[CCTK_GFINDEX3D(cctkGH,2,1,1)]);

    //
    // Free memory
    //
    /*
      delete ri;
      delete mui;
      delete pol_p0;
      delete pol_alphi;
      delete pol_rhoi;
      delete pol_gami;
      delete pol_omgi;
      delete pol_omgv;
    */
  }
}

void Setup_Hybrid_EOS_magnetar_bhns(double K_poly, double n, char *eos_file,
				    int neos,double *rho_tab,double *P_tab,
				    double *eps_tab,double *k_tab,double *gamma_tab, int ergo_star) {
  using namespace std;
  ifstream eosdat;
  eosdat.open(eos_file);
  if (!eosdat) {
    neos = 1;
    cout << "Can't open " << eos_file << endl;
    cout << "Assuming polytropic EOS..." << endl;
    cout << "Polytropic index = " << n << endl;
    cout << "Polytropic constant K = " << K_poly << endl;
    //rho_tab = new double[1]; P_tab = new double[1]; eps_tab = new double[1];
    //gamma_tab = new double[2]; k_tab = new double[2];
    rho_tab[0] = 0.03;
    P_tab[0] = K_poly * pow(rho_tab[0],1.0+1.0/n);
    eps_tab[0] = n*P_tab[0]/rho_tab[0];
    gamma_tab[0] = 1.0 + 1.0/n; k_tab[0] = K_poly;
    gamma_tab[1] = 1.0 + 1.0/n; k_tab[1] = K_poly;
  } else {
    char buf[100],c;
    int iphase;
    double det,rho_ref,p_ref;

    eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> neos;
    cout << "Assuming polytropic EOS..." << endl;

    // read gammas
    for (int i=0; i<neos+1; i++) {
      eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> gamma_tab[i];}



    // read the intermedium pieces
   
    for (int i=0; i<neos+1; i++) {
      eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> rho_tab[i];}
   
    // Lunan: if ergo_star == 0, i should go from 0 to neos. Otherwise from 0 to neos+1.                                          
    // Also, you should set neos-1 pieces of density in your eos.tab if ergo_star == 0.                    
    
    /*
    if (ergo_star == 0){
      cout << "Not an ergo star!" << endl;
      for (int i=0; i<neos; i++) {
      eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> rho_tab[i];}
    }
    else{
      cout << "Ergo star!" << endl;
      for (int i=0; i<neos+1; i++) {
      eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> rho_tab[i];}
    }
    */

    // read reference restmass density, pressure
      eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> rho_ref;
      eosdat.get(buf,100,'='); eosdat.get(c); eosdat >> p_ref;

      // compute kappa for ref press and density

      iphase = 1;
      for(int i=1; i<=neos; i++)
	{det = (rho_ref-rho_tab[i])*(rho_ref-rho_tab[i-1]);
	  if (det <= 0.0){
	      iphase = i;
	      break ;}}

      k_tab[iphase] = p_ref/pow(rho_ref,gamma_tab[iphase]);

      // compute kappa's
      if (iphase > 0){
	  for(int j=iphase-1; j>=0; j--){
	      k_tab[j] = pow(rho_tab[j], gamma_tab[j+1]-gamma_tab[j])*k_tab[j+1];}}
      if (iphase < neos){
	  for(int j=iphase+1; j<=neos; j++){
	      k_tab[j] = pow(rho_tab[j-1], gamma_tab[j-1]-gamma_tab[j]) *k_tab[j-1];}}

      //    cout << "convert variables into c=G=Msun=1" << endl;
    
    // kappas
    //for(int i=0; i<neos+1; i++){
    //  double fkappa = 1.0/c2*pow(c6/(g3*msol2), gamma_tab[i]-1.0);
    //  k_tab[i] = k_tab[i]*fkappa;}
    
    // densities
    //for(int i=0; i<neos+1; i++){
    //  double frho = g3*msol2/c6; 
      // double fpre = g3*msol2/c8;
      //rho_tab[i] = rho_tab[i]*frho;}
      //P_tab[i]   = P_tab[i]*fpre;}
    
    // compute pressure 
    for(int i=0; i<neos+1; i++){
      P_tab[i]   = k_tab[i]*pow(rho_tab[i],gamma_tab[i]);}
    
    // internal energy
    eps_tab[0] = P_tab[0]/rho_tab[0]/(gamma_tab[0]-1.0);
    
    for(int i=1; i<neos+1;i++){
      eps_tab[i] = eps_tab[i-1] + (P_tab[i]/rho_tab[i] - P_tab[i-1]/rho_tab[i-1])/(gamma_tab[i]-1.0);
    }

    /*    cout << "neos = " <<  neos <<endl;
    for(int i=0; i<neos+1; i++){
      cout << "gamma_tab = " <<  gamma_tab[i] <<endl;
      cout << "rho_tab = "   <<  rho_tab[i]   <<endl;
      cout << "P_tab   = "   <<  P_tab[i]     <<endl;
      cout << "k_tab   = "   <<  k_tab[i]     <<endl;
      cout << "eps_tab = "   <<  eps_tab[i]   <<endl;} 
    */
  } 
 eosdat.close();
}

double bicub_int_bhns(double **f, int NS, int NU,
		      double RADEQUAT,
		      double riso, double theta,
		      double **f1, double **f2, double **f12)
{
  int l,m;
  double cost = cos(theta);
  l = int(riso*NS/(RADEQUAT+riso));
  m = int(cost*NU);
  if(m > NU-1) m = NU-1;
  if (l>NS-2) {
    printf("Outer boundary too far away!\n");
    printf("Set a smaller outer boundary so that a proper interpolation can be done from data in rns.dat.\n");
    exit(1);
  }
  double R1,R2,thet1,thet2;
  R1 = RADEQUAT*l/(NS-l);
  R2 = RADEQUAT*(l+1)/(NS-l-1);
  thet1 = acos(double(m)/double(NU));
  thet2 = acos(double(m+1)/double(NU));
  double y[5], y1[5], y2[5], y12[5];
  // function
  y[1] = f[l][m];       y[2] = f[l+1][m];
  y[3] = f[l+1][m+1];   y[4] = f[l][m+1];
  // r derivatives
  y1[1] = f1[l][m];       y1[2] = f1[l+1][m];
  y1[3] = f1[l+1][m+1];   y1[4] = f1[l][m+1];
  // theta derivatives
  y2[1] = f2[l][m];       y2[2] = f2[l+1][m];
  y2[3] = f2[l+1][m+1];   y2[4] = f2[l][m+1];
  // mixed derivatives
  y12[1] = f12[l][m];       y12[2] = f12[l+1][m];
  y12[3] = f12[l+1][m+1];   y12[4] = f12[l][m+1];
  
  // interpolate
  double f_out,f1_out,f2_out;
  bcuint_bhns(y,y1,y2,y12,R1,R2,thet1,thet2,riso,theta,
	      f_out,f1_out,f2_out);
  return f_out;
}


void polar_derivs_bhns(double **f, double **f1, double **f2,
		       double **f12,
		  int NS, int NU, double RADEQUAT)
{
  int l,m;
  double *rg, *thg;
  rg = new double[NS];
  thg = new double[NU+1];
  for(l=0; l<NS; l++)
    rg[l] = RADEQUAT*l/(NS-l);
  for(m=0; m<=NU; m++)
    thg[m] = acos(double(m)/double(NU));
  
  // r derivatives
  for(m=0; m<=NU; m++) {
    f1[0][m] = (f[1][m] - f[0][m])/(rg[1]-rg[0]);
    for(l=1; l<NS-1; l++)
      f1[l][m] = (f[l+1][m]-f[l-1][m])/(rg[l+1]-rg[l-1]);
    f1[NS-1][m] = (f[NS-1][m]-f[NS-2][m])/(rg[NS-1]-rg[NS-2]);
  }
  // theta derivatives
  for(l=0; l<NS; l++) {
    f2[l][0] = (f[l][1] - f[l][0])/(thg[1]-thg[0]);
    for(m=1; m<NU; m++)
      f2[l][m] = (f[l][m+1]-f[l][m-1])/(thg[m+1]-thg[m-1]);
    //      f2[l][NU] = (f[l][NU]-f[l][NU-1])/(thg[NU]-thg[NU-1]);
    f2[l][NU] = 0.0;
  }
  // mixed derivatives
  for(l=0; l<NS; l++) {
    f12[l][0] = (f1[l][1] - f1[l][0])/(thg[1]-thg[0]);
    for(m=1; m<NU; m++)
      f12[l][m] = (f1[l][m+1]-f1[l][m-1])/(thg[m+1]-thg[m-1]);
    //      f12[l][NU] = (f1[l][NU]-f1[l][NU-1])/(thg[NU]-thg[NU-1]);
    f12[l][NU] = 0.0;
  }
}

void bcuint_bhns(double *y, double *y1, double *y2, double *y12, 
		 double x1l,double x1u, double x2l, double x2u, 
		 double x1, double x2,double & ansy, double & ansy1, double & ansy2) {
  /*   void bcucof_bhns(double *y, double *y1, double *y2, double *y12, */
  /*               double d1, double d2, double **c); */
  int i;
  double t,u,d1,d2,**c;
  
  c=dmatrix(1,4,1,4);
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof_bhns(y,y1,y2,y12,d1,d2,c);
  if (x1u == x1l || x2u == x2l) { printf("ERROR: %e %e %e %e\n",x1u,x1l,x2u,x2l); nrerror("Bad input in routine bcuint"); }
  t=(x1-x1l)/d1;
  u=(x2-x2l)/d2;
  ansy = ansy2 = ansy1 = 0.0;
  for (i=4;i>=1;i--) {
    ansy=t*ansy + ((c[i][4]*u+c[i][3])*u+c[i][2])*u + c[i][1];
    ansy2=t*ansy2 + (3.0*c[i][4]*u+2.0*c[i][3])*u + c[i][2];
    ansy1=u*ansy1 + (3.0*c[4][i]*t+2.0*c[3][i])*t + c[2][i];
  }
  ansy1 /= d1;
  ansy2 /= d2;
  free_dmatrix(c,1,4,1,4);
}

void bcucof_bhns(double *y, double *y1, double *y2, double *y12,
		 double d1, double d2, double **c) {
  static int wt[16][16]=
    { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
      -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
      2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
      0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
      0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
      0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
      -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
      9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
      -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
      2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
      -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
      4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};
  int l,k,j,i;
  double xx,d1d2,cl[16],x[16];
  
  d1d2=d1*d2;
  for (i=1;i<=4;i++) {
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) {
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++)
    for (j=1;j<=4;j++) c[i][j]=cl[l++];
}


double findmin_bhns(double *array, int &size) {

  double smallest;
  smallest = array[0];
  for (int i=1; i<=size-1; i++) 
      if (array[i] < smallest) 
	{
	  smallest=array[i];
	}
  return smallest;
}

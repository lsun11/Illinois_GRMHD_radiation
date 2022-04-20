//-----------------------------------------------------------------------
// $Id$
//-----------------------------------------------------------------------
// Read binary files and do fancy things with them...
//-----------------------------------------------------------------------
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

using namespace std;

extern "C" void CCTK_FCALL read_inputfile_bbh_cookpfeiffer_rot_
  (int *genID_cmdline_output_enable,int *gridnumber,
   int *Nradial, int *Ntheta, int *Nphi,
   double *radmin, double *radmax, double *xoffset,
   double *K_rr_rot,double *K_rth_rot,double *K_rp_rot,double *K_thth_rot,double *K_thp_rot, double *K_pp_rot,
   double *shiftr_rot,double *shiftth_rot,double *shiftp_rot,double *phi_rot,double *lapm1_rot);

extern "C" void read_inputfile_bbh_cookpfeiffer_rot(int genID_cmdline_output_enable,int gridnumber,
						    int Nradial, int Ntheta, int Nphi,
						    double radmin, double radmax, double xoffset,
						    double *K_rr_rot,double *K_rth_rot,double *K_rp_rot,double *K_thth_rot,double *K_thp_rot,double *K_pp_rot,
						    double *shiftr_rot,double *shiftth_rot,double *shiftp_rot,double *phi_rot,double *lapm1_rot){
  
  //  char filename[100];
  char filename_gxx[100];
  char filename_lapse[100];
  char filename_shiftx[100];
  char filename_shifty[100];
  char filename_shiftz[100];
  char filename_Kxx[100];
  char filename_Kxy[100];
  char filename_Kxz[100];
  char filename_Kyy[100];
  char filename_Kyz[100];
   
  //sprintf(filename,"outfields_bbhcp_grid-%d.dat",ext[1]+ext[2]+ext[0]*CCTK_MyGrid(cctkGH)+(int)(119831*dx));
  sprintf(filename_gxx,"outfields_bbhcp_grid_%d-Nid_gxx.dat",gridnumber);
  sprintf(filename_lapse,"outfields_bbhcp_grid_%d-Nid_N.dat",gridnumber);
  sprintf(filename_shiftx,"outfields_bbhcp_grid_%d-Nid_Nx.dat",gridnumber);
  sprintf(filename_shifty,"outfields_bbhcp_grid_%d-Nid_Ny.dat",gridnumber);
  sprintf(filename_shiftz,"outfields_bbhcp_grid_%d-Nid_Nz.dat",gridnumber);
  sprintf(filename_Kxx,"outfields_bbhcp_grid_%d-Nid_Kxx.dat",gridnumber);
  sprintf(filename_Kxy,"outfields_bbhcp_grid_%d-Nid_Kxy.dat",gridnumber);
  sprintf(filename_Kxz,"outfields_bbhcp_grid_%d-Nid_Kxz.dat",gridnumber);
  sprintf(filename_Kyy,"outfields_bbhcp_grid_%d-Nid_Kyy.dat",gridnumber);
  sprintf(filename_Kyz,"outfields_bbhcp_grid_%d-Nid_Kyz.dat",gridnumber);
	  
  
  //  sprintf(filename,"outfields_bbhcp_proc-%d.dat",ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(ymin));
  //sprintf(filename,"outfields_bbhcp_grid-%d.dat",gridnumber);
  
  //if(genID_cmdline_output_enable==1) {
  //  printf("./getbbhID-fish-parallel %e %e %e %e %e %e %d %d %d %d\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],ext[1]*3015+ext[2]*9121+ext[0]*CCTK_MyProc(cctkGH)+(int)(1198311*dx)-(int)(ymin));
  //} else {
  //    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);
   
    
    ifstream infile_gxx;
    ifstream infile_lapse;
    ifstream infile_shiftx;
    ifstream infile_shifty;
    ifstream infile_shiftz;
    ifstream infile_Kxx;
    ifstream infile_Kxy; 
    ifstream infile_Kxz; 
    ifstream infile_Kyy; 
    ifstream infile_Kyz;

    infile_gxx.open(filename_gxx);  
    infile_lapse.open(filename_lapse); 
    infile_shiftx.open(filename_shiftx); 
    infile_shifty.open(filename_shifty);  
    infile_shiftz.open(filename_shiftz);
    infile_Kxx.open(filename_Kxx);  
    infile_Kxy.open(filename_Kxy); 
    infile_Kxz.open(filename_Kxz);
    infile_Kyy.open(filename_Kyy);
    infile_Kyz.open(filename_Kyz);

    
    if(!infile_gxx) {
      cerr << "\a Can't open " << filename_gxx << " for input." << endl;
      exit(1);
    }
    

    
    double PI=acos(-1.0);
    
    double thetamin=-PI/1000.;
    double thetamax=PI + PI/1000.; 
    double phimin=-PI/1000.; 
    double phimax=2.0*PI+PI/1000.; 
    
    double logradmin=log(radmin);
    double logradmax=log(radmax);

    int logarithmic;

    if (gridnumber==3){
      logarithmic=1;
    } else{
      logarithmic=0;
    }

    double radial_coord;
    for (int k = 0; k < Nradial; k++) for (int j = 0; j < Nphi; j++) for (int i = 0; i < Ntheta; i++) {
      if (logarithmic==1){
	radial_coord = exp(logradmin+(logradmax-logradmin)/(Nradial-1.0)*k);
      } else {
	radial_coord = radmin+(radmax-radmin)/(Nradial-1.0)*k;
      }
      double theta = thetamin+(thetamax-thetamin)/(Ntheta-1.0)*i; 
      double phi = phimin+(phimax-phimin)/(Nphi-1.0)*j; 

      double xx = radial_coord*sqrt(1.0-cos(theta)*cos(theta))*cos(phi)+xoffset; 
      double yy = radial_coord*sqrt(1.0-cos(theta)*cos(theta))*sin(phi); 
      double zz = radial_coord*cos(theta); 
      
      double radius=sqrt(xx*xx+yy*yy+zz*zz)+1.e-10; 
      double rad_cyl=sqrt(xx*xx+yy*yy)+1.e-10; 
      
      //      cout << "radius: " << radius << " gridnumber: " << gridnumber << endl;
      double Lxp=-yy; 
      double Lyp=xx; 
      double Lzp=0.0; 
      //watch out for divide by zero 
      double Lxth=xx*zz/rad_cyl; 
      double Lyth=yy*zz/rad_cyl; 
      double Lzth=-rad_cyl; 
      double Lxr=xx/radius; 
      double Lyr=yy/radius;
      double Lzr=zz/radius; 

      double Lrx = xx/radius;
      double Lry = yy/radius;
      double Lrz = zz/radius;
      double Lpx = -yy/rad_cyl/rad_cyl;
      double Lpy = xx/rad_cyl/rad_cyl;
      double Lpz = 0.0;
      double Lthx = xx*zz/(radius*radius*rad_cyl);
      double Lthy = yy*zz/(radius*radius*rad_cyl);
      double Lthz = -rad_cyl/radius/radius;
        

      double gxxL;
      double KxxL;
      double KxyL;
      double KxzL;
      double KyyL;
      double KyzL;
      double shiftxL;
      double shiftyL;
      double shiftzL;
      double lapseL;
    
      
      //make sure this is set correctly
      int index = i+Ntheta*(j+Nphi*k);
      
      infile_Kxx.read((char *) &KxxL, sizeof(double));
      infile_Kxy.read((char *) &KxyL, sizeof(double));
      infile_Kxz.read((char *) &KxzL, sizeof(double));
      infile_Kyy.read((char *) &KyyL, sizeof(double));
      infile_Kyz.read((char *) &KyzL, sizeof(double));
      infile_shiftx.read((char *) &shiftxL, sizeof(double)); 
      infile_shifty.read((char *) &shiftyL, sizeof(double));
      infile_shiftz.read((char *) &shiftzL, sizeof(double));
      infile_gxx.read((char *) &gxxL, sizeof(double));
      infile_lapse.read((char *) &lapseL, sizeof(double));
      lapm1_rot[index] = lapseL - 1.0;
      
      double psi=pow(gxxL,0.25); 
      phi_rot[index] = log(psi);
      
      double KzzL = - (KxxL + KyyL);
      
      
      K_rr_rot[index] =  Lxr*(Lxr*KxxL+Lyr*KxyL+Lzr*KxzL)
	+Lyr*(Lxr*KxyL+Lyr*KyyL+Lzr*KyzL) 
	+Lzr*(Lxr*KxzL+Lyr*KyzL+Lzr*KzzL);
        
      K_rth_rot[index] = Lxth*(Lxr*KxxL+Lyr*KxyL+Lzr*KxzL) 
	+Lyth*(Lxr*KxyL+Lyr*KyyL+Lzr*KyzL) 
	+Lzth*(Lxr*KxzL+Lyr*KyzL+Lzr*KzzL);
        
      K_rp_rot[index] =  Lxp*(Lxr*KxxL+Lyr*KxyL+Lzr*KxzL)
	+Lyp*(Lxr*KxyL+Lyr*KyyL+Lzr*KyzL) 
	+Lzp*(Lxr*KxzL+Lyr*KyzL+Lzr*KzzL);
        
      K_thth_rot[index] = Lxth*(Lxth*KxxL+Lyth*KxyL+Lzth*KxzL) 
	+Lyth*(Lxth*KxyL+Lyth*KyyL+Lzth*KyzL) 
	+Lzth*(Lxth*KxzL+Lyth*KyzL+Lzth*KzzL); 
        
      K_thp_rot[index] = Lxp*(Lxth*KxxL+Lyth*KxyL+Lzth*KxzL) 
	+Lyp*(Lxth*KxyL+Lyth*KyyL+Lzth*KyzL) 
	+Lzp*(Lxth*KxzL+Lyth*KyzL+Lzth*KzzL);
        
      K_pp_rot[index] = Lxp*(Lxp*KxxL+Lyp*KxyL+Lzp*KxzL) 
	+Lyp*(Lxp*KxyL+Lyp*KyyL+Lzp*KyzL) 
	+Lzp*(Lxp*KxzL+Lyp*KyzL+Lzp*KzzL);
      
      shiftr_rot[index] = Lrx*shiftxL +Lry*shiftyL +Lrz*shiftzL;
      shiftth_rot[index] = Lthx*shiftxL +Lthy*shiftyL +Lthz*shiftzL;
      shiftp_rot[index] = Lpx*shiftxL +Lpy*shiftyL +Lpz*shiftzL;
    }
    printf("==========================================\n");
    printf("Okay, the initial data has been read in...\n");
    printf("==========================================\n");

    infile_gxx.close();                             
    infile_lapse.close();                           
    infile_shiftx.close();                          
    infile_shifty.close();                          
    infile_shiftz.close();                          
    infile_Kxx.close();                             
    infile_Kxy.close();                             
    infile_Kxz.close();                             
    infile_Kyy.close();                             
    infile_Kyz.close();   
    //}
}

extern "C" void CCTK_FCALL read_inputfile_bbh_cookpfeiffer_rot_
  (int *genID_cmdline_output_enable,int *gridnumber,
   int *Nradial, int *Ntheta, int *Nphi,
   double *radmin, double *radmax, double *xoffset,
   double *K_rr_rot,double *K_rth_rot,double *K_rp_rot,double *K_thth_rot,double *K_thp_rot,double *K_pp_rot,
   double *shiftr_rot,double *shiftth_rot,double *shiftp_rot,double *phi_rot,double *lapm1_rot)
{  
  read_inputfile_bbh_cookpfeiffer_rot(*genID_cmdline_output_enable,*gridnumber,
				      *Nradial, *Ntheta, *Nphi,
				      *radmin, *radmax, *xoffset,
				      K_rr_rot,K_rth_rot,K_rp_rot,K_thth_rot,K_thp_rot,K_pp_rot,
				      shiftr_rot,shiftth_rot,shiftp_rot,phi_rot,lapm1_rot);
}

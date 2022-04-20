//-----------------------------------------------------------------------
// Read WDNS bin. files and do fancy things with them...
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

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_wdns)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
   double &kappa3, double &rhoo1, double &rhoo2);

extern "C" void read_inputfile_wdns(const cGH *cctkGH,int genID_cmdline_output_enable,int fisheye_enable,
				    int reset_shift_lapse,double xmin,double ymin,double zmin,
				    double dx,double dy,double dz,int *ext,
				    double *lapse,double *shiftx,double *shifty,double *shiftz,
				    double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
				    double *rho_nbar,double *ux,double *uy, double *uz,
				    double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
				    double &kappa3, double &rhoo1, double &rhoo2){

  static int nStatic;
  nStatic += 1;
  int filename_checksum = ext[1]*301+ext[2]*91221+ext[0] + 1241*CCTK_MyProc(cctkGH)+(int)(99999*dx)+ nStatic;
  char filename[100];
  double xrot, Omega;
  sprintf(filename,"CTS_bin-proc%d.d",filename_checksum);

  ofstream solfile;




  if(genID_cmdline_output_enable==1) {
    printf("./wdns_initialdata_isotropic %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d $ascii $order\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],filename_checksum);
  } else {

    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);

    ifstream infile1;
    infile1.open(filename);
    if(!infile1) {
      cerr << "\a Can't open " << filename << " for input." << endl;
      exit(1);
    }

    double doubledum;
    int ntot,intdum, int_nx, int_ny, int_nz;

    // Read grid parameters
    infile1.read((char *) &int_nx, sizeof(int));
    infile1.read((char *) &int_ny, sizeof(int));
    infile1.read((char *) &int_nz, sizeof(int));
    infile1.read((char *) &ntot, sizeof(int));

    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];

    if(ntot!=nx*ny*nz) {
      printf("MAYBE PROBLEM: This data file has %d elements, but current grid has %d!\n",ntot,nx*ny*nz);
      //exit(1);
    } else {
      printf("This data file has %d elements, like we expect!\n",nx*ny*nz);
    }

    //    if(reset_shift_lapse==0) {
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &psi[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &shiftx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &shifty[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &shiftz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kxx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kxy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kxz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kyy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kyz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

      //NEXT GF IS A DUMMY.  We know Kzz from the other Kij's and the fact that trK==0!
	  infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &ux[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &uy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &uz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
      //    }
	}
      }
    }

    // Read WDNS parameters
    infile1.read((char *) &Omega, sizeof(double));
    cout<<"WDNS parameters read"<<endl;
    cout<<"Omega="<<Omega<<endl;
    infile1.read((char *) &xrot, sizeof(double));
    cout<<"xrot="<<xrot<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"M_baryon_NS="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"M_baryon_WD="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"ADM mass="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"Komar mass="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"Ang mom="<<doubledum<<endl;
    infile1.read((char *) &gamma1, sizeof(double));
    infile1.read((char *) &gamma2, sizeof(double));
    infile1.read((char *) &gamma3, sizeof(double));
    infile1.read((char *) &kappa1, sizeof(double));
    infile1.read((char *) &kappa2, sizeof(double));
    infile1.read((char *) &kappa3, sizeof(double));
    infile1.read((char *) &rhoo1, sizeof(double));
    infile1.read((char *) &rhoo2, sizeof(double));

    infile1.close();  


    //    solfile.open("lap.dat");
    for(int k=0;k<nz;k++){
      double z = zmin + k*dz;
      for(int j=0;j<ny;j++){
	double y = ymin + j*dy;
	for(int i=0;i<nx;i++){
	  double x = xmin + i*dx;
	  double rhoo = rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double conf= psi[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double lap = lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double betx = shiftx[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double bety = shifty[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double betz = shiftz[CCTK_GFINDEX3D(cctkGH,i,j,k)];    
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  //	  solfile << x <<" " << y  << " " << y << " " << lap << endl;
  



	  if (rhoo > 1.e-16){
	    double conf4 = conf*conf*conf*conf;
	    double lap2 = lap*lap;
	    double v2 = ( (Omega*y-betx)*(Omega*y-betx)+(Omega*(x-xrot)+bety)*(Omega*(x-xrot)+bety)+betz*betz )*conf4/lap2;
	    double Lore = 1./sqrt(1.-v2);
	    ux[CCTK_GFINDEX3D(cctkGH,i,j,k)] = Lore*(-Omega*y+betx)*conf4/lap;
	    uy[CCTK_GFINDEX3D(cctkGH,i,j,k)]= Lore*( Omega*(x-xrot)+bety )*conf4/lap;
	    uz[CCTK_GFINDEX3D(cctkGH,i,j,k)]= Lore*(betz)*conf4/lap;
	  }
	  else{
	    ux[CCTK_GFINDEX3D(cctkGH,i,j,k)] = 0.;
	    uy[CCTK_GFINDEX3D(cctkGH,i,j,k)]= 0.;
	    uz[CCTK_GFINDEX3D(cctkGH,i,j,k)]= 0.;
	  }	 
	}
      }
    }
    
  }
    
  //  solfile.close();

  /*  printf("an1,an2,an3 = %f %f %f\n",gamma1,gamma2,gamma3);
  printf("Gam1,Gam2,Gam3 = %f %f %f\n",1.+1./gamma1,1.+1./(gamma2),1.+1./(gamma3));
  printf("kappa1,kappa2,kappa3 =%f %f %f\n",kappa1,kappa2,kappa3);
  printf("rhoo1,rhoo2 =%f %f\n",rhoo1,rhoo2);
  exit(1);
  */
}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_wdns)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
   double &kappa3, double &rhoo1, double &rhoo2)
{  
  read_inputfile_wdns(*cctkGH,*genID_cmdline_output_enable,*fisheye_enable,
		      *reset_shift_lapse,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
		      lapse,shiftx,shifty,shiftz,
		      psi,kxx,kxy,kxz,kyy,kyz,
		      rho_nbar,ux,uy, uz,
		      gamma1,gamma2,gamma3,kappa1,kappa2, kappa3,rhoo1,rhoo2);
}


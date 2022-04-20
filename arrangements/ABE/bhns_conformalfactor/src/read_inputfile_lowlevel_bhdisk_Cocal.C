//-----------------------------------------------------------------------
// Read black-hole Cocal bin. files and do fancy things with them...
//-----------------------------------------------------------------------
#include <iostream>
#include <sstream> 
#include <string>  
#include <cstring>

#include <unistd.h>

#include <cmath>
#include <iomanip> 
#include <fstream> 
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_BH_disk)
  (const cGH **cctkGH,int *num_CO, int *piecewise,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *gxx,  double   *gxy,double *gxz,  double *gyy,  double *gyz,  double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
   double *psi, double *trK, double *rho_nbar,double *ux,double *uy, double *uz,
   int *nphase, double *abc, double *abi, double *rhoi, double *tedi);


extern "C" void read_inputfile_BH_disk(const cGH *cctkGH,
				       int num_CO, int piecewise,int genID_cmdline_output_enable,int fisheye_enable,
				       int reset_shift_lapse,double xmin,double ymin,double zmin,double dx,double dy,double dz,int *ext,
				       double *lapse,double *shiftx,double *shifty,double *shiftz,
				       double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				       double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
				       double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
				       double *psi, double *trK, double *rho_nbar,double *ux,double *uy, double *uz,
				       int *nphase, double *abc, double *abi, double *rhoi, double *tedi)
{
  static int nStatic;
  nStatic += 1;

  // define the other boundary
  double xmax=xmin+ext[0]*dx;
  double ymax=ymin+ext[1]*dy;
  double zmax=zmin+ext[2]*dz;

  // define lable for grids
    int filename_checksum = ext[0]*301+ext[1]*10+ext[2]*1000000  + 1241*CCTK_MyProc(cctkGH) + (int)(999999.*dx) + (int)(100100100.*dy) + (int)(800000.*dz) + 
      (int)(fabs(xmin)*10.)     + (int)(fabs(ymin)*101.)      + (int)(fabs(zmin)*1001.) +
      (int)(fabs(xmax)*999999.) + (int)(fabs(ymax)*30000.)    + (int)(fabs(zmax)*10100.)+
      (int)(xmax*xmax*191919. ) + (int)(ymax*ymax*30220. )    + (int)(zmax*zmax*20202.) +nStatic;    


  char filename[100];
  char dumname[100];
  char inputstring[50];
  int nph ;
  nph = nphase[0];

  sprintf(filename,"CTS_bin-proc%d.d",filename_checksum);

  ofstream solfile;

  cout << " ENTER read_inputfile_lowlevel_bhdisk_Cocal.C "<<endl;

  if(num_CO==1)
  {
    if(piecewise==0)     
    {
      cout << "Single star, with simple polytropic EOS." << endl;    
      strncpy( inputstring, "./read_bin_rns", 50);
    }
    else if (piecewise==1) 
    { 
      cout << "Single star, with piecewise polytropic EOS." << endl;    
      strncpy( inputstring, "./read_bin_rns_piecewise_antonios", 50);
    }
    else
    {
      printf("Parameter piecewise must be 0 (simple polytrope) or 1 (piecewise). Exiting...")  ;
      exit (0);
    }
  }
  else if (num_CO==2)
  {
    strncpy( inputstring, "./read_bin_nsns", 50);
  }
  else
  {
    printf("Parameter num_CO must be 1 or 2. Exiting...")  ;
    exit (0);
  }



  if(genID_cmdline_output_enable==1) 
  {
     printf("%s %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d $ascii $order\n",inputstring,xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],filename_checksum);
     
     nph = 2;   // This is not consistent with nph of COCAL but with the neos parameter.
     for(int ia=0; ia<=nph; ia++)
     {
       abc[ia]  = 1.0; 
       abi[ia]  = 1.0;
       rhoi[ia] = 1.0;
       tedi[ia] = 2.0; 
     }
     cout << " nph = " << nph << endl;   
  } 
  else 
  {
    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);

    // define number of points along each direction
    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];
    int ntot,intdum, int_nx, int_ny, int_nz;
    ntot=nx*ny*nz;   // total number of points

    // define auxiliary arrays
    double gtx[ntot], gty[ntot], gtz[ntot];

    ifstream infile1;
    infile1.open(filename);
    if(!infile1) {
      cerr << "\a Can't open " << filename << " for input." << endl;
      exit(1);
    }

    // Read grid parameters
    cout<<"Cocal parameters read"<<endl;
    infile1.read((char *) &int_nx, sizeof(int));
    infile1.read((char *) &int_ny, sizeof(int));
    infile1.read((char *) &int_nz, sizeof(int));
    infile1.read((char *) &ntot, sizeof(int));

    if(ntot!=nx*ny*nz) {
      printf("MAYBE PROBLEM: This data file has %d elements, but current grid has %d!\n",ntot,nx*ny*nz);
      exit(1);
    } else {
      printf("Cocal:: This data file has %d elements, like we expect!\n",nx*ny*nz);
    }

    // Read full metric, extrinsic curvature, rest-mass density and Eulerian velocities
    
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &psi[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &gtx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gty[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gtz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &gxx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gxy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gxz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gyy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gyz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &gzz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &kxx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kxy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kxz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kyy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kyz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &kzz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &ux[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &uy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &uz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    cout << " READING PEOS FILE "<<endl ;
    //EoS parameters
    infile1.read((char *) &nph, sizeof(int));
    cout << "Number of phases :" << nph << endl;
    for(int j=0; j<=nph; j++)
    {
      infile1.read((char *) &abc[j],sizeof(double));
      infile1.read((char *) &abi[j],sizeof(double));
      infile1.read((char *) &rhoi[j],sizeof(double));
      infile1.read((char *) &tedi[j],sizeof(double));
    }

    nphase[0] = nph ;

    cout << "*****************************  READING FROM ID FILE  ****************************************" << endl;
    cout << " Number of phases :" << nph << "  nphase[0]=" << nphase[0] << endl;
    cout << setw(20) << "Kappa: [G=c=Msun=1]" << setw(18) << "Gamma" << setw(20) << "rho0 [G=c=Msun=1]" << setw(18) << "rho [G=c=Msun=1]" << endl;
    for(int j=0; j<=nph; j++)
      cout << scientific << setprecision(10) << setw(20) << abc[j] << setw(18) <<
              abi[j] << setw(20) << rhoi[j] << setw(18) << tedi[j] << endl;
    cout << "*********************************************************************************************" << endl;
    cout << "Rest mass density=rho0, total energy density=rho" << endl;


    infile1.close();  


    cout << "Compute the gauge variables" << endl;

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  double Psil = psi[ind];
	  double gtxl = gtx[ind];
	  double gtyl = gty[ind];
	  double gtzl = gtz[ind];

	  double gxxl = gxx[ind];
	  double gxyl = gxy[ind];
	  double gxzl = gxz[ind];
	  double gyyl = gyy[ind];
	  double gyzl = gyz[ind];
	  double gzzl = gzz[ind];	 

	  // calculate the three-determinant
	  // -gxzl*gxzl*gyyl + 2.0*gxyl*gxzl*gyzl - gxxl*gyzl*gyzl - gxyl*gxyl*gzzl + gxxl*gyyl*gzzl;
	  //double Psil = pow(g,1./12.);

	  double g= pow(Psil,12.);
	  // Set the actual conformal factor
	  //psi[ind] = Psil;

	  // compute the gamma^ij
	  double gupxxl = (-gyzl*gyzl + gyyl*gzzl)/g;
	  double gupxyl =  (gxzl*gyzl - gxyl*gzzl)/g;
	  double gupxzl = (-gxzl*gyyl + gxyl*gyzl)/g;
	  double gupyyl = (-gxzl*gxzl + gxxl*gzzl)/g;
	  double gupyzl =  (gxyl*gxzl - gxxl*gyzl)/g;
	  double gupzzl = (-gxyl*gxyl + gxxl*gyyl)/g;

	  //shift
	  shiftx[ind] = gupxxl*gtx[ind] + gupxyl*gty[ind] + gupxzl*gtz[ind];
	  shifty[ind] = gupxyl*gtx[ind] + gupyyl*gty[ind] + gupyzl*gtz[ind];
	  shiftz[ind] = gupxzl*gtx[ind] + gupyzl*gty[ind] + gupzzl*gtz[ind];

	    
	  //Compute trK
	  trK[ind] = gupxxl*kxx[ind] + gupyyl*kyy[ind] + gupzzl*kzz[ind] +  
	    2.0*(gupxyl*kxy[ind] + gupxzl*kxz[ind] + gupyzl*kyz[ind]);


	  // compute its inverse
          double psiLm4=1.0/(Psil*Psil*Psil*Psil);
	  
	  // The 3 metric here is not conformally flat, so we set the BSSN metric, now
 	  gxx[ind] = gxxl*psiLm4;
	  gxy[ind] = gxyl*psiLm4;
	  gxz[ind] = gxzl*psiLm4;
	  gyy[ind] = gyyl*psiLm4;
	  gyz[ind] = gyzl*psiLm4;
	  gzz[ind] = gzzl*psiLm4;


 	  // Inverse conformally flat 3 metric 

          double psiL4=(Psil*Psil*Psil*Psil);
	    
	  gupxx[ind] = gupxxl*psiL4;
	  gupxy[ind] = gupxyl*psiL4;
	  gupxz[ind] = gupxzl*psiL4;
	  gupyy[ind] = gupyyl*psiL4;
	  gupyz[ind] = gupyzl*psiL4;
	  gupzz[ind] = gupzzl*psiL4;


	}
      }
    }

   
    cout << "Now compute the Lorentz factor and u_i" << endl;

    // Calculate the Lorentz factor and then u_low_i 
    for(int k=0;k<nz;k++){
      double z = zmin + k*dz;
      for(int j=0;j<ny;j++){
	double y = ymin + j*dy;
	for(int i=0;i<nx;i++){
	  double x = xmin + i*dx;

	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  double rhoo = rho_nbar[ind];
	  double lap  = lapse[ind];
	  double betx = shiftx[ind];
	  double bety = shifty[ind];
	  double betz = shiftz[ind];    
	  double Psil = psi[ind];

          double psil4=(Psil*Psil*Psil*Psil);

	  // full metric g_ij
	  double gxxl = gxx[ind]*psil4;
	  double gxyl = gxy[ind]*psil4;
	  double gxzl = gxz[ind]*psil4;
	  double gyyl = gyy[ind]*psil4;
	  double gyzl = gyz[ind]*psil4;
	  double gzzl = gzz[ind]*psil4;

	  double Uupx = ux[ind];
	  double Uupy = uy[ind];
	  double Uupz = uz[ind];


	  if (rhoo > 1.e-16){
	    double v2 = gxxl*Uupx*Uupx + gyyl*Uupy*Uupy  + gzzl*Uupz*Uupz + 2.*gxyl*Uupx*Uupy + 2.*gxzl*Uupx*Uupz  + 2.*gyzl*Uupy*Uupz; 
	    double Lore = 1./sqrt(1. - v2);
	    double Ulx = gxxl*Uupx + gxyl*Uupy + gxzl*Uupz;  
	    double Uly = gxyl*Uupx + gyyl*Uupy + gyzl*Uupz;  
	    double Ulz = gxzl*Uupx + gyzl*Uupy + gzzl*Uupz;

	    // Now set ux, uy and uz to the four-velocity with lower components
	    ux[ind] = Lore*Ulx;
	    uy[ind] = Lore*Uly;
	    uz[ind] = Lore*Ulz;
	  }
	  else{
	    ux[ind] = 0.;
	    uy[ind]= 0.;
	    uz[ind]= 0.;
	  }

	}
      }
    }

  }
  cout << " Done with read_inputfile_lowlevel_bhdisk_Cocal.C "<<endl;
}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_BH_disk)
  (const cGH **cctkGH,int *num_CO, int *piecewise,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *gxx,  double   *gxy,double *gxz,  double *gyy,  double *gyz,  double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
   double *psi, double *trK, double *rho_nbar,double *ux,double *uy, double *uz,
   int *nphase, double *abc, double *abi, double *rhoi, double *tedi)
{  
  read_inputfile_BH_disk
    (*cctkGH, *num_CO, *piecewise, *genID_cmdline_output_enable,*fisheye_enable,
			     *reset_shift_lapse,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
		             lapse,shiftx,shifty,shiftz,
		             gxx,gxy,gxz,gyy,gyz,gzz,
                             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
	         	     kxx,kxy,kxz,kyy,kyz,kzz,
                             psi,trK,rho_nbar,ux,uy, uz,
                             nphase, abc, abi, rhoi, tedi);
}


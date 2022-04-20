//-----------------------------------------------------------------------
// Read Cocal bin. files and do fancy things with them...
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

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_cocal)
  (const cGH **cctkGH,int *num_CO, int *piecewise, int *genID_cmdline_output_enable, int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   int *nphase, double *abc, double *abi, double *rhoi, double *tedi);
//   double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
//   double &kappa3, double &rhoo1, double &rhoo2);

extern "C" void read_inputfile_cocal(const cGH *cctkGH, 
                                     int num_CO, int piecewise, int genID_cmdline_output_enable,int fisheye_enable,
				     int reset_shift_lapse,double xmin,double ymin,double zmin,
				     double dx,double dy,double dz,int *ext,
				     double *lapse,double *shiftx,double *shifty,double *shiftz,
				     double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
				     double *rho_nbar,double *ux,double *uy, double *uz,
                                     int *nphase, double *abc, double *abi, double *rhoi, double *tedi)	   
//                                     double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
//                                     double &kappa3, double &rhoo1, double &rhoo2)
{
  static int nStatic;
  nStatic += 1;

  double xmax=xmin+ext[0]*dx;
  double ymax=ymin+ext[1]*dy;
  double zmax=zmin+ext[2]*dz;

    int filename_checksum = ext[0]*301+ext[1]*10+ext[2]*1000000  + 1241*CCTK_MyProc(cctkGH) + (int)(999999.*dx) + (int)(100100100.*dy) + (int)(800000.*dz) + 
      (int)(fabs(xmin)*10.)     + (int)(fabs(ymin)*101.)      + (int)(fabs(zmin)*1001.) +
      (int)(fabs(xmax)*999999.) + (int)(fabs(ymax)*30000.)    + (int)(fabs(zmax)*10100.)+
      (int)(xmax*xmax*191919. ) + (int)(ymax*ymax*30220. )    + (int)(zmax*zmax*20202.) +nStatic;    

    //  int filename_checksum = ext[1]*301+ext[2]*91221+ext[0] + 1241*CCTK_MyProc(cctkGH)+(int)(99999*dx)+ nStatic;
 


  char filename[100];
  char dumname[100];
  char inputstring[50];
  int nph ;
  nph = nphase[0];

  sprintf(filename,"CTS_bin-proc%d.d",filename_checksum);

  ofstream solfile;

  cout << "^^^^^^^^^^^^^^^^^^^^^^^ ENTER read_inputfile_lowlevel_Cocal.C ^^^^^^^^^^^^^^^^^^^^^^ "<<endl;

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
//    printf("./read_bin_ns %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d $ascii $order\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],filename_checksum);
  } 
  else 
  {
    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);

    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];
    double doubledum;
    int ntot,intdum, int_nx, int_ny, int_nz;
    ntot=nx*ny*nz;

    double gtt[ntot], gtx[ntot], gty[ntot], gtz[ntot];
    double gxx[ntot], gxy[ntot], gxz[ntot];
    double gyy[ntot], gyz[ntot], gzz[ntot], kzz_loc[ntot];


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
	  infile1.read((char *) &gtt[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
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
	  //NEXT GF IS A DUMMY.  We know Kzz from the other Kij's and the fact that trK==0!
	  infile1.read((char *) &kzz_loc[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));

	  infile1.read((char *) &ux[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &uy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  infile1.read((char *) &uz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ READING PEOS FILE ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "<<endl ;
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

//    infile1.read((char *) &gamma1, sizeof(double));
//    gamma1 = abi[0];  
//    cout<<"gamma1="<<gamma1<<endl;

//    infile1.read((char *) &gamma2, sizeof(double));
//    gamma2 = abi[1];  
//    cout<<"gamma2="<<gamma2<<endl;

//    infile1.read((char *) &gamma3, sizeof(double));
//    gamma3 = abi[2];  
//    cout<<"gamma3="<<gamma3<<endl;

//    infile1.read((char *) &kappa1, sizeof(double));
//    kappa1 = abc[0];  
//    cout<<"kappa1="<<kappa1<<endl;

//    infile1.read((char *) &kappa2, sizeof(double));
//    kappa2 = abc[1];  
//    cout<<"kappa2="<<kappa2<<endl;

//    infile1.read((char *) &kappa3, sizeof(double));
//    kappa3 = abc[2];  
//    cout<<"kappa3="<<kappa3<<endl;

//    infile1.read((char *) &rhoo1, sizeof(double));
//    rhoo1=rhoi[1];                           // surface
//    cout<<"rhoo1="<<rhoo1<<endl;

//    infile1.read((char *) &rhoo2, sizeof(double));
//    rhoo2=rhoi[2];
//    cout<<"rhoo2="<<rhoo2<<endl;

    infile1.close();  

    // Now calculate the conformal factor

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double g_xx = gxx[ind];
	  double g_xy = gxy[ind];
	  double g_xz = gxz[ind];
	  double g_yy = gyy[ind];
	  double g_yz = gyz[ind];
	  double g_zz = gzz[ind];
	  
	  // calculate the three-determinant
	  
	  double g= -g_xz*g_xz*g_yy + 2*g_xy*g_xz*g_yz - g_xx*g_yz*g_yz - g_xy*g_xy*g_zz + g_xx*g_yy*g_zz;
	  psi[ind] = pow(g,1./12.);
	}
      }
    }

    // Check that the metric is conformally flat
    
    double eps = 1.e-12;
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double psi4 = psi[ind]*psi[ind]*psi[ind]*psi[ind];
	  double g_xx = gxx[ind]/psi4;
	  double g_xy = gxy[ind]/psi4;
	  double g_xz = gxz[ind]/psi4;
	  double g_yy = gyy[ind]/psi4;
	  double g_yz = gyz[ind]/psi4;
	  double g_zz = gzz[ind]/psi4;
	  

	  if ( fabs(g_xx-1.) > eps && fabs(g_xy-1.) > eps && fabs(g_xz-1.) > eps && fabs(g_yy-1.) > eps && fabs(g_yz-1.) > eps && fabs(g_zz-1.) > eps){
	    printf("Stopping: metric is not conformally flat\n");
	    exit(1);
	  }

	}
      }
    }

    printf("Cocal: done with metric\n");

    // Check that tr K = 0
    
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  
	  double k_xx = kxx[ind];
	  double k_yy = kyy[ind];
	  double k_zz = kzz_loc[ind];
	  
	  double traceK = k_xx+k_yy+k_zz;

	  if ( fabs(traceK) > eps){
	    printf("Stopping: slicing is not maximal\n");
	    exit(1);
	  }

	}
      }
    }

    printf("Metric is conformally flat, as expected\n");
    printf("trK=0, as expected\n");


    // Now calculate the  gauge variables

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  
	  double g_tt = gtt[ind];
	  double g_tx = gtx[ind];
	  double g_ty = gty[ind];
	  double g_tz = gtz[ind];
	  double g_xx = gxx[ind];
	  double g_xy = gxy[ind];
	  double g_xz = gxz[ind];
	  double g_yy = gyy[ind];
	  double g_yz = gyz[ind];
	  double g_zz = gzz[ind];
	  
	  // det of full metric
	  double detab =  g_tz*g_tz*g_xy*g_xy - 2.*g_ty*g_tz*g_xy*g_xz
	    + g_ty*g_ty*g_xz*g_xz - g_tz*g_tz*g_xx*g_yy
	    + 2.*g_tx*g_tz*g_xz*g_yy - g_tt*g_xz*g_xz*g_yy
	    + 2.*g_ty*g_tz*g_xx*g_yz -2.*g_tx*g_tz*g_xy*g_yz
	    - 2.*g_tx*g_ty*g_xz*g_yz + 2.*g_tt*g_xy*g_xz*g_yz
	    + g_tx*g_tx*g_yz*g_yz - g_tt*g_xx*g_yz*g_yz
	    - g_ty*g_ty*g_xx*g_zz + 2.*g_tx*g_ty*g_xy*g_zz
	    -g_tt*g_xy*g_xy*g_zz - g_tx*g_tx*g_yy*g_zz 
	    + g_tt*g_xx*g_yy*g_zz;
	  double detabinv = 1.0/detab;	  

	  // lapse
	  double ginv_tt = (-g_xz*g_xz*g_yy + 2.*g_xy*g_xz*g_yz - g_xx*g_yz*g_yz - g_xy*g_xy*g_zz + g_xx*g_yy*g_zz) * detabinv;

	  lapse[ind] = 1.0/sqrt(-ginv_tt);

	  //shift
	  double ginv_tx = (g_tz*g_xz*g_yy - g_tz*g_xy*g_yz - g_ty*g_xz*g_yz + g_tx*g_yz*g_yz + g_ty*g_xy*g_zz - g_tx*g_yy*g_zz) * detabinv;
	  double ginv_ty = (-(g_tz*g_xy*g_xz) + g_ty*g_xz*g_xz + g_tz*g_xx*g_yz - g_tx*g_xz*g_yz - g_ty*g_xx*g_zz + g_tx*g_xy*g_zz) * detabinv;
	  double ginv_tz = (g_tz*g_xy*g_xy - g_ty*g_xy*g_xz - g_tz*g_xx*g_yy + g_tx*g_xz*g_yy + g_ty*g_xx*g_yz - g_tx*g_xy*g_yz) * detabinv;


	  shiftx[ind] = ginv_tx*lapse[ind]*lapse[ind];
	  shifty[ind] = ginv_ty*lapse[ind]*lapse[ind];
	  shiftz[ind] = ginv_tz*lapse[ind]*lapse[ind];

	}
      }
    }

    printf("Cocal: done with Gauge variables\n");


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


	  double g_xx = gxx[ind];
	  double g_xy = gxy[ind];
	  double g_xz = gxz[ind];
	  double g_yy = gyy[ind];
	  double g_yz = gyz[ind];
	  double g_zz = gzz[ind];

	  double Uupx = ux[ind];
	  double Uupy = uy[ind];
	  double Uupz = uz[ind];


	  if (rhoo > 1.e-16){
	    double v2 = g_xx*Uupx*Uupx + 2.*g_xy*Uupx*Uupy + 2.*g_xz*Uupx*Uupz + g_yy*Uupy*Uupy + 2.*g_yz*Uupz*Uupy + g_zz*Uupz*Uupz; 
	    double Lore = 1./sqrt(1.-v2);
	    double Ulx = g_xx*Uupx + g_xy*Uupy + g_xz*Uupz;  
	    double Uly = g_xy*Uupx + g_yy*Uupy + g_yz*Uupz;  
	    double Ulz = g_xz*Uupx + g_yz*Uupy + g_zz*Uupz;

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
  cout << "^^^^^^^^^^^^^^^^^^^^^^^ EXIT read_inputfile_lowlevel_Cocal.C ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "<<endl;
}

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_cocal)
  (const cGH **cctkGH, int *num_CO, int *piecewise, int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   int *nphase, double *abc, double *abi, double *rhoi, double *tedi)
//   double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
//   double &kappa3, double &rhoo1, double &rhoo2)
{  
  read_inputfile_cocal(*cctkGH, *num_CO, *piecewise, *genID_cmdline_output_enable,*fisheye_enable,
			     *reset_shift_lapse,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
			     lapse,shiftx,shifty,shiftz,
			     psi,kxx,kxy,kxz,kyy,kyz,
			     rho_nbar,ux,uy, uz,
                             nphase, abc, abi, rhoi, tedi);
//			     gamma1,gamma2,gamma3,kappa1,kappa2, kappa3,rhoo1,rhoo2);
}


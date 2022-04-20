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

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_wdns_lorene)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
   double &kappa3, double &rhoo1, double &rhoo2);

extern "C" void read_inputfile_wdns_lorene(const cGH *cctkGH,int genID_cmdline_output_enable,int fisheye_enable,
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
  char dumname[100];
  double xrot, Omega;
  sprintf(filename,"CTS_bin-proc%d.d",filename_checksum);

  ofstream solfile;




  if(genID_cmdline_output_enable==1) {
    printf("./wdns_initialdata_isotropic %.16e %.16e %.16e %.16e %.16e %.16e %d %d %d %d $ascii $order\n",xmin,ymin,zmin,dx,dy,dz,ext[0],ext[1],ext[2],filename_checksum);
  } else {

    printf("Attempting to read %s now...  \nNote that you'll need to store the initial data files so that ALL processors can see them!\n",filename);


    int nx = ext[0];
    int ny = ext[1];
    int nz = ext[2];
    double doubledum;
    int ntot,intdum, int_nx, int_ny, int_nz;
    ntot=nx*ny*nz;
    double xco[ntot], yco[ntot], zco[ntot];
    double gxx[ntot], gxy[ntot], gxz[ntot];
    double gyy[ntot], gyz[ntot], gzz[ntot];


    ifstream infile1;
    infile1.open(filename);
    if(!infile1) {
      cerr << "\a Can't open " << filename << " for input." << endl;
      exit(1);
    }

    // Read grid parameters
    infile1.read((char *) &dumname, sizeof(char));
    cout<<"WDNS parameters read"<<endl;
    cout<< " EOS name: " << dumname <<endl;
    infile1.read((char *) &gamma1, sizeof(double));
    infile1.read((char *) &kappa1, sizeof(double));
    infile1.read((char *) &dumname, sizeof(char));
    infile1.read((char *) &gamma2, sizeof(double));
    infile1.read((char *) &kappa2, sizeof(double));
    infile1.read((char *) &Omega, sizeof(double));
    cout<<"Omega="<<Omega<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"distance="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"distance between c.o.m's="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"M_baryon_NS="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"M_baryon_WD="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"ADM mass="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"Ang mom="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"x radius NS="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    cout<<"x radius WD="<<doubledum<<endl;
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &doubledum, sizeof(double));
    infile1.read((char *) &ntot, sizeof(int));  


    if(ntot!=nx*ny*nz) {
      printf("MAYBE PROBLEM: This data file has %d elements, but current grid has %d!\n",ntot,nx*ny*nz);
      exit(1);
    } else {
      printf("This data file has %d elements, like we expect!\n",nx*ny*nz);
    }

    //    if(reset_shift_lapse==0) {

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &xco[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  double x=xmin+i*dx;
	  if(fabs(xco[CCTK_GFINDEX3D(cctkGH,i,j,k)]-x)/x > 1.e-10) {
	    printf("PROBLEM: This data file has been read in an order inconsistent with our code\n");
	    printf("Coordinate is x = %15.12f, slot read in corresponds to x = %15.12f\n",x,xco[CCTK_GFINDEX3D(cctkGH,i,j,k)]); 
	    exit(1);
	  }
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  double y=ymin+j*dx;
	  infile1.read((char *) &yco[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  if(fabs(yco[CCTK_GFINDEX3D(cctkGH,i,j,k)]-y)/y > 1.e-10) {
	    printf("PROBLEM: This data file has been read in an order inconsistent with our code\n");
	    printf("Coordinate is y = %15.12f, slot read in corresponds to y = %15.12f\n",y,yco[CCTK_GFINDEX3D(cctkGH,i,j,k)]); 
	    exit(1);
	  }
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  double z=zmin+k*dx;
	  infile1.read((char *) &zco[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	  if(fabs(zco[CCTK_GFINDEX3D(cctkGH,i,j,k)]-z)/z > 1.e-10) {
	    printf("PROBLEM: This data file has been read in an order inconsistent with our code\n");
	    printf("Coordinate is z = %15.12f, slot read in corresponds to z = %15.12f\n",z,zco[CCTK_GFINDEX3D(cctkGH,i,j,k)]); 
	    exit(1);
	  }
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &shiftx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &shifty[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &shiftz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &gxx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &gxy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &gxz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &gyy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &gyz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &gzz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &kxx[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &kxy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &kxz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &kyy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &kyz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    //NEXT GF IS A DUMMY.  We know Kzz from the other Kij's and the fact that trK==0!
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    // Rest-mass density in LORENE data files is given in kg/m^3. Here we convert the rest-mass density to units of 1/km^2;
  
    double Goc2 = 7.42471382e-28; // G/c^2 [m/kg] multiplying by this factor we get the density in 1/m^2

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k); 
	  rho_nbar[ind] = rho_nbar[ind]*Goc2*10000.;
	}
      }
    }


    // Next GF is a dummy we will determine the conformal factor later
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &psi[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));    
	}
      }
    }


    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &ux[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &uy[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  infile1.read((char *) &uz[CCTK_GFINDEX3D(cctkGH,i,j,k)],sizeof(double));
	}
      }
    }

      //    }

    infile1.close();  


    gamma2 = gamma1;
    gamma3 = gamma1;
    kappa2 = kappa1;
    kappa3 = kappa1;
    rhoo1 = 1.;
    rhoo2 = 1.e-1;


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
    
    // Calculate the Lorentz factor and then u_low_i 

    //    solfile.open("lap.dat");
    for(int k=0;k<nz;k++){
      double z = zmin + k*dz;
      for(int j=0;j<ny;j++){
	double y = ymin + j*dy;
	for(int i=0;i<nx;i++){
	  double x = xmin + i*dx;


	  double rhoo = rho_nbar[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double lap = lapse[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double betx = shiftx[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double bety = shifty[CCTK_GFINDEX3D(cctkGH,i,j,k)];
	  double betz = shiftz[CCTK_GFINDEX3D(cctkGH,i,j,k)];    

	  int ind = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double g_xx = gxx[ind];
	  double g_xy = gxy[ind];
	  double g_xz = gxz[ind];
	  double g_yy = gyy[ind];
	  double g_yz = gyz[ind];
	  double g_zz = gzz[ind];
	  double Uupx = ux[ind];
	  double Uupy = uy[ind];
	  double Uupz = uy[ind];

	  //	  solfile << x <<" " << y  << " " << y << " " << lap << endl;

	  if (rhoo > 1.e-16){
	    double v2 = g_xx*Uupx*Uupx + 2.*g_xy*Uupx*Uupy + 2.*g_xz*Uupx*Uupz + g_yy*Uupy*Uupy + 2.*g_yz*Uupz*Uupy + g_zz*Uupz*Uupz; 
	    double Lore = 1./sqrt(1.-v2);
	    double Ulx = g_xx*Uupx + g_xy*Uupy + g_xz*Uupz;  
	    double Uly = g_xy*Uupx + g_yy*Uupy + g_yz*Uupz;  
	    double Ulz = g_xz*Uupx + g_yz*Uupy + g_zz*Uupz;
	    // Now set ux, uy and uz to the four-velocity with lower components
	    ux[CCTK_GFINDEX3D(cctkGH,i,j,k)] = Lore*Ulx;
	    uy[CCTK_GFINDEX3D(cctkGH,i,j,k)]= Lore*Uly;
	    uz[CCTK_GFINDEX3D(cctkGH,i,j,k)]= Lore*Ulz;
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

extern "C" void CCTK_FCALL CCTK_FNAME(read_inputfile_wdns_lorene)
  (const cGH **cctkGH,int *genID_cmdline_output_enable,int *fisheye_enable,
   int *reset_shift_lapse,double *xmin,double *ymin,double *zmin,
   double *dx,double *dy,double *dz,int *ext,
   double *lapse,double *shiftx,double *shifty,double *shiftz,
   double *psi,double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,
   double *rho_nbar,double *ux,double *uy, double *uz,
   double &gamma1,double &gamma2, double &gamma3, double &kappa1, double &kappa2, 
   double &kappa3, double &rhoo1, double &rhoo2)
{  
  read_inputfile_wdns_lorene(*cctkGH,*genID_cmdline_output_enable,*fisheye_enable,
			     *reset_shift_lapse,*xmin,*ymin,*zmin,*dx,*dy,*dz,ext,
			     lapse,shiftx,shifty,shiftz,
			     psi,kxx,kxy,kxz,kyy,kyz,
			     rho_nbar,ux,uy, uz,
			     gamma1,gamma2,gamma3,kappa1,kappa2, kappa3,rhoo1,rhoo2);
}


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
#include "cctk_Functions.h"

#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
//#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <fstream.h>
//#include "nrutil.h"

#include "Symmetry.h"

void Read_irr_Input(cGH *cctkGH,const char * infile, int I_Print,int *nx,int *ny,int *nz,int *nstar,
		    double *dx,double *dy,double *dz, 
		    double *xbd,double *ybd,double *zbd,double *rradi,
		    double *rkappa,double *pinx ,double *romega,double *radi);

extern void Setup_Irrot_Data(cGH *cctkGH,double *PhysicalRadius,
			     double *lapse,double *phi,double *shiftx,double *shifty,
			     double *shiftz, double *Axx, double *Axy, double *Axz,
                             double *Ayy, double *Ayz, double *Azz, 
			     double *q,double *ux,double *uy,double *uz,double *ut,
			     double *rho,double *S,double *Sx,double *Sy,double *Sz, 
			     double *Sxx,double *Sxy,double *Sxz,double *Syy,
			     double *Syz,double *Szz,
			     double dx_adm,double dY_adm,double dZ_adm,
			     int & nx, int & ny, int & nz, int *ext,
			     double *xdest,double *ydest,double *zdest,
			     int Symmetry,int new_interpolation_choice, int kt_id, 
			     double x_max, double y_max, double z_max,
                             double x_min, double y_min, double z_min);
//int *lb, int *ub);
/*
extern void ReadtoArray(const char * filename, double ** data,
			double & xmin, double & xmax, double & ymin, 
			double & ymax, double & zmin, double & zmax, 
			int & nx, int & ny, int & nz);

void Read_irr_Input(cGH *cctkGH,char * infile, int I_Print,int *nx,int *ny,int *nz,int *nstar,
		    double *dx,double *dy,double *dz, 
		    double *xbd,double *ybd,double *zbd,double *rradi,
		    double *rkappa,double *pinx ,double *romega,double *radi);

void setup_irr(double dx_adm, double dy_adm, double dz_adm,
	       double dx_ref, double dy_ref, double dz_ref, 
	       double *lapse_s, double *phi_s,
	       double *shiftx_s, double *shifty_s, double *shiftz_s,
	       double *q_s,
	       double *ux_s, double *uy_s, double *uz_s,
	       int nx_adm, int ny_adm, int nz_adm,
	       int Nxh, int Nyh, int Nzh,
	       double *lapse_i, double *phi_i, double *shiftx_i, double *shifty_i, double *shiftz_i,
	       double *q_i, double *ux_i, double *uy_i, double *uz_i);

void push(int Nx, int Ny, int Nz, 
	  double symx, double symy, double symz, 
	  double *array);

void adjust_z_grid(int Nx, int Ny, int Nz, int *Nz_t, int nz_adm,
		   double dz_adm, double dz_ref, double symz, 
		   double *array_s, double ** array_t);

double interp_3d(double *f, int il, int jl, int kl,
		 int Nx, int Ny, int Nz_t,
		 double Dx, double Dy, double Dz, int skin);
void interp4(int ntot,int Nx,int Ny,int nx_adm,int ny_adm,double *f,double &f_int,int il,
	     int iu,int jl,int ju,int kl,double Dx,double Dy,double SymX,double SymY);

void translate_irr(cGH *cctkGH,int *ex,int Nx,int Ny,int Nz,int nx_adm,int ny_adm,int ntot,
		   double *Xsource,double *Ysource,double *Zsource,
		   double *Xdest, double *Ydest,double *Zdest,
		   double *shiftx_old,double *shifty_old,double *shiftz_old,double *shiftx_new,
		   double * shifty_new,double *shiftz_new,
		   double *v2_old,double *v2_new,double *lapse_old,double *lapse_new,
		   double * phi_old,double *phi_new,
		   double *ux_old,double *uy_old,double *uz_old,
		   double *ux_new,double *uy_new,double *uz_new);
*/

//extern void Setup_Hybrid_EOS(double K_poly, double n, char * eos_file,int neos,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab);
void Setup_Hybrid_EOS_old(double K_poly, double n, char * eos_file,int neos,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, const double * gama_th);

static char *rcsid = "$Meow...$";
CCTK_FILEVERSION(BNS_initialdata_readfiles)

  extern "C" void BNS_initialdata_readfiles(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  //======================================
  // Read Input file for cbwh_adm
  //======================================
  int setphi = 1;
  double n = 1.0;
  /*
    int nx = cctk_lsh[0];
    int ny = cctk_lsh[1];
    int nz = cctk_lsh[2];
  */
  int ext[3];
  ext[0] = cctk_lsh[0];
  ext[1] = cctk_lsh[1];
  ext[2] = cctk_lsh[2];

  Setup_Hybrid_EOS_old(K_poly, n, "eos.tab",neos,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, &gamma_th);
  
  int nx=0,ny=0,nz=0,nstar=0;
  double dx=0,dy=0,dz=0;
  double xbd=0,ybd=0,zbd=0,rradi=0,rkappa=0,pinx =0,romega=0,radi=0;


  if(CCTK_Equals(rot_type,"corotational")) {
    printf("corotational initial data not yet supported!\n");
    CCTK_Abort(cctkGH, EXIT_FAILURE);
    exit(1);
  }
  else {
     Read_irr_Input(cctkGH,input_filename,CCTK_MyProc(cctkGH),&nx,&ny,&nz,&nstar,
		     &dx,&dy,&dz, 
		     &xbd,&ybd,&zbd,&rradi,
		     &rkappa,&pinx ,&romega,&radi);
/*
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&nx,&nx,CCTK_VARIABLE_INT);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&ny,&ny,CCTK_VARIABLE_INT);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&nz,&nz,CCTK_VARIABLE_INT);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&nstar,&nstar,CCTK_VARIABLE_INT);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&dx,&dx,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&dy,&dy,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&dz,&dz,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&xbd,&xbd,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&ybd,&ybd,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&zbd,&zbd,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&rradi,&rradi,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&rkappa,&rkappa,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&pinx,&pinx,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&romega,&romega,CCTK_VARIABLE_REAL);
    CCTK_ReduceLocScalar (cctkGH,-1,CCTK_ReductionHandle("sum"),&radi,&radi,CCTK_VARIABLE_REAL);
*/
  }

  //===================================================================
  // Setup the Initial Data
  //===================================================================
  if(CCTK_Equals(rot_type,"corotational")) {
  }
  else if(CCTK_Equals(rot_type,"irrotational")) {
    double dx_adm = cctk_delta_space[0];
    double dy_adm = cctk_delta_space[1];
    double dz_adm = cctk_delta_space[2];
    double x_max,x_min,y_max,y_min,z_max,z_min;

    int ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&x_max,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&y_max,1,CCTK_VarIndex("grid::Y"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&z_max,1,CCTK_VarIndex("grid::Z"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&x_min,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&y_min,1,CCTK_VarIndex("grid::Y"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&z_min,1,CCTK_VarIndex("grid::Z"));

    if (Symmetry==1) {
        // Equatorial symmetry
       z_min = z_min + cctk_nghostzones[2]*dz_adm;
    } else if (Symmetry !=0) {
      cerr << "Symmetry option not supported" << endl;
      CCTK_Abort(cctkGH, EXIT_FAILURE);
      exit(1);
    }
    //f_Test();
    int new_interpolation_choice=1; //Use new 3d interpolation scheme, should work on multiple procs
    int Nx = (int)(floor( (x_max-x_min)/dx_adm + 1.e-5 )) + 1;
    int Ny = (int)(floor( (y_max-y_min)/dy_adm + 1.e-5 )) + 1;
    int Nz = (int)(floor( (z_max-z_min)/dz_adm + 1.e-5 )) + 1;

    Setup_Irrot_Data(cctkGH,PhysicalRadius,lapm1,phi,shiftx,shifty,shiftz,
		     Axx,Axy,Axz,Ayy,Ayz,Azz,rho_br,vxr,vyr,vzr,Pr,rho,S,
		     Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,
		     dx_adm,dy_adm,dz_adm,Nx,Ny,Nz,ext,
		     x,y,z,Symmetry,new_interpolation_choice,kt_id, 
		     x_max, y_max, z_max, x_min, y_min, z_min);
  }

}


/*
  extern void ReadtoArray(char *lapse_fname, double **lapse_s,
  double xmin,double xmax,double ymin,double ymax,double zmin,double zmax, 
  int &Nxh,int &Nyh,int &Nzh);
*/


void Setup_Hybrid_EOS_old(double K_poly, double n, char * eos_file,int neos,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, 
const double * gamma_th) {
  rho_tab[0] = 1000.0;
  P_tab[0] = K_poly * pow(rho_tab[0],1.0+1.0/n);
  eps_tab[0] = n*P_tab[0]/rho_tab[0]; 
  gamma_tab[0] = 1.0 + 1.0/n; k_tab[0] = K_poly;
  gamma_tab[1] = 1.0 + 1.0/n; k_tab[1] = K_poly;
  // *gamma_th = 1.0 + 1.0/n;
}

void Read_irr_Input(cGH *cctkGH,const char * infile, int I_Print,int *nx,int *ny,int *nz,int *nstar,
		    double *dx1,double *dy1,double *dz1, 
		    double *xbd,double *ybd,double *zbd,double *rradi,
		    double *rkappa,double *pinx, double *romega,double *radi)
{
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS
  /*
    ifstream file0;
    file0.open(infile);
    if (!file0) {
    cerr << "Can't open " << infile << " for input." << endl;
    CCTK_Abort(cctkGH, EXIT_FAILURE);
    exit(1);
    } 
    char buf[600],c;
    file0.get(buf,600,'='); file0.get(c); file0 >> rho_fact;
    file0.get(buf,600,'='); file0.get(c); file0 >> EM_evolve;
    if(file0.eof()) {
    cerr << "Error reading input file " << infile << endl;
    CCTK_Abort(cctkGH, EXIT_FAILURE);
    exit(1);
    }
    file0.close();
  */

  ifstream file;
  file.open("params");
  if (!file) {
    cerr << "Can't open params for input." << endl;
    CCTK_Abort(cctkGH, EXIT_FAILURE);
    exit(1);
  }
  
  if (kt_id==1) {
    // These are for Keisuke's initial data. Check to see if the input grid
    // agrees with the initial data grid, and if the data files are here.
    double x_max,x_min,y_max,y_min,z_max,z_min;
    int ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&x_max,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&y_max,1,CCTK_VarIndex("grid::Y"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&z_max,1,CCTK_VarIndex("grid::Z"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&x_min,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&y_min,1,CCTK_VarIndex("grid::Y"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&z_min,1,CCTK_VarIndex("grid::Z"));
    double x_min_id, x_max_id, y_min_id, y_max_id, z_min_id, z_max_id;
    double dx_id, dy_id, dz_id;
    double dX = cctk_delta_space[0];
    double dY = cctk_delta_space[1];
    double dZ = cctk_delta_space[2];
    int fish_id;
    if (Symmetry==1) {
        // Equatorial symmetry
       z_min = z_min + cctk_nghostzones[2]*dZ;
    } else if (Symmetry !=0) {
      cerr << "Symmetry option not supported" << endl;
      CCTK_Abort(cctkGH, EXIT_FAILURE);
      exit(1);
    }
    double tiny = dX*1.e-5;
    x_min_id = -999.9; x_max_id = -999.9;
    y_min_id = -999.9; y_max_id = -999.9;
    z_min_id = -999.9; z_max_id = -999.9;
    dx_id = -999.9; dy_id= -999.9; dz_id = -999.9;
    file >> x_min_id >> x_max_id;
    file >> y_min_id >> y_max_id;
    file >> z_min_id >> z_max_id;
    file >> dx_id;
    file >> dy_id;
    file >> dz_id;
    file >> fish_id;
    int exit_flag = 0;
    if (fabs(x_min_id - x_min) > tiny) {
       exit_flag = 1;
       cerr << "xmin of the current grid = " << x_min << endl;
       cerr << "xmin of the initial data grid = " << x_min_id << endl;
       cerr << "xmin doesn't agree with the initial data file!" << endl;
       x_min = x_min_id;
    }
    if (fabs(x_max_id - x_max) > tiny) {
       exit_flag = 1;
       cerr << "xmax of the current grid = " << x_max << endl;
       cerr << "xmax of the initial data grid = " << x_max_id << endl;
       cerr << "xmax doesn't agree with the initial data file!" << endl;
       x_max = x_max_id;
    }
    if (fabs(y_min_id - y_min) > tiny) {
       exit_flag = 1;
       cerr << "ymin of the current grid = " << y_min << endl;
       cerr << "ymin of the initial data grid = " << y_min_id << endl;
       cerr << "ymin doesn't agree with the initial data file!" << endl;
       y_min = y_min_id;
    }
    if (fabs(y_max_id - y_max) > tiny) {
       exit_flag = 1;
       cerr << "ymax of the current grid = " << y_max << endl;
       cerr << "ymax of the initial data grid = " << y_max_id << endl;
       cerr << "ymax doesn't agree with the initial data file!" << endl;
       y_max = y_max_id;
    }
    if (fabs(z_min_id - z_min) > tiny) {
       exit_flag = 1;
       cerr << "zmin of the current grid = " << z_min << endl;
       cerr << "zmin of the initial data grid = " << z_min_id << endl;
       cerr << "zmin doesn't agree with the initial data file!" << endl;
       z_min = z_min_id;
    }
    if (fabs(z_max_id - z_max) > tiny) {
       exit_flag = 1;
       cerr << "zmax of the current grid = " << z_max << endl;
       cerr << "zmax of the initial data grid = " << z_max_id << endl;
       cerr << "zmax doesn't agree with the initial data file!" << endl;
       z_max = z_max_id;
    }
    if (fabs(dx_id - dX) > tiny) {
       exit_flag = 1;
       cerr << "dX of the current grid = " << dX << endl;
       cerr << "dX of the initial data grid = " << dx_id << endl;
       cerr << "dX doesn't agree with the initial data file!" << endl;
    }
    if (fabs(dy_id - dY) > tiny) {
       exit_flag = 1;
       cerr << "dY of the current grid = " << dY << endl;
       cerr << "dY of the initial data grid = " << dy_id << endl;
       cerr << "dY doesn't agree with the initial data file!" << endl;
    }
    if (fabs(dz_id - dZ) > tiny) {
       exit_flag = 1;
       cerr << "dZ of the current grid = " << dZ << endl;
       cerr << "dZ of the initial data grid = " << dz_id << endl;
       cerr << "dZ doesn't agree with the initial data file!" << endl;
    }

    // Check if both the current run and the initial data are using fisheye
    if (fisheye_enable != fish_id) {
       exit_flag = 1;
       cerr << "Fisheye flag of the current run = " << fisheye_enable << endl;
       cerr << "Fisheye flag of the initial data = " << fish_id << endl;
       cerr << "Fisheye flag doesn't agree with the initial data file!" << endl;
    } 

    // Now check if the fisheye parameters agree
    if (fisheye_enable==1 && fish_id==1) {
       ifstream fisheyefile;
       fisheyefile.open("RadialCoordinate_Input");
       char buf[100],c;
       int nFE1,nFE2; 
       double tmp1,tmp2;
       fisheyefile.get(buf,100,'='); fisheyefile.get(c); fisheyefile >> nFE1; 
       file.get(buf,100,'='); file.get(c); file >> nFE2;
       if (nFE1 != nFE2) {
	  exit_flag = 1;
	  cerr << "Number of fisheye transitions of the current run = " << 
			nFE1 << endl;
          cerr << "Number of fisheye transitions of the initial data = " << 
			nFE2 << endl;
	  cerr << "Number of fisheye transitions doesn't agree with the initial data file!" << endl;
       } else {
	  // Check the fisheye "a" parameters
          file.get(buf,100,'='); file.get(c); file >> tmp2;
	  fisheyefile.get(buf,100,'='); fisheyefile.get(c); fisheyefile >> tmp1;
          if (fabs((tmp1-tmp2)/tmp1) > 1.e-10) {
	     exit_flag = 1;
	     cerr << "One of the fisheye parameter disagrees with the initial data" << endl;
	  }
	  if (nFE1 > 0) { 
	     for(int i=1; i<=nFE1; i++) { 
		file >> tmp1; fisheyefile >> tmp2;
                if (fabs((tmp1-tmp2)/tmp1) > 1.e-10) {
                   exit_flag = 1;
                   cerr << "One of the fisheye parameter disagrees with the initial data" << endl;
	        }
	     } // end of for(int i=1; i<=nFE1; i++)

	     // Check the fisheye "r0" parameters 
	     file.get(buf,100,'='); file.get(c); file >> tmp2;
 	     fisheyefile.get(buf,100,'='); fisheyefile.get(c); fisheyefile >> tmp1;
             if (fabs((tmp1-tmp2)/tmp1) > 1.e-10) {
                exit_flag = 1;
                cerr << "One of the fisheye parameter disagrees with the initial data" << endl;
	     }
	     for(int i=2; i<=nFE1; i++) { 
		file >> tmp2; fisheyefile >> tmp1;
                if (fabs((tmp1-tmp2)/tmp1) > 1.e-10) {
                   exit_flag = 1;
                   cerr << "One of the fisheye parameter disagrees with the initial data" << endl;
                }
             } // end of for(int i=2; i<=nFE1; i++)

	     // Check the fisheye "s" parameters
	     file.get(buf,100,'='); file.get(c); file >> tmp2;
	     fisheyefile.get(buf,100,'='); fisheyefile.get(c); fisheyefile >> tmp1;
	     if (fabs((tmp1-tmp2)/tmp1) > 1.e-10) {
                exit_flag = 1;
                cerr << "One of the fisheye parameter disagrees with the initial data" << endl;
             }
             for(int i=2; i<=nFE1; i++) {
                file >> tmp2; fisheyefile >> tmp1;
                if (fabs((tmp1-tmp2)/tmp1) > 1.e-10) {
                   exit_flag = 1;
                   cerr << "One of the fisheye parameter disagrees with the initial data" << endl;
                }
             } // end of for(int i=2; i<=nFE1; i++)

          } // end of if (nFE1 > 0)
       } // end of if (nFE1 != nFE2)
       fisheyefile.close();
    } // end of if (fisheye_enable==1 && fish_id==1)

    ifstream Psiidfile,lapseidfile,rhobidfile,uxidfile,uyidfile,uzidfile;
    ifstream sxidfile,syidfile,szidfile,kxxidfile,kxyidfile,kxzidfile;
    ifstream kyyidfile,kyzidfile,kzzidfile;
    Psiidfile.open("Psi.id");
    if (!Psiidfile) {
       exit_flag = 1;
       cerr << "File Psi.id is missing!" << endl;
    }
    Psiidfile.close();
    lapseidfile.open("lapse.id");
    if (!lapseidfile) {
       exit_flag = 1;
       cerr << "File lapse.id is missing!" << endl;
    }
    lapseidfile.close();
    rhobidfile.open("rho_b.id");
    if (!rhobidfile) {
       exit_flag = 1;
       cerr << "File rho_b.id is missing!" << endl;
    }
    rhobidfile.close();
    uxidfile.open("u_x.id");
    if (!uxidfile) {
       exit_flag = 1;
       cerr << "File u_x.id is missing!" << endl;
    }
    uxidfile.close();
    uyidfile.open("u_y.id");
    if (!uyidfile) {
       exit_flag = 1;
       cerr << "File u_y.id is missing!" << endl;
    }
    uyidfile.close();
    uzidfile.open("u_z.id");
    if (!uzidfile) {
       exit_flag = 1;
       cerr << "File u_z.id is missing!" << endl;
    }
    uzidfile.close();
    sxidfile.open("shiftx.id");
    if (!sxidfile) {
       exit_flag = 1;
       cerr << "File shiftx.id is missing!" << endl;
    }
    sxidfile.close();
    syidfile.open("shifty.id");
    if (!syidfile) {
       exit_flag = 1;
       cerr << "File shifty.id is missing!" << endl;
    }
    syidfile.close();
    szidfile.open("shiftz.id");
    if (!szidfile) {
       exit_flag = 1;
       cerr << "File shiftz.id is missing!" << endl;
    }
    szidfile.close();
    kxxidfile.open("Kxx.id");
    if (!kxxidfile) {
       exit_flag = 1;
       cerr << "File Kxx.id is missing!" << endl;
    }
    kxxidfile.close();
    kxyidfile.open("Kxy.id");
    if (!kxyidfile) {
       exit_flag = 1;
       cerr << "File Kxy.id is missing!" << endl;
    }
    kxyidfile.close();
    kxzidfile.open("Kxz.id");
    if (!kxzidfile) {
       exit_flag = 1;
       cerr << "File Kxz.id is missing!" << endl;
    }
    kxzidfile.close();
    kyyidfile.open("Kyy.id");
    if (!kyyidfile) {
       exit_flag = 1;
       cerr << "File Kyy.id is missing!" << endl;
    }
    kyyidfile.close();
    kyzidfile.open("Kyz.id");
    if (!kyzidfile) {
       exit_flag = 1;
       cerr << "File Kyz.id is missing!" << endl;
    }
    kyzidfile.close();
    kzzidfile.open("Kzz.id");
    if (!kzzidfile) {
       exit_flag = 1;
       cerr << "File Kzz.id is missing!" << endl;
    }
    kzzidfile.close();

    if (exit_flag==1) { 
       CCTK_Abort(cctkGH, EXIT_FAILURE);
       exit(1);
    }

  } else {
    // These are for Uryu's initial data
    file >> *nx >> *ny >> *nz >> *nstar;
    file >> *dx1 >> *dy1 >> *dz1;
    file >> *xbd >> *ybd >> *zbd >> *rradi;
    file >> *rkappa >> *pinx >> *romega >> *radi;  
  }
  file.close();

// Convert the parameters xpc1 and xpc2 to fisheye coordinate
  if (xpc1==0.0) {
    cerr << "xpc1 hasn't been set properly in the par file" << endl;
    CCTK_Abort(cctkGH, EXIT_FAILURE);
    exit(1);
  }
  if (xpc2==0.0) {
    cerr << "xpc2 hasn't been set properly in the par file" << endl;
    CCTK_Abort(cctkGH, EXIT_FAILURE);
    exit(1);
  }
 
  double xpc1_fish_tmp,xpc2_fish_tmp; 
  if (fisheye_enable==0) {
    xpc1_fish_tmp = xpc1;
    xpc2_fish_tmp = xpc2;
  } else {
    ifstream file;
    file.open("RadialCoordinate_Input");
    int nFE; 
    double * aFE; double * r0FE; double * sFE; double * kappaFE;
    char buf[100],c;
    file.get(buf,100,'='); file.get(c); file >> nFE;
    aFE = new double[nFE+1];
    r0FE = new double[nFE+1];
    sFE = new double[nFE+1];
    kappaFE = new double[nFE+1];

    r0FE[0]=0.0; sFE[0]=0.0; kappaFE[0]=0.0;

    file.get(buf,100,'='); file.get(c); file >> aFE[0];
    if(nFE>0) {
      for(int i=1; i<=nFE; i++) file >> aFE[i];
      file.get(buf,100,'='); file.get(c); file >> r0FE[1];
      for(int i=2; i<=nFE; i++) file >> r0FE[i];
      file.get(buf,100,'='); file.get(c); file >> sFE[1];
      for(int i=2; i<=nFE; i++) file >> sFE[i];
    }
    file.close();
    for(int i=1; i<=nFE; i++) {     
       kappaFE[i] = (aFE[i-1]-aFE[i])*sFE[i]/(2.0*tanh(r0FE[i]/sFE[i]));
    }

    double x_max,x_min;
    int ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                            CCTK_VARIABLE_REAL,&x_max,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                            CCTK_VARIABLE_REAL,&x_min,1,CCTK_VarIndex("grid::X"));

    double Rpc1 = fabs(xpc1); 
    double Rpc2 = fabs(xpc2);

    // Looks like I have to do this stupid thing to set up xpc1_fish and xpc2_fish.
    double *rpc1, *rpc2;
    rpc1 = new double[1]; 
    rpc2 = new double[1];
    f_R_to_rfish(Rpc1, rpc1, x_max, aFE, kappaFE, r0FE, sFE, nFE);
    f_R_to_rfish(Rpc2, rpc2, fabs(x_min), aFE, kappaFE, r0FE, sFE, nFE);
    xpc1_fish_tmp = xpc1/Rpc1 * rpc1[0];
    xpc2_fish_tmp = xpc2/Rpc2 * rpc2[0];
    delete [] rpc1; delete [] rpc2;
    delete [] aFE; delete[] r0FE; delete[] sFE; delete [] kappaFE;
  } 

// Store xpc1_fish_tmp and xpc2_fish_tmp to GF vxl and vyl respectively so 
//  that they can be retrieved later. Yeah, I know this is ugly. 
  int ext[3];
  ext[0] = cctk_lsh[0];
  ext[1] = cctk_lsh[1];
  ext[2] = cctk_lsh[2];
  f_Set_xpc_fish(ext,vxl,vyl,xpc1_fish_tmp,xpc2_fish_tmp);

// *** TEST ***
//  cout << "xpc1, xpc2 = " << xpc1 << "  " << xpc2 << endl;
//  cout << "xpc1_fish, xpc2_fish = " << *xpc1_fish << "   " << *xpc2_fish << endl;
// ************
}

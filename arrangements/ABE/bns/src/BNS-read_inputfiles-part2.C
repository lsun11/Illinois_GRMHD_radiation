#include "cctk.h"
#include <stdio.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
//#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <fstream.h>

extern void ReadtoArray(const char * filename, double ** data,
			double & xmin, double & xmax, double & ymin, 
			double & ymax, double & zmin, double & zmax, 
			int & nx, int & ny, int & nz);

// translate matter for the irrotational bns problem

void Read_irr_Input(cGH *cctkGH,char * infile, int I_Print,int *nx,int *ny,int *nz,int *nstar,
		    double *dx,double *dY,double *dZ, 
		    double *xbd,double *ybd,double *zbd,double *rradi,
		    double *rkappa,double *pinx ,double *romega,double *radi);
//void weight(double Xdest,double Ydest,double Zdest, 
//    double Xgmin,double Ygmin,double Zgmin, 
//    double dX,double dY,double dZ,int i,int j,int k, 
//    int imin,int imax,int jmin,int jmax,int Nx,int Ny, 
//    int &il,int &iu,int &jl,int &ju,double &Delx,double &Dely);

void setup_irr(double dx_adm, double dY_adm, double dZ_adm,
	       double dx_ref, double dY_ref, double dZ_ref, 
	       double *lapse_s, double *phi_s,
	       double *shiftx_s, double *shifty_s, double *shiftz_s,
	       double *q_s,
	       double *ux_s, double *uy_s, double *uz_s,
	       int nx_adm, int ny_adm, int nz_adm,
	       int Nxh, int Nyh, int Nzh);

void push(int Nx, int Ny, int Nz, 
	  double symx, double symy, double symz, 
	  double *array);

void adjust_z_grid(int Nx, int Ny, int Nz, int *Nz_t, int nz_adm,
		   double dZ_adm, double dZ_ref, double symz, 
		   double *array_s, double ** array_t);

double interp_3d(double *f, int il, int jl, int kl,
		 int Nx, int Ny, int Nz_t,
		 double Dx, double DY, double DZ, int skin);
//void interp4(int ntot,int Nx,int Ny,int nx_adm,int ny_adm,double *f,double &f_int,int il,
//int iu,int jl,int ju,int kl,double &Dx,double &DY,double SymX,double SymY);
/*
  void translate_irr(char *cctkGH,int *ex,int Nx,int Ny,int Nz,int nx_adm,int ny_adm,int ntot,
  double *Xsource,double *Ysource,double *Zsource,
  double *Xdest, double *Ydest,double *Zdest,
  double *shiftx_new,double * shifty_new,double *shiftz_new,
  double *v2_new,double *lapse_new,
  double *phi_new,
  double *ux_new,double *uy_new,double *uz_new);
*/
double *q_i,*ux_i,*uy_i,*uz_i,*lapse_i,*phi_i;
double *shiftx_i,*shifty_i,*shiftz_i;

void Setup_Irrot_Data(cGH *cctkGH,double *PhysicalRadius,
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
		      double x_min, double y_min, double z_min) 
{
  int index,i,j,k;
  
  int Nxh=0,Nyh=0,Nzh=0;
  double xmin=0,xmax=0,ymin=0,ymax=0,zmin=0,zmax=0;

if (kt_id==1) {
  // Use Keisuke's initial data. Simply read the data from files.
/*
    double x_max,y_max,z_max,x_min,y_min,z_min;
    int ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                          CCTK_VARIABLE_REAL,x_max,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                          CCTK_VARIABLE_REAL,y_max,1,CCTK_VarIndex("grid::Y"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("maximum"), 1,
                          CCTK_VARIABLE_REAL,z_max,1,CCTK_VarIndex("grid::Z"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                          CCTK_VARIABLE_REAL,x_min,1,CCTK_VarIndex("grid::X"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                          CCTK_VARIABLE_REAL,y_min,1,CCTK_VarIndex("grid::Y"));
    ierr = CCTK_Reduce(cctkGH, -1,CCTK_ReductionHandle("minimum"), 1,
                          CCTK_VARIABLE_REAL,z_min,1,CCTK_VarIndex("grid::Z"));

   double dX = cctk_delta_space[0];
   double dY = cctk_delta_space[1];
   double dZ = cctk_delta_space[2];
*/

    f_read_KTid(ext, nx,ny,nz, dx_adm,dY_adm,dZ_adm, x_min,y_min,z_min,
                xdest,ydest,zdest,lapse,shiftx,shifty,shiftz,phi,
                Axx,Axy,Axz,Ayy,Ayz,Azz, q,ux,uy,uz, Symmetry);
} else {
  // Use Uryu's initial data
  
  //=====================================================
  // Read in the initial values of the fields from files
  //=====================================================
  char lapse_fname[64];  char shiftx_fname[64];
  char shifty_fname[64]; char shiftz_fname[64];
  char psi_fname[64];    char rho_fname[64];
  char ux_fname[64];     char uy_fname[64];
  char uz_fname[64];

  sprintf(lapse_fname,"alpha.dump");
  sprintf(shiftx_fname,"wx.dump");
  sprintf(shifty_fname,"wy.dump");
  sprintf(shiftz_fname,"wz.dump");
  sprintf(psi_fname,"Psi.dump");
  sprintf(rho_fname,"rho.dump");
  sprintf(ux_fname,"ux.dump");
  sprintf(uy_fname,"uy.dump");
  sprintf(uz_fname,"uz.dump");

  cout << "Reading in initial data from files" << endl; 

  double *lapse_s, *phi_s, *shiftx_s, *shifty_s, *shiftz_s;
  double *q_s, *ux_s, *uy_s, *uz_s;

  ReadtoArray(lapse_fname, &lapse_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(psi_fname, &phi_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(shiftx_fname, &shiftx_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(shifty_fname, &shifty_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);    
  ReadtoArray(shiftz_fname, &shiftz_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(rho_fname, &q_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(ux_fname, &ux_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(uy_fname, &uy_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);
  ReadtoArray(uz_fname, &uz_s,
	      xmin, xmax, ymin, ymax, zmin, zmax, 
	      Nxh, Nyh, Nzh);

  for(i=0; i<Nxh*Nyh*Nzh; i++)
    {
      phi_s[i] = log(phi_s[i]+1.0);
      shiftx_s[i] = -shiftx_s[i];
      shifty_s[i] = -shifty_s[i];
      shiftz_s[i] = -shiftz_s[i];
    }

  if(new_interpolation_choice == 1) {
    double dX = (xmax - xmin)/((double)Nxh-1.0);
    double dY = (ymax - ymin)/((double)Nyh-1.0);
    double dZ = (zmax - zmin)/((double)Nzh-1.0);
    double symx = 1.0, symy = 1.0, symz = 1.0;
    int m1 = 0, m2 = 0, m3 = 0, nn=1;
    int lbz;
    //FOLLOWING LINES WORK ONLY IN EQUATORIAL SYMMETRY!!!!!!!
    if(zdest[CCTK_GFINDEX3D(cctkGH,0,0,0)] < 0.0 && Symmetry==1) {
      lbz = -1;
      printf("I SHOULD BE HERE!\n");
    }
    else lbz = 0;

    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ, 
    			lapse_s, lapse,xdest,ydest,zdest,
			PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        phi_s, phi,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        q_s, q,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    symx = 1.0; symy = -1.0; symz = 1.0;
    m1 = 0; m2 = 1; m3 = 0; nn = 3;
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        shiftx_s, shiftx,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    symx = -1.0; symy = 1.0; symz = 1.0;
    m1 = 1; m2 = 0; m3 = 0; nn = 3;
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        shifty_s, shifty,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    symx = -1.0; symy = -1.0; symz = -1.0;
    m1 = 1; m2 = 1; m3 = 1; nn = 7;
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        shiftz_s, shiftz,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    symx = 1.0; symy = -1.0; symz = 1.0;
    m1 = 0; m2 = 0; m3 = 0; nn = 0;
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        ux_s, ux,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    symx = -1.0; symy = 1.0; symz = 1.0;
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        uy_s, uy,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);
    symx = -1.0; symy = -1.0; symz = -1.0;
    f_new_Translate_irr(ext,Nxh,Nyh,Nzh,xmax,ymax,zmax,dX,dY,dZ,
                        uz_s, uz,xdest,ydest,zdest,
                        PhysicalRadius,symx,symy,symz,m1,m2,m3,nn,lbz);    
  }
  else if(new_interpolation_choice == 0) {
    printf("WARNING: OLD INTERPOLATION SCHEME DOES NOT WORK ON MULTIPLE PROCS IN CACTUS!\n");


    push(Nxh,Nyh,Nzh,1.0,1.0,1.0,lapse_s);
    push(Nxh,Nyh,Nzh,1.0,1.0,1.0,phi_s);
    push(Nxh,Nyh,Nzh,1.0,-1.0,1.0,shiftx_s);
    push(Nxh,Nyh,Nzh,-1.0,1.0,1.0,shifty_s);
    push(Nxh,Nyh,Nzh,-1.0,-1.0,-1.0,shiftz_s);
    push(Nxh,Nyh,Nzh,1.0,1.0,1.0,q_s);
    push(Nxh,Nyh,Nzh,1.0,-1.0,1.0,ux_s);
    push(Nxh,Nyh,Nzh,-1.0,1.0,1.0,uy_s);
    push(Nxh,Nyh,Nzh,-1.0,-1.0,-1.0,uz_s);

    //  Read_irr_Input(infile, I_Print);
   
    //===============================================================
    // Note:  Nxh,Nyh,Nzh are the dimensions of the reference arrays.
    //===============================================================
    //UNUSED:
    //n_tot = Nxh*Nyh*Nzh*4;

    //=============================================================
    // Fill in the coordinate arrays for the reference matter grid.
    //=============================================================
    cout << "setting up reference coordinate grid" << endl;

    /*
      MPI_Comm_size(MPI_COMM_WORLD, &nproc);
      MPI_Comm_rank(MPI_COMM_WORLD, &thisproc);
    */
    /*
      int nproc=1;
      int thisproc=0;
      nproc = CCTK_nProcs(cctkGH);
      thisproc = CCTK_MyProc(cctkGH);
    */
    // THERE IS NO WAY THESE LINES WILL WORK IN GENERAL:
    /*
      int lbx = thisproc*(Nxh/nproc);
      int ubx = (thisproc+1)*(Nxh/nproc);
      double dX = (xmax - xmin)/((double)Nxh-1.0);
      int lby = thisproc*(Nyh/nproc);
      int uby = (thisproc+1)*(Nyh/nproc);
      double dY = (ymax - ymin)/((double)Nyh-1.0);
      int lbz = thisproc*(Nzh/nproc);
      int ubz = (thisproc+1)*(Nzh/nproc);
      double dZ = (zmax - zmin)/((double)Nzh-1.0);
    */
    //set dX,dY,dZ for data to-be-interpolated onto actual grid.
    double dX = (xmax - xmin)/((double)Nxh-1.0);
    double dY = (ymax - ymin)/((double)Nyh-1.0);
    double dZ = (zmax - zmin)/((double)Nzh-1.0);

    double * Xh; Xh = new double[Nxh*2];
    double * Yh; Yh = new double[Nyh*2];
    double * Zh; Zh = new double[Nzh*2];

    cout << "dY = " << dY << "    dZ = " << dZ << endl;

    for(i=0; i<Nxh; i++) {
      Xh[i+Nxh] = xmin + ((double)i+0.5) * dX;
      Xh[Nxh-i-1] = xmin - ((double)i+0.5) * dX;
    }
    for(i=0; i<Nyh; i++) {
      Yh[i+Nyh] = ymin + ((double)i+0.5) * dY;
      Yh[Nyh-i-1] = ymin - ((double)i+0.5) * dY;
    }
    for(i=0; i<Nzh; i++) {
      Zh[i+Nzh] = zmin + ((double)i+0.5) * dZ;
      Zh[Nzh-i-1] = zmin - ((double)i+0.5) * dZ;
    }
 
    lapse_i = new double[Nxh*Nyh*Nzh];  
    shiftx_i = new double[Nxh*Nyh*Nzh];  
    shifty_i = new double[Nxh*Nyh*Nzh];  
    shiftz_i = new double[Nxh*Nyh*Nzh];  
    phi_i = new double[Nxh*Nyh*Nzh];  
    q_i = new double[Nxh*Nyh*Nzh];  
    ux_i = new double[Nxh*Nyh*Nzh];  
    uy_i = new double[Nxh*Nyh*Nzh];  
    uz_i = new double[Nxh*Nyh*Nzh];  

    //============================================================================
    // Allocate memory for the destination arrays -- size determined by adm grid!
    //============================================================================
  
    int nx_adm = nx; int ny_adm = ny;  int nz_adm = nz;

    int n_tot_adm = nx_adm * ny_adm * nz_adm;
    /*
      int lx = (ub(0)-lb(0))+1;
      int ly = (ub(1)-lb(1))+1;
      int lz = (ub(2)-lb(2))+1;
    */
    int lx = ext[0];
    int ly = ext[1];
    int lz = ext[2];
    //  double * Xdest; Xdest = new double[(lx)*(ly)*(lz)]; 
    //double * Ydest; Ydest = new double[(lx)*(ly)*(lz)]; 
    //double * Zdest; Zdest = new double[(lx)*(ly)*(lz)]; 
    
    cout << "setting up initial conditions" << endl;
    //================================================================
    // Make _i arrays -- sized to fit the adm grid
    //================================================================

    cout << "Read_Irr_Input: Nxh, Nzh " << Nxh << " " << Nzh << endl;
    cout << "Read_Irr_Input: nx_adm,nz_adm  " << nx_adm << " " << nz_adm << endl;

    /*
      double dx_adm = cctk_delta_space[0];
      double dY_adm = cctk_delta_space[1];
      double dZ_adm = cctk_delta_space[2];
    */

    /*
      void setup_irr(char **cctkGH,double dx_adm, double dY_adm, double dZ_adm,
      double dX_ref, double dY_ref, double dZ_ref, 
      double *lapse_s, double *phi_s,
      double *shiftx_s, double *shifty_s, double *shiftz_s,
      double *q_s,
      double *ux_s, double *uy_s, double *uz_s,
      int nx_adm, int ny_adm, int nz_adm,
      int Nxh, int Nyh, int Nzh,
      double *lapse_i, double *phi_i, double *shiftx_i, double *shifty_i, double *shiftz_i,
      double *q_i, double *ux_i, double *uy_i, double *uz_i);
    */
    for(int i=0;i<Nxh*Nyh*Nzh;i++) {
      shiftx_i[i]=-1900;
    }
    setup_irr(dx_adm, dY_adm, dZ_adm,
	      dX, dY, dZ, 
	      lapse_s, phi_s,
	      shiftx_s, shifty_s, shiftz_s,
	      q_s,
	      ux_s, uy_s, uz_s,
	      nx_adm, ny_adm, nz_adm,
	      Nxh, Nyh, Nzh);
    /*
      setup_irr(cctkGH,dx_adm,dY_adm,dZ_adm,dX,dY,dZ,
      lapse_s,phi_s,shiftx_s,shifty_s,shiftz_s,
      q_s,ux_s,uy_s,uz_s,nx_adm,ny_adm,nz_adm,Nxh,Nyh,Nzh,
      lapse_i, phi_i, shiftx_i, shifty_i, shiftz_i,
      q_i, ux_i, uy_i, uz_i);
    */
    int Nzh_tot = Nzh;
    int Nxh_tot = Nxh;
    int Nyh_tot = Nyh;
    printf("hlkasjdiiiii %e\n",q_i[64]);
    //  Test();
    int lbz;
    //FOLLOWING LINES WORK ONLY IN EQUATORIAL SYMMETRY!!!!!!!
    if(zdest[CCTK_GFINDEX3D(cctkGH,0,0,0)] < 0.0 && Symmetry==1) {
      lbz = -1;
      printf("I SHOULD BE HERE!\n");
    }
    else lbz = 0;
    f_Translate_irr(ext,Nxh_tot,Nyh_tot,Nzh_tot,nx_adm,ny_adm,n_tot_adm,
		    Xh,Yh,Zh,xdest,ydest,zdest,
		    shiftx_i,shifty_i,shiftz_i,
		    shiftx,shifty,shiftz,
		    q_i,q,lapse_i,lapse, 
		    phi_i,phi,
		    ux_i,uy_i,uz_i,ux,uy,uz,lbz);
    //		lbx,lby,lbz); 

    /*
      Translate_irr(cctkGH,ext,Nxh_tot,Nyh_tot,Nzh_tot,nx_adm,ny_adm,n_tot_adm,
      Xh,Yh,Zh,xdest,ydest,zdest,
      shiftx,shifty,shiftz,
      q,lapse, 
      phi,
      ux,uy,uz);
    */
    printf("hlkasjdiiiii2 %e\n",q_i[64]);
    printf("hlkasjdiiiii5 %e\n",q[CCTK_GFINDEX3D(cctkGH,1,1,1)]);
    //exit(0);
    //		  lbx,lby,lbz); 

    /*
      int ir, jr, kr;
      double symx, symy;
      for(i=0; i<lx; i++)
      for(j=0; j<ly; j++)
      for(k=0; k<lz; k++)
      {
      //int index2 = lx*ly*(k-lb(2)) + lx*(j-lb(1)) + i-lb(0);
      index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      Xdest[index2] = x[index];  Ydest[index2] = y(index);  Zdest[index2] = z(index);
      }
    
      int lbx = lb(0), lby = lb(1), lbz = lb(2);
      int Nzh_tot = Nzh;
      int Nxh_tot = Nxh;
      int Nyh_tot = Nyh;
    
      //=========================================================
      // Copy reference arrays to initial conditions  
      //=========================================================
      f_translate_irr(ext,Nxh_tot,Nyh_tot,Nzh_tot,nx_adm,ny_adm,n_tot_adm,
      Xh,Yh,Zh,Xdest,Ydest,Zdest,
      shiftx_i,shifty_i,shiftz_i,
      shiftx,shifty,shiftz,
      q_i,q,lapse_i,lapse, 
      phi_i,phi,
      ux_i,uy_i,uz_i,ux,uy,uz,
      lbx,lby,lbz); 

      double q_atmos = q->GF_maxval(t0,level,IDENT);
      rho_b_atm = q_atmos * rho_fact;
      //    q_atmos = q_atmos * rho_fact;
      //    ambient = q_atmos;
      //    atm_radius = 0.e0; 
    
      //  TEMPORARY: Polytropic Index hardwired here!
      double n = 1.0;
    
      f_compute_ibwh(ext,Omega_Frame,q_atmos,n, 
      X,Y,Z,
      q,lapse,phi,
      shiftx,shifty,shiftz,
      ux,uy,uz,ut,
      rho,S,Sx,Sy,Sz,
      Sxx,Sxy,Sxz,Syy,Syz,Szz,
      rho_star,tau,st_x);
    */
  } //end if(new_interpolation_choice
 } // end if (kt_id==1)
}
/*
  void Read_bin_Input(char * infile, int I_Print, 
  double & re, double & C, double & Omega,
  double & rho_fact, double & C_vis)
  {
  if (I_Print) {
  ifstream file;
  file.open(infile);
  if (!file) {
  cerr << "Can't open " << infile << " for input." << endl;
  CCTK_Abort(cctkGH, EXIT_FAILURE);
  exit(1);
  }
  char buf[600],c;
  file.get(buf,600,'='); file.get(c); file >> re;
  file.get(buf,600,'='); file.get(c); file >> C;
  file.get(buf,600,'='); file.get(c); file >> Omega;
  file.get(buf,600,'='); file.get(c); file >> rho_fact;  
  file.get(buf,600,'='); file.get(c); file >> C_vis;  
  if(file.eof()) {
  cerr << "Error reading input file " << infile << endl;
  CCTK_Abort(cctkGH, EXIT_FAILURE);
  exit(1);
  }
  //
  // echo out input data:
  //
  cout << "  re:  " << re << endl;
  cout << "  C:  " << C << endl;
  cout << "  Omega:  " << Omega << endl;
  cout << "  rho_fact:  " << rho_fact << endl;
  cout << "  C_vis:  " << C_vis << endl;
  }
  MPI_Bcast(&re,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&C,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&Omega,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&rho_fact,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&C_vis,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
*/


//=============================================================================
//  This subroutine converts the arrays read in from files into the reference 
//  arrays to be used by the evolution code.  It performs two tasks.  First, 
//  the arrays are interpolated in the z-direction so that the z = const
//  planes in the interpolated reference arrays correspond to the z = const 
//  planes in the adm grid.  Second, the code extrapolates the reference 
//  arrays in all three directions so that it covers the whole range 
//  [0,xmax_adm*sqrt(2)] [0,ymax_adm*sqrt(2)] [0,zmax_adm]
//  The extrapolation is accomplished using the asymptotic forms of the 
//  metric functions.
//=============================================================================
void setup_irr(double dx_adm, double dY_adm, double dZ_adm,
	       double dX_ref, double dY_ref, double dZ_ref, 
	       double *lapse_s, double *phi_s,
	       double *shiftx_s, double *shifty_s, double *shiftz_s,
	       double *q_s,
	       double *ux_s, double *uy_s, double *uz_s,
	       int nx_adm, int ny_adm, int nz_adm,
	       int Nxh, int Nyh, int Nzh)
{

  int i,j,k;
  double xmax_adm = dx_adm * (nx_adm/2 - 1) + 0.5*dx_adm;
  double ymax_adm = dY_adm * (ny_adm/2 - 1) + 0.5*dY_adm;
  double zmax_adm = dZ_adm * (nz_adm - 1) + 0.5*dZ_adm;
  int Nz_t;
  double *lapse_t, *phi_t, *shiftx_t, *shifty_t, *shiftz_t;
  double *q_t, *ux_t, *uy_t, *uz_t; 
  int I_Print = 1;//( comm_service::proc_me() == 0 );
  //==============================================
  // Option 1a:  The xy planes alreadY line up.
  //==============================================
  if (fabs(dZ_adm-dZ_ref)<1e-4*dZ_adm) {
    if(I_Print) cout << "No z-interpolation necessary." << endl;
    Nz_t = Nzh;
    lapse_t  = new double[Nxh*Nyh*Nz_t];
    phi_t    = new double[Nxh*Nyh*Nz_t];
    shiftx_t = new double[Nxh*Nyh*Nz_t];
    shifty_t = new double[Nxh*Nyh*Nz_t];
    shiftz_t = new double[Nxh*Nyh*Nz_t];
    q_t      = new double[Nxh*Nyh*Nz_t];
    ux_t     = new double[Nxh*Nyh*Nz_t];
    uy_t     = new double[Nxh*Nyh*Nz_t];
    uz_t     = new double[Nxh*Nyh*Nz_t];
    for(i=0; i<Nxh*Nyh*Nz_t; i++) {
      lapse_t[i]  = lapse_s[i];
      phi_t[i]    = phi_s[i];
      shiftx_t[i] = shiftx_s[i];
      shifty_t[i] = shifty_s[i];
      shiftz_t[i] = shiftz_s[i];
      q_t[i]      = q_s[i];
      ux_t[i]     = ux_s[i];
      uy_t[i]     = uy_s[i];
      uz_t[i]     =  uz_s[i];
      if (q_t[i]!=0)
	cout << "i,q_t[i]:"<<i<<" "<<q_t[i]<<endl;

    }
  }
  //==============================================
  // Option 2a:  The xy planes do not line up.
  //==============================================
  else {
    if(I_Print) cout << "hrl\t" << fabs(dZ_adm-dZ_ref)<< '\t' <<1e-4*dZ_adm << endl;

    //    if(I_Print) 
    cout << "dZ_adm not equal to dZ_ref:  interpolating reference grid" << endl;
    cout << "dZ_adm = " << dZ_adm << ",  dZ_ref = " << dZ_ref << endl;

    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,lapse_s,&lapse_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,phi_s,&phi_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,shiftx_s,&shiftx_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,shifty_s,&shifty_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,-1.0,shiftz_s,&shiftz_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,q_s,&q_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,ux_s,&ux_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,1.0,uy_s,&uy_t);
    adjust_z_grid(Nxh,Nyh,Nzh,&Nz_t,nz_adm,dZ_adm,dZ_ref,-1.0,uz_s,&uz_t);
    dZ_ref = dZ_adm;
  }
  double xmax_ref = dX_ref * (Nxh-1) + 0.5*dX_ref;
  double ymax_ref = dY_ref * (Nyh-1) + 0.5*dY_ref;
  double zmax_ref = dZ_ref * (Nz_t-1) + 0.5*dZ_adm;
  double x_tmax = sqrt(2)*xmax_adm + 2.0*dx_adm;
  double y_tmax = sqrt(2)*ymax_adm + 2.0*dY_adm;
  double z_tmax = zmax_adm + 2.0*dZ_adm;
  
  int i_rt = floor((xmax_ref-0.5*dx_adm)/dx_adm);
  int j_rt = floor((ymax_ref-0.5*dY_adm)/dY_adm);
  int k_rt = floor((zmax_ref-0.5*dZ_adm)/dZ_adm);
  int Nxi,Nyi,Nzi;
  if(I_Print) {
    cout << "xmax_ref = " << xmax_ref 
	 << ", x_tmax = " << x_tmax 
	 << ", zmax_ref = " << zmax_ref
	 << ", z_tmax = " << z_tmax << endl;
    cout << "xmax_adm = " << xmax_adm 
	 << ", dx_adm = " << dx_adm
	 << ", nx_adm = " << nx_adm << endl;
    cout << "zmax_adm = " << zmax_adm 
	 << ", dZ_adm = " << dZ_adm
	 << ", nz_adm = " << nz_adm << endl;
  }
  //==============================================
  // Option 1b:  The grid is alreadY big enough.
  //============================================== 
  if((xmax_ref > x_tmax)&&(ymax_ref > y_tmax)&&(zmax_ref > z_tmax))
    {
      if(I_Print) cout << "Grid does not need to be extended." << endl;

      lapse_i = new double[Nxh*Nyh*Nz_t];
      phi_i = new double[Nxh*Nyh*Nz_t];
      shiftx_i = new double[Nxh*Nyh*Nz_t];
      shifty_i = new double[Nxh*Nyh*Nz_t];
      shiftz_i = new double[Nxh*Nyh*Nz_t];
      q_i = new double[Nxh*Nyh*Nz_t];
      ux_i = new double[Nxh*Nyh*Nz_t];
      uy_i = new double[Nxh*Nyh*Nz_t];
      uz_i = new double[Nxh*Nyh*Nz_t];
      Nxi = Nxh;  Nyi = Nyh;  Nzi = Nz_t;

      for(i=0; i<Nxh*Nyh*Nz_t; i++) {
	lapse_i[i]  = lapse_t[i];
	phi_i[i]    = phi_t[i];
	shiftx_i[i] = shiftx_t[i];
	shifty_i[i] = shifty_t[i];
	shiftz_i[i] = shiftz_t[i];
	q_i[i]      = q_t[i];
	ux_i[i]     = ux_t[i];
	uy_i[i]     = uy_t[i];
	uz_i[i]     =  uz_t[i];
      }
    }
  //==============================================
  // Option 2b:  The grid is not yet big enough.
  //==============================================
  else
    {
      if(I_Print) cout << "Grid is being extended." << endl;
      Nxi = floor(x_tmax/dX_ref + 0.5);
      Nyi = floor(y_tmax/dY_ref + 0.5);
      Nzi = floor(z_tmax/dZ_ref + 0.5);
      if(I_Print) cout << "Nxi = " << Nxi << ", Nyi = " << Nyi 
		       << ", Nzi = " << Nzi 
		       << ", Nxh = " << Nxh << endl;

      lapse_i  = new double[Nxi*Nyi*Nzi];
      phi_i = new double[Nxi*Nyi*Nzi];
      shiftx_i = new double[Nxi*Nyi*Nzi];
      shifty_i = new double[Nxi*Nyi*Nzi];
      shiftz_i = new double[Nxi*Nyi*Nzi];
      q_i = new double[Nxi*Nyi*Nzi];
      ux_i = new double[Nxi*Nyi*Nzi];
      uy_i = new double[Nxi*Nyi*Nzi];
      uz_i = new double[Nxi*Nyi*Nzi];

      //===========================================
      // Fill in the grids
      //===========================================
      double r1,x1,y1,z1,rx,ry,rz,r2;
      double x2,y2,z2,DX,DY,DZ,xl,yl,zl,xu,yu,zu;
      double lapse_int,phi_int,shiftx_int,shifty_int,shiftz_int;
      double shiftx_int2,shifty_int2,shiftz_int2,r22,x22,y22,z22;
      double r1_3,r1_5,r1_7,r2_3,r2_5,r2_7,r22_3,r22_5,r22_7;
      double m11,m12,m21,m22,det,gx,gy,gz,bx,by,bz;
      int index,index2,il,jl,kl;
      for(i=0; i<Nxi; i++)
	for(j=0; j<Nyi; j++)
	  for(k=0; k<Nzi; k++) {
	    index = Nxi*Nyi*k + Nxi*j + i;
	    index2 = Nxh*Nyh*k + Nxh*j + i;
	    //
	    // Fill in the interior grid.
	    if((i<Nxh)&&(j<Nyh)&&(k<Nz_t)) {
	      lapse_i[index]  = lapse_t[index2];
	      phi_i[index]    = phi_t[index2];
	      shiftx_i[index] = shiftx_t[index2];
	      shifty_i[index] = shifty_t[index2];
	      shiftz_i[index] = shiftz_t[index2];
	      q_i[index]      = q_t[index2];
	      ux_i[index]     = ux_t[index2];
	      uy_i[index]     = uy_t[index2];
	      uz_i[index]     = uz_t[index2];
	    }
	    else {
	      //
	      // Fill in the exterior grid.
	      x1 = (i+0.5)*dX_ref;
	      y1 = (j+0.5)*dY_ref;
	      z1 = (k+0.5)*dZ_ref; // Note: dZ_ref = dZ_adm, by earlier construction.
	      r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
	      rx = r1*xmax_ref/x1;  ry = r1*ymax_ref/y1;  rz = r1*zmax_ref/z1;
	      //
	      // Option i:  Radial line to r 
	      //            intersects yz-face of original reference cube.
	      if((rx<=ry)&&(rx<=rz)) {
		r2 = rx;
		x2 = (r2/r1)*x1;  y2 = (r2/r1)*y1;  z2 = (r2/r1)*z1; 
		il = Nxh-1;
		jl = floor(y2/dY_ref - 0.5); kl = floor(z2/dZ_ref - 0.5);
		if(jl>Nyh-3) jl = Nyh-3;
		if(kl>Nz_t-3) kl = Nz_t-3;
		yl = (0.5 + jl)*dY_ref;  zl = (0.5 + kl)*dZ_ref;
		yu = yl + dY_ref;  zu = zl + dZ_ref;
		DZ = (z2 - zl)/(zu - zl);
		DY = (y2 - yl)/(yu - yl);
		DX = 0.0;
		lapse_int  = interp_3d(lapse_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		phi_int    = interp_3d(phi_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shiftx_int = interp_3d(shiftx_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shifty_int = interp_3d(shifty_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shiftz_int = interp_3d(shiftz_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		x22 = x2 - dx_adm;  y22 = y2;  z22 = z2;
		r22 = sqrt(x22*x22 + y22*y22 + z22*z22);
		shiftx_int2 = interp_3d(shiftx_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
		shifty_int2 = interp_3d(shifty_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
		shiftz_int2 = interp_3d(shiftz_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
	      }
	      //
	      // Option ii:  Radial line to r 
	      //            intersects xz-face of original reference cube.
	      if((ry<rx)&&(ry<=rz)) {
		r2 = ry;
		x2 = (r2/r1)*x1; y2 = (r2/r1)*y1;  z2 = (r2/r1)*z1; 
		jl = Nyh-1;
		il = floor(x2/dX_ref - 0.5); kl = floor(z2/dZ_ref - 0.5);
		if(kl>Nz_t-3) kl = Nz_t-3;
		if(il>Nxh-3) il = Nxh-3;
		xl = (0.5 + il)*dX_ref;  zl = (0.5 + kl)*dZ_ref;
		xu = xl + dX_ref;  zu = zl + dZ_ref;
		DZ = (z2 - zl)/(zu - zl);
		DX = (x2 - xl)/(xu - xl);
		DY = 0.0;
		lapse_int  = interp_3d(lapse_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		phi_int    = interp_3d(phi_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shiftx_int = interp_3d(shiftx_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shifty_int = interp_3d(shifty_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shiftz_int = interp_3d(shiftz_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);     
		x22 = x2;  y22 = y2 - dY_adm;  z22 = z2;
		r22 = sqrt(x22*x22 + y22*y22 + z22*z22);
		shiftx_int2 = interp_3d(shiftx_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
		shifty_int2 = interp_3d(shifty_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
		shiftz_int2 = interp_3d(shiftz_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
	      }
	      //
	      // Option iii:  Radial line to r 
	      //            intersects xy-face of original reference cube.
	      if((rz<rx)&&(rz<ry)) {
		r2 = rz;
		x2 = (r2/r1)*x1; y2 = (r2/r1)*y1;  z2 = (r2/r1)*z1; 
		kl = Nz_t-1;
		il = floor(x2/dX_ref - 0.5); jl = floor(y2/dY_ref - 0.5);
		if(jl>Nyh-3) jl = Nyh-3;
		if(il>Nxh-3) il = Nxh-3;
		xl = (0.5 + il)*dX_ref;  yl = (0.5 + jl)*dY_ref;
		xu = xl + dX_ref;  yu = yl + dY_ref;
		DY = (y2 - yl)/(yu - yl);
		DX = (x2 - xl)/(xu - xl);
		DZ = 0.0;
		lapse_int  = interp_3d(lapse_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		phi_int    = interp_3d(phi_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shiftx_int = interp_3d(shiftx_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shifty_int = interp_3d(shifty_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);
		shiftz_int = interp_3d(shiftz_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,0);     
		x22 = x2;  y22 = y2;  z22 = z2 - dZ_adm;
		r22 = sqrt(x22*x22 + y22*y22 + z22*z22);
		shiftx_int2 = interp_3d(shiftx_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
		shifty_int2 = interp_3d(shifty_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
		shiftz_int2 = interp_3d(shiftz_t,il,jl,kl,Nxh,Nyh,Nz_t,DX,DY,DZ,1);
	      }
	      //====================================================================
	      //  Assumed Asymptotic Behaviors:
	      //    lapse  ~ 1/r              phi   ~ 1/r
	      //    shiftx ~  y/r^3          shifty ~ x/r^3
	      //    shiftz ~ xyz/r^7
	      //====================================================================
	      q_i[index] = ux_i[index] = uy_i[index] = uz_i[index] = 0.0;
	      lapse_i[index] = (r2/r1)*lapse_int;
	      phi_i[index] = (r2/r1)*phi_int;
	      // Quick and dirty estimate of shift:
  	      // shiftx_i[index] = pow((r2/r1),3)*(y1/y2)*shiftx_int;
  	      // shifty_i[index] = pow((r2/r1),3)*(x1/x2)*shifty_int;
  	      // shiftz_i[index] = pow((r2/r1),7)*((x1*y1*z1)/(x2*y2*z2))*shiftz_int;
	      // More rigorous estimate:
	      r2_3 = pow(r2,3);  r2_5 = r2_3*r2*r2;  r2_7 = r2_5*r2*r2;
	      r22_3 = pow(r22,3);  r22_5 = r22_3*r22*r22;  r22_7 = r22_5*r22*r22;
	      m11 = y2/r2_3;  m12 = x2*x2*y2/r2_5;
	      m21 = y22/r22_3;  m22 = x22*x22*y22/r22_5;
	      det = m11*m22 - m12*m21;
	      gx = (shiftx_int*m22 - shiftx_int2*m12)/det;
	      bx = (shiftx_int2*m11 - shiftx_int*m21)/det;
	      m11 = x2/r2_3;  m12 = y2*y2*x2/r2_5;
	      m21 = x22/r22_3;  m22 = y22*y22*x22/r22_5;
	      det = m11*m22 - m12*m21;
	      gy = (shifty_int*m22 - shifty_int2*m12)/det;
	      by = (shifty_int2*m11 - shifty_int*m21)/det;
	      m11 = x2*y2*z2/r2_5;  m12 = x2*y2*z2/r2_7;
	      m21 = x22*y22*z22/r22_5;  m22 = x22*y22*z22/r22_7;
	      det = m11*m22 - m12*m21;
	      gz = (shiftz_int*m22 - shiftz_int2*m12)/det;
	      bz = (shiftz_int2*m11 - shiftz_int*m21)/det;
	      r1_3 = r1*r1*r1;  r1_5 = r1_3*r1*r1;  r1_7 = r1_5*r1*r1;
	      shiftx_i[index] = gx*y1/r1_3 + bx*x1*x1*y1/r1_5;
	      shifty_i[index] = gy*x1/r1_3 + by*x1*y1*y1/r1_5;
	      shiftz_i[index] = gz*x1*y1*z1/r1_7 + bz*x1*y1*z1/r1_5;
	    }
	  }
    }
  //======================================================
  // Make final reference coordinate arrays.
  //======================================================
  /*
    Nxh = Nxi;  Nyh = Nyi;  Nzh = Nzi;
    double * Xh; Xh = new double[Nxh*2];
    double * Yh; Yh = new double[Nyh*2];
    double * Zh; Zh = new double[Nzh*2];

    for(i=0; i<Nxh; i++) {
    Xh[i+Nxh] = ((double)i+0.5) * dX_ref;
    Xh[Nxh-i-1] = - ((double)i+0.5) * dX_ref;
    }
    for(i=0; i<Nyh; i++) {
    Yh[i+Nyh] = ((double)i+0.5) * dY_ref;
    Yh[Nyh-i-1] = - ((double)i+0.5) * dY_ref;
    }
    for(i=0; i<Nzh; i++) {
    Zh[i+Nzh] = ((double)i+0.5) * dZ_ref;
    Zh[Nzh-i-1] = - ((double)i+0.5) * dZ_ref;
    }
  */
  cout << "reference grid:  dY = " << dY_ref << "    dZ = " << dZ_ref << endl;
}
//============================================================
//  "Push" data onto a cell-centered grid.
//   Note:  for vertex-centered data:
//                                    f[-1] = f[1], etc.
//============================================================
void push(int Nx, int Ny, int Nz, 
	  double symx, double symy, double symz, 
	  double *array)
{
  double *temp;
  temp = new double[Nx*Ny*Nz];
  int tx, ty, tz, i;
  //
  // z push
  //
  for (i = 0; i<Nx*Ny*Nz; i++)
    {
      tz = floor(i/(Nx*Ny));
      ty = floor((i-tz*Nx*Ny)/Nx);
      tx = i - Nx*Ny*tz - Nx*ty;
      if(tz==Nz-1)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  1.5*array[Nx*Ny*tz + Nx*ty + tx] - 0.5*array[Nx*Ny*(tz-1) + Nx*ty + tx];
      else if(tz==Nz-2)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  0.5*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*(tz+1) + Nx*ty + tx]); 	
      else if(tz==0)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (9.0/16.0)*(array[Nx*ty + tx] + array[Nx*Ny + Nx*ty + tx])
	  - (1.0/16.0)*(symz*array[Nx*Ny + Nx*ty + tx]
			+ array[Nx*Ny*2 + Nx*ty + tx]);	
      else
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (9.0/16.0)*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*(tz+1) + Nx*ty + tx])
	  - (1.0/16.0)*(array[Nx*Ny*(tz-1) + Nx*ty + tx]
			+ array[Nx*Ny*(tz+2) + Nx*ty + tx]);
    }
  for (i = 0; i<Nx*Ny*Nz; i++)  
    array[i] = temp[i];
  //
  // y push
  //
  for (i = 0; i<Nx*Ny*Nz; i++)
    {
      tz = floor(i/(Nx*Ny));
      ty = floor((i-tz*Nx*Ny)/Nx);
      tx = i - Nx*Ny*tz - Nx*ty;
      if(ty==Ny-1)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  1.5*array[Nx*Ny*tz + Nx*ty + tx] - 0.5*array[Nx*Ny*tz + Nx*(ty-1) + tx];
      else if(ty==Ny-2)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  0.5*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*tz + Nx*(ty+1) + tx]); 	
      else if(ty==0)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (9.0/16.0)*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*tz + Nx*(ty+1) + tx])
	  - (1.0/16.0)*(symy*array[Nx*Ny*tz + Nx*(ty+1) + tx]
			+ array[Nx*Ny*tz + Nx*(ty+2) + tx]);	
      else
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (9.0/16.0)*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*tz + Nx*(ty+1) + tx])
	  - (1.0/16.0)*(symy*array[Nx*Ny*tz + Nx*(ty-1) + tx]
			+ array[Nx*Ny*tz + Nx*(ty+2) + tx]);
    }
  for (i = 0; i<Nx*Ny*Nz; i++)  
    array[i] = temp[i];
  //
  // x push
  //
  for (i = 0; i<Nx*Ny*Nz; i++)
    {
      tz = floor(i/(Nx*Ny));
      ty = floor((i-tz*Nx*Ny)/Nx);
      tx = i - Nx*Ny*tz - Nx*ty;
      if(tx==Nx-1)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  1.5*array[Nx*Ny*tz + Nx*ty + tx] - 0.5*array[Nx*Ny*tz + Nx*ty + tx-1];
      else if(tx==Nx-2)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  0.5*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*tz + Nx*ty + tx+1]); 	
      else if(tx==0)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (9.0/16.0)*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*tz + Nx*ty + tx+1])
	  - (1.0/16.0)*(symx*array[Nx*Ny*tz + Nx*ty + tx+1]
			+ array[Nx*Ny*tz + Nx*ty + tx+2]);	
      else
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (9.0/16.0)*(array[Nx*Ny*tz + Nx*ty + tx] + array[Nx*Ny*tz + Nx*ty + tx+1])
	  - (1.0/16.0)*(array[Nx*Ny*tz + Nx*ty + tx-1]
			+ array[Nx*Ny*tz + Nx*ty + tx+2]);
    }
  for (i = 0; i<Nx*Ny*Nz; i++)  
    array[i] = temp[i];
}  

//=============================================================
// Do 2D cubic interpolation, 
// not knowing a priori onto which plane one is interpolating.
//=============================================================
double interp_3d(double *f, int il, int jl, int kl,
		 int Nx, int Ny, int Nz_t,
		 double DX, double DY, double DZ, int skin)
{
  double array_int;
  int iu,iuu,ill,ju,juu,jll,ku,kuu,kll;
  int ll_ll,ll_l,ll_u,ll_uu,l_ll,l_l;
  int l_u,l_uu,u_ll,u_l,u_u,u_uu,uu_ll,uu_l,uu_u,uu_uu;
  double f_1,f_2,f_3,f_4,f_int;
  //
  // yz-interpolation
  if (fabs(DX)<1e-5) {
    ku = kl + 1; kuu = ku + 1;  kll = kl - 1;  
    ju = jl + 1; juu = ju + 1;  jll = jl - 1;
    if (kl<0) kl = -kl - 1;
    if (jl<0) jl = -jl - 1;
    if (kll<0) kll = -kll - 1;
    if (jll<0) jll = -jll - 1;
    ll_ll  = Ny*Nx*(kll) + Nx*(jll) + (Nx-1-skin); 
    ll_l   = Ny*Nx*(kll) + Nx*jl    + (Nx-1-skin); 
    ll_u   = Ny*Nx*(kll) + Nx*ju    + (Nx-1-skin); 
    ll_uu  = Ny*Nx*(kll) + Nx*(juu) + (Nx-1-skin); 
    l_ll   = Ny*Nx*(kl) + Nx*(jll) + (Nx-1-skin);  
    l_l    = Ny*Nx*(kl) + Nx*jl    + (Nx-1-skin);  
    l_u    = Ny*Nx*(kl) + Nx*ju    + (Nx-1-skin);  
    l_uu   = Ny*Nx*(kl) + Nx*(juu) + (Nx-1-skin);  
    u_ll   = Ny*Nx*(ku) + Nx*(jll) + (Nx-1-skin);  
    u_l    = Ny*Nx*(ku) + Nx*jl    + (Nx-1-skin);  
    u_u    = Ny*Nx*(ku) + Nx*ju    + (Nx-1-skin);  
    u_uu   = Ny*Nx*(ku) + Nx*(juu) + (Nx-1-skin);   
    uu_ll  = Ny*Nx*(kuu) + Nx*(jll) + (Nx-1-skin); 
    uu_l   = Ny*Nx*(kuu) + Nx*jl    + (Nx-1-skin); 
    uu_u   = Ny*Nx*(kuu) + Nx*ju    + (Nx-1-skin); 
    uu_uu  = Ny*Nx*(kuu) + Nx*(juu) + (Nx-1-skin);
    f_1 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[ll_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[ll_l] +
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[ll_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[ll_uu];

    f_2 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[l_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[l_l] + 
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[l_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[l_uu];
    
    f_3 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[u_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[u_l] + 
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[u_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[u_uu];

    f_4 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[uu_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[uu_l] + 
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[uu_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[uu_uu];

    array_int = (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*f_1 + 
      ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*f_2 + 
      (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*f_3 + 
      (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*f_4;
  }
  //
  // xz-interpolation
  else if (fabs(DY)<1e-5) {
    iu = il + 1; iuu = iu + 1;  ill = il - 1;  
    ku = kl + 1; kuu = ku + 1;  kll = kl - 1;
    if (ill<0) ill = -ill - 1;
    if (kll<0) kll = -kll - 1;
    if (il<0) il = -il - 1;
    if (kl<0) kl = -kl - 1;
    ll_ll  = Ny*Nx*(kll) + Nx*(Ny-1-skin) + (ill); 
    ll_l   = Ny*Nx*(kl)  + Nx*(Ny-1-skin) + (ill); 
    ll_u   = Ny*Nx*(ku)  + Nx*(Ny-1-skin) + (ill); 
    ll_uu  = Ny*Nx*(kuu) + Nx*(Ny-1-skin) + (ill); 
    l_ll   = Ny*Nx*(kll) + Nx*(Ny-1-skin) + (il);  
    l_l    = Ny*Nx*(kl)  + Nx*(Ny-1-skin) + (il);  
    l_u    = Ny*Nx*(ku)  + Nx*(Ny-1-skin) + (il);  
    l_uu   = Ny*Nx*(kuu) + Nx*(Ny-1-skin) + (il);  
    u_ll   = Ny*Nx*(kll) + Nx*(Ny-1-skin) + (iu);  
    u_l    = Ny*Nx*(kl)  + Nx*(Ny-1-skin) + (iu);  
    u_u    = Ny*Nx*(ku)  + Nx*(Ny-1-skin) + (iu);  
    u_uu   = Ny*Nx*(kuu) + Nx*(Ny-1-skin) + (iu);   
    uu_ll  = Ny*Nx*(kll) + Nx*(Ny-1-skin) + (iuu); 
    uu_l   = Ny*Nx*(kl)  + Nx*(Ny-1-skin) + (iuu); 
    uu_u   = Ny*Nx*(ku)  + Nx*(Ny-1-skin) + (iuu); 
    uu_uu  = Ny*Nx*(kuu) + Nx*(Ny-1-skin) + (iuu);
    f_1 = (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*f[ll_ll] + 
      ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*f[ll_l] +
      (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*f[ll_u] + 
      (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*f[ll_uu];

    f_2 = (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*f[l_ll] + 
      ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*f[l_l] + 
      (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*f[l_u] + 
      (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*f[l_uu];
    
    f_3 = (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*f[u_ll] + 
      ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*f[u_l] + 
      (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*f[u_u] + 
      (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*f[u_uu];

    f_4 = (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*f[uu_ll] + 
      ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*f[uu_l] + 
      (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*f[uu_u] + 
      (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*f[uu_uu];

    array_int = (DX*(DX-1.0)*(2.0-DX)/6.0)*f_1 + 
      ((1.0+DX)*(1.0-DX)*(2.0-DX)/2.0)*f_2 + 
      (DX*(1.0+DX)*(2.0-DX)/2.0)*f_3 + 
      (DX*(1.0+DX)*(DX-1.0)/6.0)*f_4;
  }
  //
  // xy-interpolation
  else if (fabs(DZ)<1e-5) {
    iu = il + 1; iuu = iu + 1;  ill = il - 1;  
    ju = jl + 1; juu = ju + 1;  jll = jl - 1;
    if (il<0) il = -il - 1;
    if (jl<0) jl = -jl - 1;
    if (ill<0) ill = -ill - 1;
    if (jll<0) jll = -jll - 1;
    ll_ll  = Ny*Nx*(Nz_t-1-skin) + Nx*(jll) + (ill); 
    ll_l   = Ny*Nx*(Nz_t-1-skin) + Nx*jl    + (ill); 
    ll_u   = Ny*Nx*(Nz_t-1-skin) + Nx*ju    + (ill); 
    ll_uu  = Ny*Nx*(Nz_t-1-skin) + Nx*(juu) + (ill); 
    l_ll   = Ny*Nx*(Nz_t-1-skin) + Nx*(jll) + (il);  
    l_l    = Ny*Nx*(Nz_t-1-skin) + Nx*jl    + (il);  
    l_u    = Ny*Nx*(Nz_t-1-skin) + Nx*ju    + (il);  
    l_uu   = Ny*Nx*(Nz_t-1-skin) + Nx*(juu) + (il);  
    u_ll   = Ny*Nx*(Nz_t-1-skin) + Nx*(jll) + (iu);  
    u_l    = Ny*Nx*(Nz_t-1-skin) + Nx*jl    + (iu);  
    u_u    = Ny*Nx*(Nz_t-1-skin) + Nx*ju    + (iu);  
    u_uu   = Ny*Nx*(Nz_t-1-skin) + Nx*(juu) + (iu);   
    uu_ll  = Ny*Nx*(Nz_t-1-skin) + Nx*(jll) + (iuu); 
    uu_l   = Ny*Nx*(Nz_t-1-skin) + Nx*jl    + (iuu); 
    uu_u   = Ny*Nx*(Nz_t-1-skin) + Nx*ju    + (iuu); 
    uu_uu  = Ny*Nx*(Nz_t-1-skin) + Nx*(juu) + (iuu);
    f_1 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[ll_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[ll_l] +
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[ll_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[ll_uu];

    f_2 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[l_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[l_l] + 
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[l_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[l_uu];
    
    f_3 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[u_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[u_l] + 
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[u_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[u_uu];

    f_4 = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[uu_ll] + 
      ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[uu_l] + 
      (DY*(1.0+DY)*(2.0-DY)/2.0)*f[uu_u] + 
      (DY*(1.0+DY)*(DY-1.0)/6.0)*f[uu_uu];

    array_int = (DX*(DX-1.0)*(2.0-DX)/6.0)*f_1 + 
      ((1.0+DX)*(1.0-DX)*(2.0-DX)/2.0)*f_2 + 
      (DX*(1.0+DX)*(2.0-DX)/2.0)*f_3 + 
      (DX*(1.0+DX)*(DX-1.0)/6.0)*f_4;
  }
  else {
    cout << "error!" << endl;
    return 0.0;
  }
  return array_int;
}

//=============================================================
// Shift the xy-planes of the reference grid so that they line 
// up with those of the adm grid.
//=============================================================
void adjust_z_grid(int Nx, int Ny, int Nz, int *Nz_t, int nz_adm,
		   double dZ_adm, double dZ_ref, double symz, 
		   double *array_s, double ** array_t)
{
  double z_max = (Nz - 0.5) * dZ_ref;
  *Nz_t = floor(z_max/dZ_adm + 0.5) + 1.0;
  //int I_Print = ( comm_service::proc_me() == 0 );
  //  if(I_Print) cout << "Nz = " << Nz << ", Nz_t = " << *Nz_t 
  //		   << " adjust_z_grid:  Nx = " << Nx 
  //		   << ", dZ_ref = " << dZ_ref << endl;
  int nzt = floor(z_max/dZ_adm + 0.5) + 1.0;
  *array_t = new double[Nx*Ny*nzt];

  int i,j,k,tz,ty,tx,kl;
  double z,zl,zu,DZ;
  double *temp;
  temp = new double[Nx*Ny*nzt];
  for (i = 0; i<Nx*Ny*nzt; i++)
    {
      tz = floor(i/(Nx*Ny));
      ty = floor((i-tz*Nx*Ny)/Nx);
      tx = i - Nx*Ny*tz - Nx*ty;
      z = (tz + 0.5)*dZ_adm;
      kl = floor((z - 0.5*dZ_ref)/dZ_ref);
      if(kl>=Nz-2)
	kl = Nz-3;
      zl = (kl+0.5)*dZ_ref;
      zu = zl + dZ_ref;
      DZ = (z - zl)/(zu - zl);
      //
      // Note:  for cell-centered data:
      //                                f[-1] = f[0]
      //                                f[-2] = f[1], etc.
      //
      if(kl==0)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*array_s[Nx*ty + tx] + 
	  ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*array_s[Nx*ty + tx] + 
	  (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*array_s[Nx*Ny + Nx*ty + tx] + 
	  (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*array_s[Nx*Ny*2 + Nx*ty + tx];
      else if(kl==-1)
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*symz*array_s[Nx*Ny + Nx*ty + tx] + 
	  ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*symz*array_s[Nx*ty + tx] + 
	  (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*array_s[Nx*ty + tx] + 
	  (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*array_s[Nx*Ny + Nx*ty + tx];
      else
	temp[Nx*Ny*tz + Nx*ty + tx] = 
	  (DZ*(DZ-1.0)*(2.0-DZ)/6.0)*symz*array_s[Nx*Ny*(kl-1) + Nx*ty + tx] + 
	  ((1.0+DZ)*(1.0-DZ)*(2.0-DZ)/2.0)*array_s[Nx*Ny*kl + Nx*ty + tx] + 
	  (DZ*(1.0+DZ)*(2.0-DZ)/2.0)*array_s[Nx*Ny*(kl+1) + Nx*ty + tx] + 
	  (DZ*(1.0+DZ)*(DZ-1.0)/6.0)*array_s[Nx*Ny*(kl+2) + Nx*ty + tx];
    }
  for (i = 0; i<Nx*Ny*nzt; i++)  
    (*array_t)[i] = temp[i];
}

//-----------------------------------------------------------------------------
// Translate irrotational binary longitudinal terms.
//-----------------------------------------------------------------------------
/*
  void translate_irr(char *cctkGH,int *ex,int Nx,int Ny,int Nz,int nx_adm,int ny_adm,int ntot,
  double *Xsource,double *Ysource,double *Zsource,
  double *Xdest, double *Ydest,double *Zdest,
  double *shiftx_new,double * shifty_new,double *shiftz_new,
  double *v2_new,double *lapse_new,
  double *phi_new,
  double *ux_new,double *uy_new,double *uz_new)
  {
  //
  // Other variables:
  //
  int i,j,k,index;
  int imin, jmin, kmin, imax, jmax, kmax;
  int il,jl,kl,iu,ju,ku;
  double Courant,dX,dY,dZ,Delx,Dely,Delz,f_int,fac;
  double HALF, ONE, TWO;
  double symx, symy;
  // parameter(HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
  //
  // where the action is...
  //
  imin = 0;
  jmin = 0;
  kmin = 0;
  imax = ex[0];
  jmax = ex[1];
  kmax = ex[2];

  dX = Xsource[2] - Xsource[1];
  dY = Ysource[2] - Ysource[1];
  dZ = Zsource[2] - Zsource[1];
  printf("lasdf %e\t%e\n",dX,dZ);
  //
  //
  int lbz = 1; //cheap guess for EQUATORIAL symmetry
  for(int i=imin;i<imax;i++) 
  for(int j=jmin;j<jmax;j++) 
  for(int k=kmin;k<kmax;k++) {
  //index = ex(1)*ex(2)*(k-1) + ex(1)*(j-1) + i;
  index = CCTK_GFINDEX3D(cctkGH,i,j,k);
  weight(Xdest[index],Ydest[index],Zdest[index], 
  Xsource[1],Ysource[1],Zsource[1], 
  dX,dY,dZ,i,j,k,imin,imax,jmin,jmax,Nx,Ny, 
  il,iu,jl,ju,Delx,Dely);
  kl = k + lbz;
  //
  // Interpolate 
  //
  symx = 1.0;
  symy = -1.0;
  //printf("%d\t%d\t%d\t%d\t%d\t%e\t%e\t%e\n",ntot,Nx,Ny,Nz,nx_adm,Delx,Dely,shiftx_i[2]);
  //printf("bugga: %d\t%d\t%d\t%d\t%d\n",il,iu,jl,ju,kl);
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,shiftx_i,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy);
  //printf("iliublash: %d %d %d %d %d\n",il,iu,jl,ju,kl);
  //printf("hlkasjdeeeeeee %e %e %e %e %e\n",shiftx_i[il],shiftx_i[iu],shiftx_i[jl],shiftx_i[ju],shiftx_i[kl]);
  //printf("f_int: %d %d %d\t%.15f\t%e\n",i,j,k,f_int,Delx);
  //printf("xyz: %e %e %e\n",Xdest[index],Ydest[index],Zdest[index]);
  shiftx_new[index] = f_int;
  symx = -1.0;
  symy = 1.0;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,shifty_i,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy);
  shifty_new[index] = f_int;
  symx = -1.0;
  symy = -1.0;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,shiftz_i,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy);
  shiftz_new[index] = f_int;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,lapse_i,f_int,il,iu,jl,ju,kl,Delx,Dely,1.0,1.0);
  lapse_new[index] = f_int;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,phi_i,f_int,il,iu,jl,ju,kl,Delx,Dely,1.0,1.0);
  phi_new[index] = f_int;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,q_i,f_int,il,iu,jl,ju,kl,Delx,Dely,1.0,1.0);
  v2_new[index] = f_int;
  //if(f_int !=0) printf("xyz: %d %d %d\t %e\n",i,j,k,v2_new[index]);
  symx = 1.0;
  symy = -1.0;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,ux_i,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy);
  ux_new[index] = f_int;
  symx = -1.0;
  symy = 1.0;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,uy_i,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy);
  uy_new[index] = f_int;
  symx = -1.0;
  symy = -1.0;
  interp4(ntot,Nx,Ny,nx_adm,ny_adm,uz_i,f_int,il,iu,jl,ju,kl,Delx,Dely,symx,symy);
  uz_new[index] = f_int;
  }
  printf("v2_new = %.15f\n",v2_new[CCTK_GFINDEX3D(cctkGH,59,46,1)]);
  }
*/


//-----------------------------------------------------------------------------
//
// Interpolate function f using weights Delx and Dely
//
//-----------------------------------------------------------------------------
////
/*
  void interp4(int ntot,int Nx,int Ny,int nx_adm,int ny_adm,double *f,double &f_int,int il,
  int iu,int jl,int ju,int kl,double &DX,double &DY,double SymX,double SymY)
  {
  //
  // Other variables:
  //
  int iuc, ilc, juc, jlc, klc;
  int iuuc, illc, juuc, jllc;
  int ll_ll, ll_l, ll_u, ll_uu;
  int l_ll, l_l, l_u, l_uu;
  int u_ll, u_l, u_u, u_uu;
  int uu_ll, uu_l, uu_u, uu_uu;
  int Lx, Ly;
  double f_1x1y, f_2x1y, f_3x1y, f_4x1y;
  double f_1x2y, f_2x2y, f_3x2y, f_4x2y;
  double f_1x3y, f_2x3y, f_3x3y, f_4x3y;
  double f_1x4y, f_2x4y, f_3x4y, f_4x4y;
  double f_1x, f_2x, f_3x, f_4x;
  double f_int1, f_int2;
  double sxu,syu,sxuu,syuu,sxl,syl,sxll,syll;
  
  //
  // First find the appropriate indicies as they were assigned in C.
  //
  iuc  = iu - 1;
  iuuc = iuc + 1;
  ilc  = il - 1;
  illc = ilc - 1;
  juc  = ju - 1;
  juuc = juc + 1;
  jlc  = jl - 1;
  jllc = jlc - 1;
  klc  = kl - 1;
  Lx = Nx;
  Ly = Ny;
  //printf("mmmmmmmmmmm %d\t%d\n",Lx,Ly);
  //printf("dasdaieo: %d\t%d\t%d\t%d\t%d\t%d\t%d\n",iuc,iuuc,ilc,illc,juc,juuc,jlc);
  //
  // Now we reflect onto the appropriate quadrant
  //  
  if (iuc<Lx) {
  iuc = Lx - iuc - 1;
  sxu = SymX;
  } else {
  iuc = iuc - Lx;
  sxu = 1.0;
  }
  if (ilc<Lx) {
  ilc = Lx - ilc - 1;
  sxl = SymX;
  } else {
  ilc = ilc - Lx;
  sxl = 1.0;
  } 
  if (illc<Lx) {
  illc = Lx - illc - 1;
  sxll = SymX;
  } else {
  illc = illc - Lx;
  sxll = 1.0;
  }
  if (iuuc<Lx) {
  iuuc = Lx - iuuc - 1;
  sxuu = SymX;
  } else {
  iuuc = iuuc - Lx;
  sxuu = 1.0;
  }  
  if (juc<Ly) {
  juc = Ly - juc - 1;
  syu = SymY;
  } else {
  juc = juc - Ly;
  syu = 1.0;
  }
  if (jlc<Ly) {
  jlc = Ly - jlc - 1;
  syl = SymY;
  } else {
  jlc = jlc - Ly;
  syl = 1.0;
  } 
  if (jllc<Ly) {
  jllc = Ly - jllc - 1;
  syll = SymY;
  } else {
  jllc = jllc - Ly;
  syll = 1.0;
  }
  if (juuc<Ly) {
  juuc = Ly - juuc - 1;
  syuu = SymY;
  } else {
  juuc = juuc - Ly;
  syuu = 1.0;
  }
  //printf("aslfkkkkkkk: %e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",sxu,syu,sxuu,syuu,sxl,syl,sxll,syll);
  //
  // No point in interpolating in the vacuum
  //
  if (f[Ny*Nx*klc + Nx*jlc + ilc + 1]==0.0 
  && f[Ny*Nx*klc + Nx*jlc + iuc + 1]==0.0 
  && f[Ny*Nx*klc + Nx*juc + ilc + 1]==0.0 
  && f[Ny*Nx*klc + Nx*juc + iuc + 1]==0.0) {
  f_int = 0.0;
  }
  else {
  //
  //// Now for some 2D cubic interpolation
  //
  ll_ll  = Ny*Nx*klc + Nx*(jllc) + (illc) + 1;
  ll_l   = Ny*Nx*klc + Nx*jlc    + (illc) + 1;
  ll_u   = Ny*Nx*klc + Nx*juc    + (illc) + 1;
  ll_uu  = Ny*Nx*klc + Nx*(juuc) + (illc) + 1;
  l_ll   = Ny*Nx*klc + Nx*(jllc) + (ilc)  + 1;
  l_l    = Ny*Nx*klc + Nx*jlc    + (ilc)  + 1;
  l_u    = Ny*Nx*klc + Nx*juc    + (ilc)  + 1;
  l_uu   = Ny*Nx*klc + Nx*(juuc) + (ilc)  + 1;
  u_ll   = Ny*Nx*klc + Nx*(jllc) + (iuc)  + 1;
  u_l    = Ny*Nx*klc + Nx*jlc    + (iuc)  + 1;
  u_u    = Ny*Nx*klc + Nx*juc    + (iuc)  + 1;
  u_uu   = Ny*Nx*klc + Nx*(juuc) + (iuc)  + 1; 
  uu_ll  = Ny*Nx*klc + Nx*(jllc) + (iuuc) + 1;
  uu_l   = Ny*Nx*klc + Nx*jlc    + (iuuc) + 1;
  uu_u   = Ny*Nx*klc + Nx*juc    + (iuuc) + 1;
  uu_uu  = Ny*Nx*klc + Nx*(juuc) + (iuuc) + 1;      
  f_1x = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[ll_ll]*sxll*syll + 
  ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[ll_l]*sxll*syl + 
  (DY*(1.0+DY)*(2.0-DY)/2.0)*f[ll_u]*sxll*syu + 
  (DY*(1.0+DY)*(DY-1.0)/6.0)*f[ll_uu]*sxll*syuu;
  f_2x = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[l_ll]*sxl*syll + 
  ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[l_l]*sxl*syl + 
  (DY*(1.0+DY)*(2.0-DY)/2.0)*f[l_u]*sxl*syu + 
  (DY*(1.0+DY)*(DY-1.0)/6.0)*f[l_uu]*sxl*syuu;
  f_3x = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[u_ll]*sxu*syll + 
  ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[u_l]*sxu*syl + 
  (DY*(1.0+DY)*(2.0-DY)/2.0)*f[u_u]*sxu*syu + 
  (DY*(1.0+DY)*(DY-1.0)/6.0)*f[u_uu]*sxu*syuu;
  f_4x = (DY*(DY-1.0)*(2.0-DY)/6.0)*f[uu_ll]*sxuu*syll + 
  ((1.0+DY)*(1.0-DY)*(2.0-DY)/2.0)*f[uu_l]*sxuu*syl + 
  (DY*(1.0+DY)*(2.0-DY)/2.0)*f[uu_u]*sxuu*syu + 
  (DY*(1.0+DY)*(DY-1.0)/6.0)*f[uu_uu]*sxuu*syuu;
  //    if(f_4x==0) {
  //printf("blah %e %e %e %e %e\n",DY,f_1x,f_2x,f_3x,f_4x);
  //printf("\nblah2 %e %e %e %e\n\n",f[uu_ll],f[uu_l],f[uu_u],f[uu_uu]);
  //printf("blah3 %d %d %d %d\n",uu_ll,uu_l,uu_u,uu_uu);
  //printf("dasdaieo: %d\t%d\t%d\t%d\t%d\t%d\t%d\n",iuc,iuuc,ilc,illc,juc,juuc,jlc);
  //printf("blah3 %e %e\n",sxl,syl);
  }
  //printf("blah %e %e %e %e %e %e %e %e\n",sxuu,syll,sxuu,syl,sxuu,syu,sxuu,syuu);
  f_int1 = (DX*(DX-1.0)*(2.0-DX)/6.0)*f_1x + 
  ((1.0+DX)*(1.0-DX)*(2.0-DX)/2.0)*f_2x + 
  (DX*(1.0+DX)*(2.0-DX)/2.0)*f_3x + 
  (DX*(1.0+DX)*(DX-1.0)/6.0)*f_4x;
  f_int = f_int1;
  }
  }

 
  //
  //-----------------------------------------------------------------------------
  //
  // Interpolation for point i,j,k
  //
  //-----------------------------------------------------------------------------
  //
  void weight(double Xdest,double Ydest,double Zdest, 
  double Xgmin,double Ygmin,double Zgmin, 
  double dX,double dY,double dZ,int i,int j,int k, 
  int imin,int imax,int jmin,int jmax,int Nx,int Ny, 
  int &il,int &iu,int &jl,int &ju,double &Delx,double &Dely) {
  double xl,xu,yl,yu,ONE,TWO;
  ONE = 2.0;
  TWO = 2.0;
  
  il = floor((Xdest-Xgmin)/dX) + 1;
  jl = floor((Ydest-Ygmin)/dY) + 1;
  if (il < 2) {
  il = 2;
  }
  if (jl < 2) {
  jl = 2;
  }
  if (il > 2*Nx-1) {
  il = 2*Nx-1;
  }
  if (jl > 2*Ny-1) {
  jl = 2*Ny-1;
  }
  //
  iu = il + 1;
  ju = jl + 1;
  xl = Xgmin + dX*(il-1) ;
  xu = xl + dX;
  yl = Ygmin + dY*(jl - 1);
  yu = yl + dY;
  //
  // get 'upper' values
  //
  // find coefficients for interpolation
  //
  Delx = (Xdest - xl)/(xu - xl);
  Dely = (Ydest - yl)/(yu - yl);
  
  }
*/

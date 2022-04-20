//----------------------------------------------------------------------------------
// Advection routines for hydro variables, valid for arbitrary numbers of ghostzones
//----------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define f1o3 0.3333333333333333333333333333333333333333
#define f1o6 0.1666666666666666666666666666666666666666
// Set B_LIMITER_COEFF = 2.0  for MC, 1 for minmod
#define B_LIMITER_COEFF 2.0

//FORTRAN-callable wrappers
extern "C" void CCTK_FCALL compute_du_everywhere_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,int *Symmetry,
   int *cell_centering_enabled,int *U_syms,double *U,double *dU);

extern "C" void CCTK_FCALL mc_slope_limiter_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,int *Symmetry,
   double *dU,double *nabla_U);

extern "C" void CCTK_FCALL compute_ur_ul_general_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,
   double *X,double *Y,double *Z,int *Symmetry,int *cell_centering_enabled,
   int *U_syms,double *U,double *nabla_U,double *Ur,double *Ul);

extern "C" void CCTK_FCALL compute_ur_ul_sppm_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,
   double *X,double *Y,double *Z,int *Symmetry,int *cell_centering_enabled,
   int *U_syms,double *U,double *dU,double *Ur,double *Ul);

//Compute dU = U(i) - U(i-1) everywhere on grid, where i is an array index whose spatial direction is given by flux_direction
// Since dU is computed everywhere on grid, we need U to be set EVERYWHERE on the grid as well (even ghostzones!)
extern "C" void compute_dU_everywhere(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
				      int cell_centering_enabled,int *U_syms,double *U,double *dU) {
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  //Stencil width below is 1.
  if(flux_direction==1) istart++;
  else if(flux_direction==2) jstart++;
  else if(flux_direction==3) kstart++;

  //Only update y=0 plane in axisymmetry!
  if(Symmetry==4) { jstart=1; jend=2; }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int indexm1,index;
    // Yes, it is bad form to put a bunch of if() statements inside a loop like this, 
    //      but I'm hoping the compiler is smart enough to handle them efficiently.
    if(flux_direction==1) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
      index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
    } else if(flux_direction==2) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
      index   = CCTK_GFINDEX3D(cctkGH,i,j  ,k);
    } else if(flux_direction==3){
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
      index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
    }
    //if(isnan(U[index]) || isnan(U[indexm1])) {
    //  printf("FOUND NAN AT %d %d %d, direction %d\n",i,j,k,flux_direction);
    //  exit(1);
    //}
    dU[index] = U[index] - U[indexm1];
  }

  //Next fill in imin,jmin,kmin values for dU, accounting for arbitrary numbers of ghostzones.
  // Note that if imin,jmin,kmin is NOT a symmetry boundary, these values are really not correct; 
  //   you should set symi = 1 in this case (i.e., treat U as a nearly flat function near the outer boundary).
  if(flux_direction==1) {
    //x-direction:
    int nsymghostx;
    if(Symmetry==4) {
      //In axisymmetry, nghostzones=2, though [# symmetry ghostzones] = 1.
      nsymghostx=1;
    } else {
      nsymghostx=nghostzones[0];
    }
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) {
      int index   = CCTK_GFINDEX3D(cctkGH,0,j,k);
      int indexm1 = CCTK_GFINDEX3D(cctkGH,nsymghostx*2,j,k);
      if(cell_centering_enabled==0) indexm1 = CCTK_GFINDEX3D(cctkGH,nsymghostx*2+1,j,k);
      if (Symmetry==2 || Symmetry==4) {
	dU[index] = U[index] - U_syms[0]*U[indexm1];
      } else {
	dU[index] = dU[CCTK_GFINDEX3D(cctkGH,1,j,k)];
      }
    }
  } else if(flux_direction==2) {
    //y-direction:
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int i=istart;i<iend;i++) {
      int index   = CCTK_GFINDEX3D(cctkGH,i,0,k);
      int indexm1 = CCTK_GFINDEX3D(cctkGH,i,nghostzones[1]*2,k);
      if(cell_centering_enabled==0) indexm1 = CCTK_GFINDEX3D(cctkGH,i,nghostzones[1]*2+1,k);
      if (Symmetry==2) {
	dU[index] = U[index] - U_syms[1]*U[indexm1];
      } else { 
	dU[index] = dU[CCTK_GFINDEX3D(cctkGH,i,1,k)];
      }
    }
  } else if(flux_direction==3){
    //z-direction:
#pragma omp parallel for
    for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index   = CCTK_GFINDEX3D(cctkGH,i,j,0);
      int indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,nghostzones[2]*2);
      if(cell_centering_enabled==0) indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,nghostzones[2]*2+1);
      if (Symmetry != 0) {
	dU[index] = U[index] - U_syms[2]*U[indexm1];
      } else {
	dU[index] = dU[CCTK_GFINDEX3D(cctkGH,i,j,1)];
      } 
    }
  }
}

extern "C" void CCTK_FCALL compute_du_everywhere_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,int *Symmetry,
   int *cell_centering_enabled,int *U_syms,double *U,double *dU)
{
  compute_dU_everywhere(*flux_direction,*cctkGH,cctk_lsh,nghostzones,*Symmetry,*cell_centering_enabled,U_syms,U,dU);
}

//Apply MC slope limiter to dU = U(i) - U(i-1), where the direction of i is given by flux_direction.
extern "C" void mc_slope_limiter(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
				 double *dU,double *nabla_U) {
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  //Stencil width below is 1.
  if(flux_direction==1) iend--;
  else if(flux_direction==2) jend--;
  else if(flux_direction==3) kend--;

  //Only update y=0 plane in axisymmetry!
  if(Symmetry==4) { jstart=1; jend=2; }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int indexm1,index,indexp1;
    // Yes, it is bad form to put a bunch of if() statements inside a loop like this, 
    //      but I'm hoping the compiler is smart enough to handle them efficiently.
    if(flux_direction==1) {
      index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
    } else if(flux_direction==2) {
      index   = CCTK_GFINDEX3D(cctkGH,i,j  ,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
    } else if(flux_direction==3){
      index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
    }

    double dUL     = dU[index];
    double dUL_ip1 = dU[indexp1];

    //Eq. 60 in JOURNAL OF COMPUTATIONAL PHYSICS 123, 1–14 (1996) [note the factor of 2 missing in the |a_{j+1} - a_{j}| term]. Recall that dU = U_{i} - U_{i-1}.
    if(dUL*dUL_ip1 > 0.0) {
      //m_dU is the average (u_(i+1) - u_(i-1))/2 -- first derivative, second-order; this should happen most of the time (smooth flows)
      double m_dU = 0.5*(dUL + dUL_ip1);
      //Decide whether to use 2nd order derivative or first-order derivative, limiting slope.
      if ( (fabs(dUL) <= fabs(dUL_ip1)) && (B_LIMITER_COEFF*fabs(dUL) < fabs(m_dU)) ) {
	nabla_U[index] = B_LIMITER_COEFF*dUL; // B_LIMITER_COEFF times the smallest amplitude 1st-order slope
      } else if ( (fabs(dUL_ip1) <= fabs(dUL)) && (B_LIMITER_COEFF*fabs(dUL_ip1) < fabs(m_dU)) ) {
	nabla_U[index] = B_LIMITER_COEFF*dUL_ip1; // B_LIMITER_COEFF times the smallest amplitude 1st-order slope
      } else {
	nabla_U[index] = m_dU;
      }
    } else {
      nabla_U[index] = 0.0;
    }
  }

  //Ur and Ul depend on nabla_U being set everywhere on grid, so we'd better set nabla_U everywhere!
  if(flux_direction==1) {
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) {
      int index = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,j,k);
      nabla_U[index] = dU[index];
    }
  } else if(flux_direction==2) {
    //y-direction:
    for(int k=kstart;k<kend;k++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,cctk_lsh[1]-1,k);
      nabla_U[index] = dU[index];
    }
  } else if(flux_direction==3){
    //z-direction:
    for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,cctk_lsh[2]-1);
      nabla_U[index] = dU[index];
    }
  }
  //  printf("HIASDHFOIASDF: %.15e\n",nabla_U[CCTK_GFINDEX3D(cctkGH,53,1,0)]);
}
extern "C" void CCTK_FCALL mc_slope_limiter_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,int *Symmetry,
   double *dU,double *nabla_U)
{
  mc_slope_limiter(*flux_direction,*cctkGH,cctk_lsh,nghostzones,*Symmetry,dU,nabla_U);
}

//Set Ur and Ul:
extern "C" void compute_Ur_Ul_general(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,
				      double *X,double *Y,double *Z,int Symmetry,int cell_centering_enabled,
				      int *U_syms,double *U,double *nabla_U,double *Ur,double *Ul) {
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  //Stencil width below is 1.
  if(flux_direction==1) istart++;
  else if(flux_direction==2) jstart++;
  else if(flux_direction==3) kstart++;

  //Only update y=0 plane in axisymmetry!
  if(Symmetry==4) { jstart=1; jend=2; }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    // Yes, it is bad form to put a bunch of if() statements inside a loop like this, 
    //      but I'm hoping the compiler is smart enough to handle them efficiently.
    int index,indexm1;
    if(flux_direction==1) {
      index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
    } else if(flux_direction==2) {
      index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
    } else if(flux_direction==3){
      index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
    }

    Ur[index] = U[index]   - 0.5*nabla_U[index];
    Ul[index] = U[indexm1] + 0.5*nabla_U[indexm1];
  }

  //Next fill in i=0,j=0,k=0 values for dU, handling arbitrary numbers of ghostzones.
  //Gotta be careful here!  
  //  Note that vxr(2,:,:) = vxr(x=0+,:,:) and vxl(2,:,:) = vxl(x=0-,:,:), so 
  //             when x is a symmetry axis, vxl(2,:,:) = [Symmetry argument] * vxr(2,:,:)!
  // Also, note that if imin,jmin,kmin is NOT a symmetry boundary, these values are really not correct; 
  //   you should set symi = 1 in this case (i.e., treat U as a nearly flat function near the outer boundary).
  if(flux_direction==1) {
    //x-direction:
    int nsymghostx=nghostzones[0];
    //In axisymmetry, nghostzones=2, though [# symmetry ghostzones] = 1.
    if(Symmetry==4) nsymghostx=1;
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) {
      int index = CCTK_GFINDEX3D(cctkGH,0,j,k);
      int indexsym = CCTK_GFINDEX3D(cctkGH,nsymghostx*2,j,k);
      if(cell_centering_enabled==0) indexsym = CCTK_GFINDEX3D(cctkGH,nsymghostx*2+1,j,k);
      if(indexsym!=index) Ur[index] = U_syms[0]*Ul[indexsym];
      else                Ur[index] = U[index] - 0.5*nabla_U[index];
      Ul[index] = U_syms[0]*Ur[indexsym];
      if (Symmetry==0) Ul[index] = Ur[index];
      
      if (isnan(Ul[index])){
	printf("x-direction!!! \n");
	  printf("Ul is NaN!!!!, i=0,j=%d,k=%d \n", j,k);
        printf("Ur=%e, U=%e, nabla_U=%e \n", Ur[index], U[index], nabla_U[index]);
	printf("Ul[indexsym]=%e,Ur[indexsym]=%e", Ul[indexsym],Ur[indexsym]);
	  }
      }
  } else if(flux_direction==2) {
    //y-direction (flux_direction!=2 in AXISYM!):
    for(int k=kstart;k<kend;k++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,0,k);
      int indexsym = CCTK_GFINDEX3D(cctkGH,i,nghostzones[1]*2,k);
      if(cell_centering_enabled==0) indexsym = CCTK_GFINDEX3D(cctkGH,i,nghostzones[1]*2+1,k);
      if(indexsym!=index) Ur[index] = U_syms[1]*Ul[indexsym];
      else                Ur[index] = U[index] - 0.5*nabla_U[index];
      Ul[index] = U_syms[1]*Ur[indexsym];
      if (Symmetry==0) Ul[index] = Ur[index];

      if (isnan(Ul[index])){
	printf("y-direction!!! \n");
	  printf("Ul is NaN!!!!, i=%d,j=0,k=%d \n", i,k);
        printf("Ur=%e, U=%e, nabla_U=%e \n", Ur[index], U[index], nabla_U[index]);
	printf("Ul[indexsym]=%e,Ur[indexsym]=%e", Ul[indexsym],Ur[indexsym]);
	  }
    }
  } else if(flux_direction==3){
    //z-direction:
    for(int j=jstart;j<jend;j++)  for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,0);
      //FIXME: AXISYMMETRY UNSUPPORTED:
      int indexsym = CCTK_GFINDEX3D(cctkGH,i,j,nghostzones[2]*2);
      if(cell_centering_enabled==0) indexsym = CCTK_GFINDEX3D(cctkGH,i,j,nghostzones[2]*2+1);
      if(indexsym!=index) Ur[index] = U_syms[2]*Ul[indexsym];
      else                Ur[index] = U[index] - 0.5*nabla_U[index];
      Ul[index] = U_syms[2]*Ur[indexsym];
      if (Symmetry==0) Ul[index] = Ur[index];
      
      if (isnan(Ul[index])){
	printf("z-direction!!! \n");
	printf("Ul is NaN!!!!, i=%d,j=%d,k=0 \n", i,j);
	printf("Ur=%e, U=%e, nabla_U=%e \n", Ur[index], U[index], nabla_U[index]);
	printf("Ul[indexsym]=%e,Ur[indexsym]=%e", Ul[indexsym],Ur[indexsym]);
      }

      }
  }

  
}



extern "C" void CCTK_FCALL compute_ur_ul_general_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,
   double *X,double *Y,double *Z,int *Symmetry,int *cell_centering_enabled,
   int *U_syms,double *U,double *nabla_U,double *Ur,double *Ul)
{
  compute_Ur_Ul_general(*flux_direction,*cctkGH,cctk_lsh,nghostzones,X,Y,Z,*Symmetry,*cell_centering_enabled,U_syms,U,nabla_U,Ur,Ul);
}

extern "C" void compute_Ur_Ul_SPPM(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,
				   double *X,double *Y,double *Z,int Symmetry,int cell_centering_enabled,
				   int *U_syms,double *U,double *dU,double *Ur,double *Ul) {
  // Shibata's b parameter:
  double b = 3.0;

  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  //Stencil width at outer bdry is 1.
  if(flux_direction==1)      iend--;
  else if(flux_direction==2) jend--;
  else if(flux_direction==3) kend--;

  //Only update y=0 plane in axisymmetry!
  if(Symmetry==4) { jstart=1; jend=2; }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    // Yes, it is bad form to put a bunch of if() statements inside a loop like this, 
    //      but I'm hoping the compiler is smart enough to handle them efficiently.
    int index,indexp1;
    if(flux_direction==1) {
      index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
    } else if(flux_direction==2) {
      index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
    } else if(flux_direction==3){
      index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
    }
    
    double dpU=0.0;
    double dmU=0.0;
   
 
    if(dU[index]*dU[indexp1] > 0.0) {
      double dUl = dU[index];
      double dUlp1 = dU[indexp1];

      double dUfrac = b*dUlp1/dUl;
      if(dUl==0.0) dUfrac = b*dUlp1/1e-100;

      double dUp1frac = b*dUl/dUlp1;
      if(dUlp1==0.0) dUfrac = b*dUl/1e-100;
      
      dpU = dUl;
      dmU = dUlp1;

      if(dUfrac<1.0) dpU = b*dUlp1;
      if(dUp1frac<1.0) dmU = b*dUl;
    }

    Ur[index] = U[index] - f1o3*dpU - f1o6*dmU;
    Ul[index] = U[index] + f1o6*dpU + f1o3*dmU;
  }

  // Now we fill in the boundary points "i".  
  //    We reuse the "i-1" point values of dpU and dmU.
  if(flux_direction==1) {
    int index,indexm1;
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) {
      index   = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,j,k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-2,j,k);

      Ur[index] = (Ur[indexm1] - U[indexm1]) + U[index];
      Ul[index] = (Ul[indexm1] - U[indexm1]) + U[index];
    }
  } else if(flux_direction==2) {
    int index,indexm1;
    for(int k=kstart;k<kend;k++) for(int i=istart;i<iend;i++) {
      index   = CCTK_GFINDEX3D(cctkGH,i,cctk_lsh[1]-1,k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,cctk_lsh[1]-2,k);

      Ur[index] = (Ur[indexm1] - U[indexm1]) + U[index];
      Ul[index] = (Ul[indexm1] - U[indexm1]) + U[index];
    }
  } else if(flux_direction==3) {
    int index,indexm1;
    for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      index   = CCTK_GFINDEX3D(cctkGH,i,j,cctk_lsh[2]-1);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,cctk_lsh[2]-2);

      Ur[index] = (Ur[indexm1] - U[indexm1]) + U[index];
      Ul[index] = (Ul[indexm1] - U[indexm1]) + U[index];
    }
  }

}

extern "C" void CCTK_FCALL compute_ur_ul_sppm_
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,
   double *X,double *Y,double *Z,int *Symmetry,int *cell_centering_enabled,
   int *U_syms,double *U,double *dU,double *Ur,double *Ul)
{
  compute_Ur_Ul_SPPM(*flux_direction,*cctkGH,cctk_lsh,nghostzones,X,Y,Z,*Symmetry,*cell_centering_enabled,U_syms,U,dU,Ur,Ul);
}

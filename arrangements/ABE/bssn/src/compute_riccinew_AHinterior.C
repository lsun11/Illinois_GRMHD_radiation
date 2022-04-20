/*
  compute Ricci tensor and Hamiltonian+momentum constraints
  Optimized from FORTRAN version by Zach Etienne: Apr, 2007
  (3x speed increase with better scalability)
  Speed increase due primarily to
  1) Fewer calls to main memory 
  -> can do calculations entirely in (uber-fast) processor cache!
  2) Replace derivative function calls (eww) 
  with "#define" finite differencing (yay).
  3) No more allocating entire temporary grid functions!
*/

#define KRANC_C
#define FD_SET_BY_USER
#define FD_C2

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

#define F2o3 0.666666666666666666666666666666666

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "GenericFD.h"

static char *rcsid="$mew. $";
CCTK_FILEVERSION(BSSN_ricci_AHinterior)


  extern "C" void CCTK_FCALL CCTK_FNAME(BSSN_ricci_AHinterior)
  (const cGH **cctkGH,double *dT, double *dx, double *dy, double *dz,
   int *nghostzones,int *cctk_lsh,double *r,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,
   double *trRtilde,
   double *Gammax,double *Gammay,double *Gammaz,
   double *gxxx, double *gxxy, double *gxxz,
   double *gxyx, double *gxyy, double *gxyz,
   double *gxzx, double *gxzy, double *gxzz,
   double *gyyx, double *gyyy, double *gyyz,
   double *gyzx, double *gyzy, double *gyzz,
   double *gzzx, double *gzzy, double *gzzz,
   double *Gamxxx, double *Gamxxy, double *Gamxxz, double *Gamxyy, double *Gamxyz, double *Gamxzz,
   double *Gamyxx, double *Gamyxy, double *Gamyxz, double *Gamyyy, double *Gamyyz, double *Gamyzz,
   double *Gamzxx, double *Gamzxy, double *Gamzxz, double *Gamzyy, double *Gamzyz, double *Gamzzz,
   double *Sx, double *Sy, double *Sz,
   double *Aupxx,double *Aupxy,double *Aupxz,double *Aupyy,double *Aupyz,double *Aupzz,
   double *phi, double *trK, double *MResx,double *MResy,double *MResz,double *MNorm, 
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *psi, double *rho, double *PsiRes, double *PsiNorm,int *compute_constraint_flag);
  
extern "C" void BSSN_ricci_AHinterior(const cGH *cctkGH, double dT, double dx, double dy, double dz,
						      int *nghostzones,int *cctk_lsh,double *r,
					   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
					   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
					   double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,
					   double *trRtilde,
					   double *Gammax,double *Gammay,double *Gammaz,
					   double *gxxx, double *gxxy, double *gxxz,
					   double *gxyx, double *gxyy, double *gxyz,
					   double *gxzx, double *gxzy, double *gxzz,
					   double *gyyx, double *gyyy, double *gyyz,
					   double *gyzx, double *gyzy, double *gyzz,
					   double *gzzx, double *gzzy, double *gzzz,
					   double *Gamxxx, double *Gamxxy, double *Gamxxz, double *Gamxyy, double *Gamxyz, double *Gamxzz,
					   double *Gamyxx, double *Gamyxy, double *Gamyxz, double *Gamyyy, double *Gamyyz, double *Gamyzz,
					   double *Gamzxx, double *Gamzxy, double *Gamzxz, double *Gamzyy, double *Gamzyz, double *Gamzzz,
					   double *Sx, double *Sy, double *Sz,
					   double *Aupxx,double *Aupxy,double *Aupxz,double *Aupyy,double *Aupyz,double *Aupzz,
					   double *phi, double *trK, double *MResx,double *MResy,double *MResz,double *MNorm, 
					   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
					   double *psi, double *rho, double *PsiRes, double *PsiNorm,int compute_constraint_flag) {
 
  DECLARE_CCTK_PARAMETERS;

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = nghostzones[0];
  int jstart = nghostzones[1];
  int kstart = nghostzones[2];
  int iend = cctk_lsh[0] - nghostzones[0];
  int jend = cctk_lsh[1] - nghostzones[1];
  int kend = cctk_lsh[2] - nghostzones[2];

  //Following lines needed since nghostzones[0] = ORDER, and
  //   not ORDER-1 in axisymmetry
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }

#pragma omp parallel 
  {
    // MAIN LOOP: Computes Ricci & constraints_AHinterior (if desired).
    //#pragma omp parallel for
#pragma omp for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

      if(r[index]<0.25) {
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];

      double GammaxL,GammayL,GammazL;
      if(Gammai_enforcement_type!=2) {
	GammaxL = Gammax[index];
	GammayL = Gammay[index];
	GammazL = Gammaz[index];
      } else {
	// Gamma^i constraint enforcing a la arXiv:gr-qc/0701123 (Marronetti et al.) and others
	GammaxL = -(D1gf(gupxx,i,j,k) + D2gf(gupxy,i,j,k) + D3gf(gupxz,i,j,k));
	GammayL = -(D1gf(gupxy,i,j,k) + D2gf(gupyy,i,j,k) + D3gf(gupyz,i,j,k));
	GammazL = -(D1gf(gupxz,i,j,k) + D2gf(gupyz,i,j,k) + D3gf(gupzz,i,j,k));
      }

      double gxxxL = D1gf(gxx,i,j,k);
      double gxyxL = D1gf(gxy,i,j,k);
      double gxzxL = D1gf(gxz,i,j,k);
      double gyyxL = D1gf(gyy,i,j,k);
      double gyzxL = D1gf(gyz,i,j,k);
      double gzzxL = D1gf(gzz,i,j,k);

      double gxxyL = D2gf(gxx,i,j,k);
      double gxyyL = D2gf(gxy,i,j,k);
      double gxzyL = D2gf(gxz,i,j,k);
      double gyyyL = D2gf(gyy,i,j,k);
      double gyzyL = D2gf(gyz,i,j,k);
      double gzzyL = D2gf(gzz,i,j,k);

      double gxxzL = D3gf(gxx,i,j,k);
      double gxyzL = D3gf(gxy,i,j,k);
      double gxzzL = D3gf(gxz,i,j,k);
      double gyyzL = D3gf(gyy,i,j,k);
      double gyzzL = D3gf(gyz,i,j,k);
      double gzzzL = D3gf(gzz,i,j,k);

      double Gamxx = D1gf(Gammax,i,j,k);
      double Gamyx = D1gf(Gammay,i,j,k);
      double Gamzx = D1gf(Gammaz,i,j,k);

      double Gamxy = D2gf(Gammax,i,j,k);
      double Gamyy = D2gf(Gammay,i,j,k);
      double Gamzy = D2gf(Gammaz,i,j,k);

      double Gamxz = D3gf(Gammax,i,j,k);
      double Gamyz = D3gf(Gammay,i,j,k);
      double Gamzz = D3gf(Gammaz,i,j,k);
	
      double gupxxL = gupxx[index];
      double gupxyL = gupxy[index];
      double gupxzL = gupxz[index];
      double gupyyL = gupyy[index];
      double gupyzL = gupyz[index];
      double gupzzL = gupzz[index];

      double GamxxxL = 0.5*(gupxxL *gxxxL +gupxyL *(2.0*gxyxL -gxxyL )+gupxzL *(2.0*gxzxL -gxxzL ));
      double GamyxxL = 0.5*(gupxyL *gxxxL +gupyyL *(2.0*gxyxL -gxxyL )+gupyzL *(2.0*gxzxL -gxxzL ));
      double GamzxxL = 0.5*(gupxzL *gxxxL +gupyzL *(2.0*gxyxL -gxxyL )+gupzzL *(2.0*gxzxL -gxxzL ));          

      double GamxyyL = 0.5*(gupxxL *(2.0*gxyyL -gyyxL )+gupxyL *gyyyL +gupxzL *(2.0*gyzyL -gyyzL ));
      double GamyyyL = 0.5*(gupxyL *(2.0*gxyyL -gyyxL )+gupyyL *gyyyL +gupyzL *(2.0*gyzyL -gyyzL ));
      double GamzyyL = 0.5*(gupxzL *(2.0*gxyyL -gyyxL )+gupyzL *gyyyL +gupzzL *(2.0*gyzyL -gyyzL ));

      double GamxzzL = 0.5*(gupxxL *(2.0*gxzzL -gzzxL )+gupxyL *(2.0*gyzzL -gzzyL )+gupxzL *gzzzL );
      double GamyzzL = 0.5*(gupxyL *(2.0*gxzzL -gzzxL )+gupyyL *(2.0*gyzzL -gzzyL )+gupyzL *gzzzL );
      double GamzzzL = 0.5*(gupxzL *(2.0*gxzzL -gzzxL )+gupyzL *(2.0*gyzzL -gzzyL )+gupzzL *gzzzL );

      double GamxxyL = 0.5*( gupxxL * gxxyL + gupxyL * gyyxL + gupxzL *( gxzyL + gyzxL - gxyzL ) );
      double GamyxyL = 0.5*( gupxyL * gxxyL + gupyyL * gyyxL + gupyzL *( gxzyL + gyzxL - gxyzL ) );
      double GamzxyL = 0.5*( gupxzL * gxxyL + gupyzL * gyyxL + gupzzL *( gxzyL + gyzxL - gxyzL ) );

      double GamxxzL = 0.5*( gupxxL * gxxzL + gupxyL *( gxyzL + gyzxL - gxzyL ) + gupxzL * gzzxL );
      double GamyxzL = 0.5*( gupxyL * gxxzL + gupyyL *( gxyzL + gyzxL - gxzyL ) + gupyzL * gzzxL );
      double GamzxzL = 0.5*( gupxzL * gxxzL + gupyzL *( gxyzL + gyzxL - gxzyL ) + gupzzL * gzzxL );

      double GamxyzL = 0.5*( gupxxL *( gxyzL + gxzyL - gyzxL ) + gupxyL * gyyzL + gupxzL * gzzyL );
      double GamyyzL = 0.5*( gupxyL *( gxyzL + gxzyL - gyzxL ) + gupyyL * gyyzL + gupyzL * gzzyL );
      double GamzyzL = 0.5*( gupxzL *( gxyzL + gxzyL - gyzxL ) + gupyzL * gyyzL + gupzzL * gzzyL );

      //-------------------
      //Next compute Ricci:
      //-------------------
      //

      // Note that we store g^{lm} g_{ij,lm} into R_{ij}
      double gxxxx = D11gf(gxx,i,j,k);
      double gxxxy = D21gf(gxx,i,j,k);
      double gxxxz = D31gf(gxx,i,j,k);
      double gxxyy = D22gf(gxx,i,j,k);
      double gxxyz = D32gf(gxx,i,j,k);
      double gxxzz = D33gf(gxx,i,j,k);
      double RxxL = gupxxL * gxxxx + gupyyL * gxxyy + gupzzL * gxxzz +
	2.0 * ( gupxyL * gxxxy + gupxzL * gxxxz + gupyzL * gxxyz );

      double gxyxx = D11gf(gxy,i,j,k);
      double gxyxy = D21gf(gxy,i,j,k);
      double gxyxz = D31gf(gxy,i,j,k);
      double gxyyy = D22gf(gxy,i,j,k);
      double gxyyz = D32gf(gxy,i,j,k);
      double gxyzz = D33gf(gxy,i,j,k);
      double RxyL = gupxxL * gxyxx + gupyyL * gxyyy + gupzzL * gxyzz +
	2.0 * ( gupxyL * gxyxy + gupxzL * gxyxz + gupyzL * gxyyz );

      double gxzxx = D11gf(gxz,i,j,k);
      double gxzxy = D21gf(gxz,i,j,k);
      double gxzxz = D31gf(gxz,i,j,k);
      double gxzyy = D22gf(gxz,i,j,k);
      double gxzyz = D32gf(gxz,i,j,k);
      double gxzzz = D33gf(gxz,i,j,k);
      double RxzL = gupxxL * gxzxx + gupyyL * gxzyy + gupzzL * gxzzz +
	2.0 * ( gupxyL * gxzxy + gupxzL * gxzxz + gupyzL * gxzyz );

      double gyyxx = D11gf(gyy,i,j,k);
      double gyyxy = D21gf(gyy,i,j,k);
      double gyyxz = D31gf(gyy,i,j,k);
      double gyyyy = D22gf(gyy,i,j,k);
      double gyyyz = D32gf(gyy,i,j,k);
      double gyyzz = D33gf(gyy,i,j,k);
      double RyyL = gupxxL * gyyxx + gupyyL * gyyyy + gupzzL * gyyzz +
	2.0 * ( gupxyL * gyyxy + gupxzL * gyyxz + gupyzL * gyyyz );
	
      double gyzxx = D11gf(gyz,i,j,k);
      double gyzxy = D21gf(gyz,i,j,k);
      double gyzxz = D31gf(gyz,i,j,k);
      double gyzyy = D22gf(gyz,i,j,k);
      double gyzyz = D32gf(gyz,i,j,k);
      double gyzzz = D33gf(gyz,i,j,k);
      double RyzL = gupxxL * gyzxx + gupyyL * gyzyy + gupzzL * gyzzz +
	2.0 * ( gupxyL * gyzxy + gupxzL * gyzxz + gupyzL * gyzyz );

      double gzzxx = D11gf(gzz,i,j,k);
      double gzzxy = D21gf(gzz,i,j,k);
      double gzzxz = D31gf(gzz,i,j,k);
      double gzzyy = D22gf(gzz,i,j,k);
      double gzzyz = D32gf(gzz,i,j,k);
      double gzzzz = D33gf(gzz,i,j,k);
      double RzzL = gupxxL * gzzxx + gupyyL * gzzyy + gupzzL * gzzzz +
	2.0 * ( gupxyL * gzzxy + gupxzL * gzzxz + gupyzL * gzzyz );

      double ass_Gamxxx =  gxxL* GamxxxL+ gxyL* GamyxxL+ gxzL* GamzxxL;
      double ass_Gamxxy =  gxxL* GamxxyL+ gxyL* GamyxyL+ gxzL* GamzxyL;
      double ass_Gamxxz =  gxxL* GamxxzL+ gxyL* GamyxzL+ gxzL* GamzxzL;
      double ass_Gamxyy =  gxxL* GamxyyL+ gxyL* GamyyyL+ gxzL* GamzyyL;
      double ass_Gamxyz =  gxxL* GamxyzL+ gxyL* GamyyzL+ gxzL* GamzyzL;
      double ass_Gamxzz =  gxxL* GamxzzL+ gxyL* GamyzzL+ gxzL* GamzzzL;
	
      double ass_Gamyxx =  gxyL* GamxxxL+ gyyL* GamyxxL+ gyzL* GamzxxL;
      double ass_Gamyxy =  gxyL* GamxxyL+ gyyL* GamyxyL+ gyzL* GamzxyL;
      double ass_Gamyxz =  gxyL* GamxxzL+ gyyL* GamyxzL+ gyzL* GamzxzL;
      double ass_Gamyyy =  gxyL* GamxyyL+ gyyL* GamyyyL+ gyzL* GamzyyL;
      double ass_Gamyyz =  gxyL* GamxyzL+ gyyL* GamyyzL+ gyzL* GamzyzL;
      double ass_Gamyzz =  gxyL* GamxzzL+ gyyL* GamyzzL+ gyzL* GamzzzL;
	
      double ass_Gamzxx =  gxzL* GamxxxL+ gyzL* GamyxxL+ gzzL* GamzxxL;
      double ass_Gamzxy =  gxzL* GamxxyL+ gyzL* GamyxyL+ gzzL* GamzxyL;
      double ass_Gamzxz =  gxzL* GamxxzL+ gyzL* GamyxzL+ gzzL* GamzxzL;
      double ass_Gamzyy =  gxzL* GamxyyL+ gyzL* GamyyyL+ gzzL* GamzyyL;
      double ass_Gamzyz =  gxzL* GamxyzL+ gyzL* GamyyzL+ gzzL* GamzyzL;
      double ass_Gamzzz =  gxzL* GamxzzL+ gyzL* GamyzzL+ gzzL* GamzzzL;

	
      RxxL= - 0.5 * RxxL+  
	gxxL* Gamxx + gxyL* Gamyx+ gxzL* Gamzx +  
	GammaxL* ass_Gamxxx + GammayL* ass_Gamxxy + GammazL* ass_Gamxxz +  
	gupxxL * ( 
		  2.0 * (  
			 GamxxxL* ass_Gamxxx + GamyxxL* ass_Gamxxy + GamzxxL* ass_Gamxxz ) +  
		  GamxxxL* ass_Gamxxx + GamyxxL* ass_Gamyxx + GamzxxL* ass_Gamzxx ) +  
	gupxyL* (  
		 2.0 * (  
			GamxxxL* ass_Gamxxy + GamyxxL* ass_Gamxyy + GamzxxL* ass_Gamxyz ) +  
		 GamxxyL* ass_Gamxxx + GamyxyL* ass_Gamyxx + GamzxyL* ass_Gamzxx +  
		 2.0 * (  
			GamxxyL* ass_Gamxxx + GamyxyL* ass_Gamxxy + GamzxyL* ass_Gamxxz ) +  
		 GamxxxL* ass_Gamxxy + GamyxxL* ass_Gamyxy + GamzxxL* ass_Gamzxy ) +  
	gupxzL* (  
		 2.0 * (  
			GamxxxL* ass_Gamxxz + GamyxxL* ass_Gamxyz + GamzxxL* ass_Gamxzz ) +  
		 GamxxzL* ass_Gamxxx + GamyxzL* ass_Gamyxx + GamzxzL* ass_Gamzxx +  
		 2.0 * (  
			GamxxzL* ass_Gamxxx + GamyxzL* ass_Gamxxy + GamzxzL* ass_Gamxxz ) +  
		 GamxxxL* ass_Gamxxz + GamyxxL* ass_Gamyxz + GamzxxL* ass_Gamzxz ) +  
	gupyyL* (  
		 2.0 * (  
			GamxxyL* ass_Gamxxy + GamyxyL* ass_Gamxyy + GamzxyL* ass_Gamxyz ) +  
		 GamxxyL* ass_Gamxxy + GamyxyL* ass_Gamyxy + GamzxyL* ass_Gamzxy ) +  
	gupyzL* (  
		 2.0 * (  
			GamxxyL* ass_Gamxxz + GamyxyL* ass_Gamxyz + GamzxyL* ass_Gamxzz ) +  
		 GamxxzL* ass_Gamxxy + GamyxzL* ass_Gamyxy + GamzxzL* ass_Gamzxy +  
		 2.0 * (  
			GamxxzL* ass_Gamxxy + GamyxzL* ass_Gamxyy + GamzxzL* ass_Gamxyz ) +  
		 GamxxyL* ass_Gamxxz + GamyxyL* ass_Gamyxz + GamzxyL* ass_Gamzxz ) +  
	gupzzL* (  
		 2.0 * (  
			GamxxzL* ass_Gamxxz + GamyxzL* ass_Gamxyz + GamzxzL* ass_Gamxzz ) +  
		 GamxxzL* ass_Gamxxz + GamyxzL* ass_Gamyxz + GamzxzL* ass_Gamzxz ); 

	  
      RxyL= 0.5 * ( - RxyL+  
		    gxxL* Gamxy+ gxyL* Gamyy+ gxzL* Gamzy+  
		    gxyL* Gamxx+ gyyL* Gamyx+ gyzL* Gamzx+  
		    GammaxL* ass_Gamxxy + GammayL* ass_Gamxyy + GammazL* ass_Gamxyz +  
		    GammaxL* ass_Gamyxx + GammayL* ass_Gamyxy + GammazL* ass_Gamyxz ) +  
	gupxxL* (  
		 GamxxxL* ass_Gamyxx + GamyxxL* ass_Gamyxy + GamzxxL* ass_Gamyxz +  
		 GamxxyL* ass_Gamxxx + GamyxyL* ass_Gamxxy + GamzxyL* ass_Gamxxz +  
		 GamxxxL* ass_Gamxxy + GamyxxL* ass_Gamyxy + GamzxxL* ass_Gamzxy ) +  
	gupxyL* (  
		 GamxxxL* ass_Gamyxy + GamyxxL* ass_Gamyyy + GamzxxL* ass_Gamyyz +  
		 GamxxyL* ass_Gamxxy + GamyxyL* ass_Gamxyy + GamzxyL* ass_Gamxyz +  
		 GamxxyL* ass_Gamxxy + GamyxyL* ass_Gamyxy + GamzxyL* ass_Gamzxy +  
		 GamxxyL* ass_Gamyxx + GamyxyL* ass_Gamyxy + GamzxyL* ass_Gamyxz +  
		 GamxyyL* ass_Gamxxx + GamyyyL* ass_Gamxxy + GamzyyL* ass_Gamxxz +  
		 GamxxxL* ass_Gamxyy + GamyxxL* ass_Gamyyy + GamzxxL* ass_Gamzyy ) +  
	gupxzL* (  
		 GamxxxL* ass_Gamyxz + GamyxxL* ass_Gamyyz + GamzxxL* ass_Gamyzz +  
		 GamxxyL* ass_Gamxxz + GamyxyL* ass_Gamxyz + GamzxyL* ass_Gamxzz +  
		 GamxxzL* ass_Gamxxy + GamyxzL* ass_Gamyxy + GamzxzL* ass_Gamzxy +  
		 GamxxzL* ass_Gamyxx + GamyxzL* ass_Gamyxy + GamzxzL* ass_Gamyxz +  
		 GamxyzL* ass_Gamxxx + GamyyzL* ass_Gamxxy + GamzyzL* ass_Gamxxz +  
		 GamxxxL* ass_Gamxyz + GamyxxL* ass_Gamyyz + GamzxxL* ass_Gamzyz ) +  
	gupyyL* (  
		 GamxxyL* ass_Gamyxy + GamyxyL* ass_Gamyyy + GamzxyL* ass_Gamyyz +  
		 GamxyyL* ass_Gamxxy + GamyyyL* ass_Gamxyy + GamzyyL* ass_Gamxyz +  
		 GamxxyL* ass_Gamxyy + GamyxyL* ass_Gamyyy + GamzxyL* ass_Gamzyy ) +  
	gupyzL* (  
		 GamxxyL* ass_Gamyxz + GamyxyL* ass_Gamyyz + GamzxyL* ass_Gamyzz +  
		 GamxyyL* ass_Gamxxz + GamyyyL* ass_Gamxyz + GamzyyL* ass_Gamxzz +  
		 GamxxzL* ass_Gamxyy + GamyxzL* ass_Gamyyy + GamzxzL* ass_Gamzyy +  
		 GamxxzL* ass_Gamyxy + GamyxzL* ass_Gamyyy + GamzxzL* ass_Gamyyz +  
		 GamxyzL* ass_Gamxxy + GamyyzL* ass_Gamxyy + GamzyzL* ass_Gamxyz +  
		 GamxxyL* ass_Gamxyz + GamyxyL* ass_Gamyyz + GamzxyL* ass_Gamzyz ) +  
	gupzzL* (  
		 GamxxzL* ass_Gamyxz + GamyxzL* ass_Gamyyz + GamzxzL* ass_Gamyzz +  
		 GamxyzL* ass_Gamxxz + GamyyzL* ass_Gamxyz + GamzyzL* ass_Gamxzz +  
		 GamxxzL* ass_Gamxyz + GamyxzL* ass_Gamyyz + GamzxzL* ass_Gamzyz ); 
	
	
	
      RxzL= 0.5 * ( - RxzL+  
		    gxxL* Gamxz+ gxyL* Gamyz+ gxzL* Gamzz+  
		    gxzL* Gamxx+ gyzL* Gamyx+ gzzL* Gamzx+  
		    GammaxL* ass_Gamxxz + GammayL* ass_Gamxyz + GammazL* ass_Gamxzz +  
		    GammaxL* ass_Gamzxx + GammayL* ass_Gamzxy + GammazL* ass_Gamzxz ) +  
	gupxxL* (  
		 GamxxxL* ass_Gamzxx + GamyxxL* ass_Gamzxy + GamzxxL* ass_Gamzxz +  
		 GamxxzL* ass_Gamxxx + GamyxzL* ass_Gamxxy + GamzxzL* ass_Gamxxz +  
		 GamxxxL* ass_Gamxxz + GamyxxL* ass_Gamyxz + GamzxxL* ass_Gamzxz ) +  
	gupxyL* (  
		 GamxxxL* ass_Gamzxy + GamyxxL* ass_Gamzyy + GamzxxL* ass_Gamzyz +  
		 GamxxzL* ass_Gamxxy + GamyxzL* ass_Gamxyy + GamzxzL* ass_Gamxyz +  
		 GamxxyL* ass_Gamxxz + GamyxyL* ass_Gamyxz + GamzxyL* ass_Gamzxz +  
		 GamxxyL* ass_Gamzxx + GamyxyL* ass_Gamzxy + GamzxyL* ass_Gamzxz +  
		 GamxyzL* ass_Gamxxx + GamyyzL* ass_Gamxxy + GamzyzL* ass_Gamxxz +  
		 GamxxxL* ass_Gamxyz + GamyxxL* ass_Gamyyz + GamzxxL* ass_Gamzyz ) +  
	gupxzL* (  
		 GamxxxL* ass_Gamzxz + GamyxxL* ass_Gamzyz + GamzxxL* ass_Gamzzz +  
		 GamxxzL* ass_Gamxxz + GamyxzL* ass_Gamxyz + GamzxzL* ass_Gamxzz +  
		 GamxxzL* ass_Gamxxz + GamyxzL* ass_Gamyxz + GamzxzL* ass_Gamzxz +  
		 GamxxzL* ass_Gamzxx + GamyxzL* ass_Gamzxy + GamzxzL* ass_Gamzxz +  
		 GamxzzL* ass_Gamxxx + GamyzzL* ass_Gamxxy + GamzzzL* ass_Gamxxz +  
		 GamxxxL* ass_Gamxzz + GamyxxL* ass_Gamyzz + GamzxxL* ass_Gamzzz ) +  
	gupyyL* (  
		 GamxxyL* ass_Gamzxy + GamyxyL* ass_Gamzyy + GamzxyL* ass_Gamzyz +  
		 GamxyzL* ass_Gamxxy + GamyyzL* ass_Gamxyy + GamzyzL* ass_Gamxyz +  
		 GamxxyL* ass_Gamxyz + GamyxyL* ass_Gamyyz + GamzxyL* ass_Gamzyz ) +  
	gupyzL* (  
		 GamxxyL* ass_Gamzxz + GamyxyL* ass_Gamzyz + GamzxyL* ass_Gamzzz +  
		 GamxyzL* ass_Gamxxz + GamyyzL* ass_Gamxyz + GamzyzL* ass_Gamxzz +  
		 GamxxzL* ass_Gamxyz + GamyxzL* ass_Gamyyz + GamzxzL* ass_Gamzyz +  
		 GamxxzL* ass_Gamzxy + GamyxzL* ass_Gamzyy + GamzxzL* ass_Gamzyz +  
		 GamxzzL* ass_Gamxxy + GamyzzL* ass_Gamxyy + GamzzzL* ass_Gamxyz +  
		 GamxxyL* ass_Gamxzz + GamyxyL* ass_Gamyzz + GamzxyL* ass_Gamzzz ) +  
	gupzzL* (  
		 GamxxzL* ass_Gamzxz + GamyxzL* ass_Gamzyz + GamzxzL* ass_Gamzzz +  
		 GamxzzL* ass_Gamxxz + GamyzzL* ass_Gamxyz + GamzzzL* ass_Gamxzz +  
		 GamxxzL* ass_Gamxzz + GamyxzL* ass_Gamyzz + GamzxzL* ass_Gamzzz );
	
	
      RyyL= - 0.5 * RyyL+ 
	gxyL* Gamxy+ gyyL* Gamyy+ gyzL* Gamzy+ 
	GammaxL* ass_Gamyxy + GammayL* ass_Gamyyy + GammazL* ass_Gamyyz + 
	gupxxL* (  
		 2.0 * ( 
			GamxxyL* ass_Gamyxx + GamyxyL* ass_Gamyxy + GamzxyL* ass_Gamyxz ) + 
		 GamxxyL* ass_Gamxxy + GamyxyL* ass_Gamyxy + GamzxyL* ass_Gamzxy ) + 
	gupxyL* (  
		 2.0 * ( 
			GamxxyL* ass_Gamyxy + GamyxyL* ass_Gamyyy + GamzxyL* ass_Gamyyz ) + 
		 GamxyyL* ass_Gamxxy + GamyyyL* ass_Gamyxy + GamzyyL* ass_Gamzxy +  
		 2.0 * ( 
			GamxyyL* ass_Gamyxx + GamyyyL* ass_Gamyxy + GamzyyL* ass_Gamyxz ) + 
		 GamxxyL* ass_Gamxyy + GamyxyL* ass_Gamyyy + GamzxyL* ass_Gamzyy ) + 
	gupxzL* (  
		 2.0 * ( 
			GamxxyL* ass_Gamyxz + GamyxyL* ass_Gamyyz + GamzxyL* ass_Gamyzz ) + 
		 GamxyzL* ass_Gamxxy + GamyyzL* ass_Gamyxy + GamzyzL* ass_Gamzxy + 
		 2.0 * ( 
			GamxyzL* ass_Gamyxx + GamyyzL* ass_Gamyxy + GamzyzL* ass_Gamyxz ) + 
		 GamxxyL* ass_Gamxyz + GamyxyL* ass_Gamyyz + GamzxyL* ass_Gamzyz ) + 
	gupyyL* (  
		 2.0 * ( 
			GamxyyL* ass_Gamyxy + GamyyyL* ass_Gamyyy + GamzyyL* ass_Gamyyz ) + 
		 GamxyyL* ass_Gamxyy + GamyyyL* ass_Gamyyy + GamzyyL* ass_Gamzyy ) + 
	gupyzL* (  
		 2.0 * ( 
			GamxyyL* ass_Gamyxz + GamyyyL* ass_Gamyyz + GamzyyL* ass_Gamyzz ) + 
		 GamxyzL* ass_Gamxyy + GamyyzL* ass_Gamyyy + GamzyzL* ass_Gamzyy + 
		 2.0 * ( 
			GamxyzL* ass_Gamyxy + GamyyzL* ass_Gamyyy + GamzyzL* ass_Gamyyz ) + 
		 GamxyyL* ass_Gamxyz + GamyyyL* ass_Gamyyz + GamzyyL* ass_Gamzyz ) + 
	gupzzL* (  
		 2.0 * ( 
			GamxyzL* ass_Gamyxz + GamyyzL* ass_Gamyyz + GamzyzL* ass_Gamyzz ) + 
		 GamxyzL* ass_Gamxyz + GamyyzL* ass_Gamyyz + GamzyzL* ass_Gamzyz );
	  
	  
      RyzL= 0.5 * ( - RyzL+  
		    gxyL* Gamxz+ gyyL* Gamyz+ gyzL* Gamzz+  
		    gxzL* Gamxy+ gyzL* Gamyy+ gzzL* Gamzy+  
		    GammaxL* ass_Gamyxz + GammayL* ass_Gamyyz + GammazL* ass_Gamyzz +  
		    GammaxL* ass_Gamzxy + GammayL* ass_Gamzyy + GammazL* ass_Gamzyz ) +  
	gupxxL* (  
		 GamxxyL* ass_Gamzxx + GamyxyL* ass_Gamzxy + GamzxyL* ass_Gamzxz +  
		 GamxxzL* ass_Gamyxx + GamyxzL* ass_Gamyxy + GamzxzL* ass_Gamyxz +  
		 GamxxyL* ass_Gamxxz + GamyxyL* ass_Gamyxz + GamzxyL* ass_Gamzxz ) +  
	gupxyL* (  
		 GamxxyL* ass_Gamzxy + GamyxyL* ass_Gamzyy + GamzxyL* ass_Gamzyz +  
		 GamxxzL* ass_Gamyxy + GamyxzL* ass_Gamyyy + GamzxzL* ass_Gamyyz +  
		 GamxyyL* ass_Gamxxz + GamyyyL* ass_Gamyxz + GamzyyL* ass_Gamzxz +  
		 GamxyyL* ass_Gamzxx + GamyyyL* ass_Gamzxy + GamzyyL* ass_Gamzxz +  
		 GamxyzL* ass_Gamyxx + GamyyzL* ass_Gamyxy + GamzyzL* ass_Gamyxz +  
		 GamxxyL* ass_Gamxyz + GamyxyL* ass_Gamyyz + GamzxyL* ass_Gamzyz ) +  
	gupxzL* (  
		 GamxxyL* ass_Gamzxz + GamyxyL* ass_Gamzyz + GamzxyL* ass_Gamzzz +  
		 GamxxzL* ass_Gamyxz + GamyxzL* ass_Gamyyz + GamzxzL* ass_Gamyzz +  
		 GamxyzL* ass_Gamxxz + GamyyzL* ass_Gamyxz + GamzyzL* ass_Gamzxz +  
		 GamxyzL* ass_Gamzxx + GamyyzL* ass_Gamzxy + GamzyzL* ass_Gamzxz +  
		 GamxzzL* ass_Gamyxx + GamyzzL* ass_Gamyxy + GamzzzL* ass_Gamyxz +  
		 GamxxyL* ass_Gamxzz + GamyxyL* ass_Gamyzz + GamzxyL* ass_Gamzzz ) +  
	gupyyL* (  
		 GamxyyL* ass_Gamzxy + GamyyyL* ass_Gamzyy + GamzyyL* ass_Gamzyz +  
		 GamxyzL* ass_Gamyxy + GamyyzL* ass_Gamyyy + GamzyzL* ass_Gamyyz +  
		 GamxyyL* ass_Gamxyz + GamyyyL* ass_Gamyyz + GamzyyL* ass_Gamzyz ) +  
	gupyzL* (  
		 GamxyyL* ass_Gamzxz + GamyyyL* ass_Gamzyz + GamzyyL* ass_Gamzzz +  
		 GamxyzL* ass_Gamyxz + GamyyzL* ass_Gamyyz + GamzyzL* ass_Gamyzz +  
		 GamxyzL* ass_Gamxyz + GamyyzL* ass_Gamyyz + GamzyzL* ass_Gamzyz +  
		 GamxyzL* ass_Gamzxy + GamyyzL* ass_Gamzyy + GamzyzL* ass_Gamzyz +  
		 GamxzzL* ass_Gamyxy + GamyzzL* ass_Gamyyy + GamzzzL* ass_Gamyyz +  
		 GamxyyL* ass_Gamxzz + GamyyyL* ass_Gamyzz + GamzyyL* ass_Gamzzz ) +  
	gupzzL* (  
		 GamxyzL* ass_Gamzxz + GamyyzL* ass_Gamzyz + GamzyzL* ass_Gamzzz +  
		 GamxzzL* ass_Gamyxz + GamyzzL* ass_Gamyyz + GamzzzL* ass_Gamyzz +  
		 GamxyzL* ass_Gamxzz + GamyyzL* ass_Gamyzz + GamzyzL* ass_Gamzzz );
	
	
	
      RzzL= - 0.5 * RzzL+
	gxzL* Gamxz+ gyzL* Gamyz+ gzzL* Gamzz + 
	GammaxL* ass_Gamzxz + GammayL* ass_Gamzyz + GammazL* ass_Gamzzz +
	gupxxL* (
		 2.0 * (
			GamxxzL* ass_Gamzxx + GamyxzL* ass_Gamzxy + GamzxzL* ass_Gamzxz ) +
		 GamxxzL* ass_Gamxxz + GamyxzL* ass_Gamyxz + GamzxzL* ass_Gamzxz ) +
	gupxyL* (
		 2.0 * (
			GamxxzL* ass_Gamzxy + GamyxzL* ass_Gamzyy + GamzxzL* ass_Gamzyz ) +
		 GamxyzL* ass_Gamxxz + GamyyzL* ass_Gamyxz + GamzyzL* ass_Gamzxz +
		 2.0 * (
			GamxyzL* ass_Gamzxx + GamyyzL* ass_Gamzxy + GamzyzL* ass_Gamzxz ) +
		 GamxxzL* ass_Gamxyz + GamyxzL* ass_Gamyyz + GamzxzL* ass_Gamzyz ) +
	gupxzL* (
		 2.0 * (
			GamxxzL* ass_Gamzxz + GamyxzL* ass_Gamzyz + GamzxzL* ass_Gamzzz ) +
		 GamxzzL* ass_Gamxxz + GamyzzL* ass_Gamyxz + GamzzzL* ass_Gamzxz +
		 2.0 * (
			GamxzzL* ass_Gamzxx + GamyzzL* ass_Gamzxy + GamzzzL* ass_Gamzxz ) +
		 GamxxzL* ass_Gamxzz + GamyxzL* ass_Gamyzz + GamzxzL* ass_Gamzzz ) +
	gupyyL* (
		 2.0 * (
			GamxyzL* ass_Gamzxy + GamyyzL* ass_Gamzyy + GamzyzL* ass_Gamzyz ) +
		 GamxyzL* ass_Gamxyz + GamyyzL* ass_Gamyyz + GamzyzL* ass_Gamzyz ) +
	gupyzL* (
		 2.0 * (
			GamxyzL* ass_Gamzxz + GamyyzL* ass_Gamzyz + GamzyzL* ass_Gamzzz ) +
		 GamxzzL* ass_Gamxyz + GamyzzL* ass_Gamyyz + GamzzzL* ass_Gamzyz +
		 2.0 * (
			GamxzzL* ass_Gamzxy + GamyzzL* ass_Gamzyy + GamzzzL* ass_Gamzyz ) +
		 GamxyzL* ass_Gamxzz + GamyyzL* ass_Gamyzz + GamzyzL* ass_Gamzzz ) +
	gupzzL* (
		 2.0 * (
			GamxzzL* ass_Gamzxz + GamyzzL* ass_Gamzyz + GamzzzL* ass_Gamzzz ) +
		 GamxzzL* ass_Gamxzz + GamyzzL* ass_Gamyzz + GamzzzL* ass_Gamzzz );
	
      trRtilde[index] = gupxxL * RxxL + gupyyL * RyyL + gupzzL * RzzL +
	2.0 * ( gupxyL * RxyL + gupxzL * RxzL + gupyzL * RyzL );

      gxxx[index] = gxxxL;
      gxyx[index] = gxyxL;
      gxzx[index] = gxzxL;
      gyyx[index] = gyyxL;
      gyzx[index] = gyzxL;
      gzzx[index] = gzzxL;

      gxxy[index] = gxxyL;
      gxyy[index] = gxyyL;
      gxzy[index] = gxzyL;
      gyyy[index] = gyyyL;
      gyzy[index] = gyzyL;
      gzzy[index] = gzzyL;

      gxxz[index] = gxxzL;
      gxyz[index] = gxyzL;
      gxzz[index] = gxzzL;
      gyyz[index] = gyyzL;
      gyzz[index] = gyzzL;
      gzzz[index] = gzzzL;


      Gamxxx[index] =  GamxxxL;
      Gamyxx[index] =  GamyxxL;
      Gamzxx[index] =  GamzxxL;
	
      Gamxyy[index] =  GamxyyL;
      Gamyyy[index] =  GamyyyL;
      Gamzyy[index] =  GamzyyL;
	
      Gamxzz[index] =  GamxzzL; 
      Gamyzz[index] =  GamyzzL; 
      Gamzzz[index] =  GamzzzL; 
	
      Gamxxy[index] =  GamxxyL; 
      Gamyxy[index] =  GamyxyL; 
      Gamzxy[index] =  GamzxyL; 
	
      Gamxxz[index] =  GamxxzL; 
      Gamyxz[index] =  GamyxzL; 
      Gamzxz[index] =  GamzxzL; 
	
      Gamxyz[index] =  GamxyzL; 
      Gamyyz[index] =  GamyyzL;
      Gamzyz[index] =  GamzyzL;

      Rxx[index] = RxxL;
      Rxy[index] = RxyL;
      Rxz[index] = RxzL;
      Ryy[index] = RyyL;
      Ryz[index] = RyzL;
      Rzz[index] = RzzL;
    }
    }
  }
}


extern "C" void CCTK_FCALL CCTK_FNAME(BSSN_ricci_AHinterior)
  (const cGH **cctkGH,double *dT, double *dx, double *dy, double *dz,
   int *nghostzones,int *cctk_lsh,double *r,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,
   double *trRtilde,
   double *Gammax,double *Gammay,double *Gammaz,
   double *gxxx, double *gxxy, double *gxxz,
   double *gxyx, double *gxyy, double *gxyz,
   double *gxzx, double *gxzy, double *gxzz,
   double *gyyx, double *gyyy, double *gyyz,
   double *gyzx, double *gyzy, double *gyzz,
   double *gzzx, double *gzzy, double *gzzz,
   double *Gammaxxx, double *Gammaxxy, double *Gammaxxz, double *Gammaxyy, double *Gammaxyz, double *Gammaxzz,
   double *Gammayxx, double *Gammayxy, double *Gammayxz, double *Gammayyy, double *Gammayyz, double *Gammayzz,
   double *Gammazxx, double *Gammazxy, double *Gammazxz, double *Gammazyy, double *Gammazyz, double *Gammazzz,
   double *Sx, double *Sy, double *Sz,
   double *Aupxx,double *Aupxy,double *Aupxz,double *Aupyy,double *Aupyz,double *Aupzz,
   double *phi, double *trK, double *MResx,double *MResy,double *MResz,double *MNorm, 
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *psi, double *rho, double *PsiRes, double *PsiNorm,int *compute_constraint_flag)
{
  BSSN_ricci_AHinterior(*cctkGH,  *dT,  *dx,  *dy,  *dz,
					nghostzones, cctk_lsh,r,
			     gxx, gxy, gxz, gyy, gyz, gzz,
			     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
			     Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
			     trRtilde,
			     Gammax, Gammay, Gammaz,
			     gxxx, gxxy, gxxz,
			     gxyx, gxyy, gxyz,
			     gxzx, gxzy, gxzz,
			     gyyx, gyyy, gyyz,
			     gyzx, gyzy, gyzz,
			     gzzx, gzzy, gzzz,
			     Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz,
			     Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz,
			     Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzz,
			     Sx,Sy,Sz,
			     Aupxx,Aupxy,Aupxz,Aupyy,Aupyz,Aupzz,
			     phi, trK, MResx,MResy,MResz,MNorm, 
			     Axx,Axy,Axz,Ayy,Ayz,Azz,
			     psi, rho, PsiRes, PsiNorm,*compute_constraint_flag);
}

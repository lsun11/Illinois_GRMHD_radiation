//-----------------------------------------------------------------------------
// Convert various quantities from physical coords to fisheye
//-----------------------------------------------------------------------------

#include "cctk.h"
#include "math.h"

#define F1o3 0.3333333333333333333333333333

extern "C" void CCTK_FCALL trans_phys_fish_tensor_flat_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz);

extern "C" void CCTK_FCALL trans_phys_fish_tensor_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz);

extern "C" void CCTK_FCALL trans_phys_fish_tensor_inv_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz);


extern "C" void CCTK_FCALL trans_phys_fish_gamt_flat_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,double *d2RG,
   double *Gammax,double *Gammay,double *Gammaz);

extern "C" void CCTK_FCALL trans_phys_fish_phi_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *phi);
  
extern "C" void trans_phys_fish_tensor_flat(const cGH *cctkGH,int *cctk_lsh,int Symmetry,
						      double *xG, double *yG, double *zG, 
						      double *RpG,double *dRG,
						      double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) {

  /* Set up variables used in the grid loop for the physical grid points */
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    double Rp = RpG[index];
    double dR = dRG[index];
    double x = xG[index];
    double y = yG[index];
    double z = zG[index];
    double r = sqrt(x*x + y*y + z*z);

    double fac = dR - Rp/r;
    double Jxx = Rp/r + x/r * x/r * fac;
    double Jxy = x/r * y/r * fac;
    double Jxz = x/r * z/r * fac;
    double Jyy = Rp/r + y/r * y/r * fac;
    double Jyz = y/r * z/r * fac;
    double Jzz = Rp/r + z/r * z/r * fac;

    fac = pow((Rp/r),(-4.0/3.0)) * pow(dR,(-2.0/3.0));

    gxx[index] = (Jxx*Jxx + Jxy*Jxy + Jxz*Jxz)*fac;
    gxy[index] = (Jxx*Jxy + Jxy*Jyy + Jxz*Jyz)*fac;
    gxz[index] = (Jxx*Jxz + Jxy*Jyz + Jxz*Jzz)*fac;
    gyy[index] = (Jxy*Jxy + Jyy*Jyy + Jyz*Jyz)*fac;
    gyz[index] = (Jxy*Jxz + Jyy*Jyz + Jyz*Jzz)*fac;
    gzz[index] = (Jxz*Jxz + Jyz*Jyz + Jzz*Jzz)*fac;

  }
}



extern "C" void trans_phys_fish_tensor(const cGH *cctkGH,int *cctk_lsh,int Symmetry,
						 double *xG, double *yG, double *zG, 
						 double *RpG,double *dRG,
						 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) {
  

  /* Set up variables used in the grid loop for the physical grid points */
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    double Rp = RpG[index];
    double dR = dRG[index];
    double x = xG[index];
    double y = yG[index];
    double z = zG[index];
    double r = sqrt(x*x + y*y + z*z);

    double fac = dR - Rp/r;
    double Jxx = Rp/r + x/r * x/r * fac;
    double Jxy = x/r * y/r * fac;
    double Jxz = x/r * z/r * fac;
    double Jyy = Rp/r + y/r * y/r * fac;
    double Jyz = y/r * z/r * fac;
    double Jzz = Rp/r + z/r * z/r * fac;

    fac = pow((Rp/r),(-4.0/3.0)) * pow(dR,(-2.0/3.0));

    double gxxi=gxx[index];
    double gxyi=gxy[index];
    double gxzi=gxz[index];
    double gyyi=gyy[index];
    double gyzi=gyz[index];
    double gzzi=gzz[index];

    gxx[index] = (Jxx*Jxx*gxxi + 2.0*Jxx*Jxy*gxyi + 2.0*Jxx*Jxz*gxzi + 
		  Jxy*Jxy*gyyi + 2.0*Jxy*Jxz*gyzi + Jxz*Jxz*gzzi)*fac;
    gxy[index] = (Jxx*Jxy*gxxi + (Jxx*Jyy+Jxy*Jxy)*gxyi + 
		  (Jxx*Jyz+Jxz*Jxy)*gxzi + Jxy*Jyy*gyyi + 
		  (Jxy*Jyz+Jxz*Jyy)*gyzi + Jxz*Jyz*gzzi)*fac;
    gxz[index] = (Jxx*Jxz*gxxi + (Jxx*Jyz+Jxy*Jxz)*gxyi + 
		  (Jxx*Jzz+Jxz*Jxz)*gxzi + Jxy*Jyz*gyyi + 
		  (Jxy*Jzz+Jxz*Jyz)*gyzi + Jxz*Jzz*gzzi)*fac;
    gyy[index] = (Jxy*Jxy*gxxi + 2.0*Jxy*Jyy*gxyi + 
		  2.0*Jxy*Jyz*gxzi + Jyy*Jyy*gyyi + 
		  2.0*Jyy*Jyz*gyzi + Jyz*Jyz*gzzi)*fac;
    gyz[index] = (Jxy*Jxz*gxxi + (Jxy*Jyz+Jyy*Jxz)*gxyi + 
		  (Jxy*Jzz+Jyz*Jxz)*gxzi + Jyy*Jyz*gyyi + 
		  (Jyy*Jzz+Jyz*Jyz)*gyzi + Jyz*Jzz*gzzi)*fac;
    gzz[index] = (Jxz*Jxz*gxxi + 2.0*Jxz*Jyz*gxyi + 
		  2.0*Jxz*Jzz*gxzi + Jyz*Jyz*gyyi + 
		  2.0*Jyz*Jzz*gyzi + Jzz*Jzz*gzzi)*fac;
    
  }
}



extern "C" void trans_phys_fish_tensor_inv(const cGH *cctkGH,int *cctk_lsh,int Symmetry,
						     double *xG, double *yG, double *zG, 
						     double *RpG,double *dRG,
						     double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) {

  /* Set up variables used in the grid loop for the physical grid points */
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    double Rp = RpG[index];
    double dR = dRG[index];
    double x = xG[index];
    double y = yG[index];
    double z = zG[index];
    double r = sqrt(x*x + y*y + z*z);

    double fac = 1.0/dR - r/Rp;
    double Jxx = r/Rp + x/r * x/r * fac;
    double Jxy = x/r * y/r * fac;
    double Jxz = x/r * z/r * fac;
    double Jyy = r/Rp + y/r * y/r * fac;
    double Jyz = y/r * z/r * fac;
    double Jzz = r/Rp + z/r * z/r * fac;
    fac = pow((r/Rp),(-4.0/3.0)) * pow(dR,(2.0/3.0));

    double gxxi=gxx[index];
    double gyyi=gyy[index];
    double gzzi=gzz[index];
    double gxyi=gxy[index];
    double gxzi=gxz[index];
    double gyzi=gyz[index];
    
    gxx[index] = (Jxx*Jxx*gxxi + 2.0*Jxx*Jxy*gxyi + 2.0*Jxx*Jxz*gxzi + 
	   Jxy*Jxy*gyyi + 2.0*Jxy*Jxz*gyzi + Jxz*Jxz*gzzi)*fac;
    gxy[index] = (Jxx*Jxy*gxxi + (Jxx*Jyy+Jxy*Jxy)*gxyi + 
	   (Jxx*Jyz+Jxz*Jxy)*gxzi + Jxy*Jyy*gyyi + 
	   (Jxy*Jyz+Jxz*Jyy)*gyzi + Jxz*Jyz*gzzi)*fac;
    gxz[index] = (Jxx*Jxz*gxxi + (Jxx*Jyz+Jxy*Jxz)*gxyi + 
	   (Jxx*Jzz+Jxz*Jxz)*gxzi + Jxy*Jyz*gyyi + 
	   (Jxy*Jzz+Jxz*Jyz)*gyzi + Jxz*Jzz*gzzi)*fac;
    gyy[index] = (Jxy*Jxy*gxxi + 2.0*Jxy*Jyy*gxyi + 
	   2.0*Jxy*Jyz*gxzi + Jyy*Jyy*gyyi + 
	   2.0*Jyy*Jyz*gyzi + Jyz*Jyz*gzzi)*fac;
    gyz[index] = (Jxy*Jxz*gxxi + (Jxy*Jyz+Jyy*Jxz)*gxyi + 
	   (Jxy*Jzz+Jyz*Jxz)*gxzi + Jyy*Jyz*gyyi + 
	   (Jyy*Jzz+Jyz*Jyz)*gyzi + Jyz*Jzz*gzzi)*fac;
    gzz[index] = (Jxz*Jxz*gxxi + 2.0*Jxz*Jyz*gxyi + 
	   2.0*Jxz*Jzz*gxzi + Jyz*Jyz*gyyi + 
	   2.0*Jyz*Jzz*gyzi + Jzz*Jzz*gzzi)*fac;
    
  }
}



extern "C" void trans_phys_fish_phi(const cGH *cctkGH,int *cctk_lsh,int Symmetry,
					      double *xG, double *yG, double *zG, 
					      double *RpG,double *dRG,
					      double *phi) {
  /* Set up variables used in the grid loop for the physical grid points */
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    double x = xG[index];
    double y = yG[index];
    double z = zG[index];
    double r = sqrt(x*x + y*y + z*z);
    double Rp = RpG[index];
    double dR = dRG[index];
   
    phi[index] += 1.0/3.0*log(Rp/r)+1.0/6.0*log(dR);

  }
}


extern "C" void  trans_phys_fish_gamt_flat(const cGH *cctkGH,int *cctk_lsh,int Symmetry,
						     double *xG, double *yG, double *zG, 
						     double *RpG,double *dRG,double *d2RG,
						     double *Gammax,double *Gammay,double *Gammaz) {
  

  /* Set up variables used in the grid loop for the physical grid points */
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    double x = xG[index];
    double y = yG[index];
    double z = zG[index];
    double r = sqrt(x*x + y*y + z*z);

    double nx=x/r;
    double ny=y/r;
    double nz=z/r;
 
    double RP = RpG[index];
    double RD = dRG[index];
    double RD2 = d2RG[index];

    //   Now the Jacobian terms (we need to transform the downstairs version)
    double juxx=r/RP+nx*nx*(1.0/RD-r/RP);
    double juyy=r/RP+ny*ny*(1.0/RD-r/RP);
    double juzz=r/RP+nz*nz*(1.0/RD-r/RP);
    double juxy=nx*ny*(1.0/RD-r/RP);
    double juxz=nx*nz*(1.0/RD-r/RP);
    double juyz=ny*nz*(1.0/RD-r/RP);
       
    // (phii/phii_p)^4 = (RP/r)^4/3 * (dr_p/dr_fish)^2/3
    double p1=pow((RP/r),(4.0/3.0)) * pow(RD,(2.0/3.0));
    
    // dxp is /partial_{x-bar} (RP/r)
    double dxp=nx/r*(RD-RP/r);
    double dyp=ny/r*(RD-RP/r);
    double dzp=nz/r*(RD-RP/r);

    dxp=2.0/3.0*p1*(2.0*dxp*r/RP+RD2/RD*nx);
    dyp=2.0/3.0*p1*(2.0*dyp*r/RP+RD2/RD*ny);
    dzp=2.0/3.0*p1*(2.0*dzp*r/RP+RD2/RD*nz);

    //
    // and now the Gamma's
    //                   =del_ij       +del_ik  +del_jk  -terms
    double g2=RP/r-RD;
    double g3=1.0/RD-r/RP;

    double RD_squared = RD*RD;
  
    double jxxx=nx*r/(RP*RP)*g2+nx/r*g3+nx/r*g3-2*nx*nx*nx/r*g3 
      -nx*nx*nx*r/(RP*RP)*g2-nx*nx*nx*RD2/RD_squared;
    double jxxy=ny*r/(RP*RP)*g2                  -2*nx*nx*ny/r*g3 
      -nx*nx*ny*r/(RP*RP)*g2-nx*nx*ny*RD2/RD_squared;
    double jxxz=nz*r/(RP*RP)*g2                  -2*nx*nx*nz/r*g3 
      -nx*nx*nz*r/(RP*RP)*g2-nx*nx*nz*RD2/RD_squared;
  
    double jyyx=nx*r/(RP*RP)*g2                  -2*ny*ny*nx/r*g3 
      -ny*ny*nx*r/(RP*RP)*g2-ny*ny*nx*RD2/RD_squared;
    double jyyy=ny*r/(RP*RP)*g2+ny/r*g3+ny/r*g3-2*ny*ny*ny/r*g3 
      -ny*ny*ny*r/(RP*RP)*g2-ny*ny*ny*RD2/RD_squared;
    double jyyz=nz*r/(RP*RP)*g2                  -2*ny*ny*nz/r*g3 
      -ny*ny*nz*r/(RP*RP)*g2-ny*ny*nz*RD2/RD_squared;
  
    double jzzx=nx*r/(RP*RP)*g2                  -2*nz*nz*nx/r*g3 
      -nz*nz*nx*r/(RP*RP)*g2-nz*nz*nx*RD2/RD_squared;
    double jzzy=ny*r/(RP*RP)*g2                  -2*nz*nz*ny/r*g3 
      -nz*nz*ny*r/(RP*RP)*g2-nz*nz*ny*RD2/RD_squared;
    double jzzz=nz*r/(RP*RP)*g2+nz/r*g3+nz/r*g3-2*nz*nz*nz/r*g3 
      -nz*nz*nz*r/(RP*RP)*g2-nz*nz*nz*RD2/RD_squared;
  
    double jxyx=                       ny/r*g3-2*nx*ny*nx/r*g3 
      -nx*ny*nx*r/(RP*RP)*g2-nx*ny*nx*RD2/RD_squared;
    double jxyy=                       nx/r*g3-2*nx*ny*ny/r*g3 
      -nx*ny*ny*r/(RP*RP)*g2-nx*ny*ny*RD2/RD_squared;
    double jxyz=                               -2*nx*ny*nz/r*g3 
      -nx*ny*nz*r/(RP*RP)*g2-nx*ny*nz*RD2/RD_squared;
  
    double jxzx=                       nz/r*g3-2*nx*nz*nx/r*g3 
      -nx*nz*nx*r/(RP*RP)*g2-nx*nz*nx*RD2/RD_squared;
    double jxzy=                               -2*nx*nz*ny/r*g3 
      -nx*nz*ny*r/(RP*RP)*g2-nx*nz*ny*RD2/RD_squared;
    double jxzz=                       nx/r*g3-2*nx*nz*nz/r*g3 
      -nx*nz*nz*r/(RP*RP)*g2-nx*nz*nz*RD2/RD_squared;
  
    double jyzx=                               -2*ny*nz*nx/r*g3 
      -ny*nz*nx*r/(RP*RP)*g2-ny*nz*nx*RD2/RD_squared;
    double jyzy=                       nz/r*g3-2*ny*nz*ny/r*g3 
      -ny*nz*ny*r/(RP*RP)*g2-ny*nz*ny*RD2/RD_squared;
    double jyzz=                       ny/r*g3-2*ny*nz*nz/r*g3 
      -ny*nz*nz*r/(RP*RP)*g2-ny*nz*nz*RD2/RD_squared;
  
    Gammax[index] = -1.0*(dxp*(juxx*juxx+juxy*juxy+juxz*juxz)+ 
			  2.0*p1*(juxx*jxxx+juxy*jxyx+juxz*jxzx)+ 
			  dyp*(juxx*juxy+juxy*juyy+juxz*juyz)+ 
			  p1*(juxx*jxyy+juxy*jyyy+juxz*jyzy+ 
			      jxxy*juxy+jxyy*juyy+jxzy*juyz)+ 
			  dzp*(juxx*juxz+juxy*juyz+juxz*juzz)+ 
			  p1*(juxx*jxzz+juxy*jyzz+juxz*jzzz+ 
			      jxxz*juxz+jxyz*juyz+jxzz*juzz));
    Gammay[index] = -1.0*(dxp*(juxx*juxy+juxy*juyy+juxz*juyz)+ 
			  p1*(juxx*jxyx+juxy*jyyx+juxz*jyzx+ 
			      jxxx*juxy+jxyx*juyy+jxzx*juyz)+ 
			  dyp*(juxy*juxy+juyy*juyy+juyz*juyz)+ 
			  2.0*p1*(juxy*jxyy+juyy*jyyy+juyz*jyzy)+ 
			  dzp*(juxy*juxz+juyy*juyz+juyz*juzz)+ 
			  p1*(juxy*jxzz+juyy*jyzz+juyz*jzzz+ 
			      jxyz*juxz+jyyz*juyz+jyzz*juzz));
    Gammaz[index] = -1.0*(dxp*(juxx*juxz+juxy*juyz+juxz*juzz)+ 
			  p1*(juxx*jxzx+juxy*jyzx+juxz*jzzx+ 
			      jxxx*juxz+jxyx*juyz+jxzx*juzz)+ 
			  dyp*(juxy*juxz+juyy*juyz+juyz*juzz)+ 
			  p1*(juxy*jxzy+juyy*jyzy+juyz*jzzy+ 
			      jxyy*juxz+jyyy*juyz+jyzy*juzz)+ 
			  dzp*(juxz*juxz+juyz*juyz+juzz*juzz)+ 
			  2.0*p1*(juxz*jxzz+juyz*jyzz+juzz*jzzz)) ;

  }
}

extern "C" void CCTK_FCALL trans_phys_fish_tensor_flat_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) 
{
  trans_phys_fish_tensor_flat(*cctkGH,cctk_lsh,*Symmetry,
					xG, yG, zG, 
					RpG,dRG,
					gxx,gxy,gxz,gyy,gyz,gzz); 
}
extern "C" void CCTK_FCALL trans_phys_fish_tensor_inv_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) 
{
  trans_phys_fish_tensor_inv(*cctkGH,cctk_lsh,*Symmetry,
				       xG, yG, zG, 
				       RpG,dRG,
				       gxx,gxy,gxz,gyy,gyz,gzz);
}

extern "C" void CCTK_FCALL trans_phys_fish_tensor_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz)
{
  trans_phys_fish_tensor(*cctkGH,cctk_lsh,*Symmetry,
				   xG, yG, zG, 
				   RpG,dRG,
				   gxx,gxy,gxz,gyy,gyz,gzz);
}


extern "C" void CCTK_FCALL trans_phys_fish_phi_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,
   double *phi)
{
  trans_phys_fish_phi(*cctkGH,cctk_lsh,*Symmetry,
				xG, yG, zG, 
				RpG,dRG,
				phi);
}

extern "C" void CCTK_FCALL trans_phys_fish_gamt_flat_
  (const cGH **cctkGH,int *cctk_lsh,int *Symmetry,
   double *xG, double *yG, double *zG, 
   double *RpG,double *dRG,double *d2RG,
   double *Gammax,double *Gammay,double *Gammaz)
{
  trans_phys_fish_gamt_flat(*cctkGH,cctk_lsh,*Symmetry,
				      xG, yG, zG, 
				      RpG,dRG,d2RG,
				      Gammax,Gammay,Gammaz);
}

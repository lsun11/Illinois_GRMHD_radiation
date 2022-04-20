//--------------------------------------------------------------------
// rhs 
//--------------------------------------------------------------------

#define SQR(x) ((x) * (x))

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

void average_cpp(const cGH *cctkGH,int *ext,double *f,double *favg,
                 int m );

extern "C" void CCTK_FCALL CCTK_FNAME(mhd_A_psi6phi_rhs_cpp)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *psi6phi_rhs, double *Ax_rhs, double *Ay_rhs, double *Az_rhs,
   double *psi6phi, double *Ax, double *Ay, double *Az,
   double *lapm1, double *shiftx, double *shifty, double *shiftz, double *phi, 
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *temp1, double *temp2, double *temp3, double *temp4, double *temp5,
   double *temp6, double *temp7, double *temp8, double *temp9, double *temp10,
   double *temp11, double *dX,double *dY,double *dZ,
   double &x_posn,double &y_posn,double &z_posn,double &time,
   double *r, double *damp_lorenz);

extern "C" void mhd_A_psi6phi_rhs_cpp(const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
				      double *psi6phi_rhs, double *Ax_rhs, double *Ay_rhs, double *Az_rhs,
				      double *psi6phi, double *Ax, double *Ay, double *Az,
                                      double *lapm1, double *shiftx, double *shifty, double *shiftz, double *phi, 
                                      double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
                                      double *temp1, double *temp2, double *temp3, double *temp4, double *temp5, 
                                      double *temp6, double *temp7, double *temp8, double *temp9, double *temp10,
    			              double *temp11,
				      double dX,double dY,double dZ,
				      double &x_posn,double &y_posn,double &z_posn,double &time,
				      double *r, double damp_lorenz) {


  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  double dXm1 = 1.0/dX, dYm1 = 1.0/dY, dZm1 = 1.0/dZ;

  // compute alpha psi^2 to temp1 and psi^6 to temp11
#pragma omp parallel for
  for(int k=0;k<kend;k++) for(int j=0;j<jend;j++) for(int i=0;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double Psi2 = exp(2.0*phi[index]);
	temp1[index] = (1.0+lapm1[index])*Psi2;
	temp11[index] = Psi2*Psi2*Psi2;
      }

  // Compute gf at points i,j+1/2,k+1/2
  // temp2: gupxx
  // temp3: gupxy
  // temp4: gupxz
  // temp5: alpha psi^2
  // temp6: A_y
  // temp7: A_z
  average_cpp(cctkGH,cctk_lsh,gupxx,temp2,23);
  average_cpp(cctkGH,cctk_lsh,gupxy,temp3,23);
  average_cpp(cctkGH,cctk_lsh,gupxz,temp4,23);
  average_cpp(cctkGH,cctk_lsh,temp1,temp5,23);
  average_cpp(cctkGH,cctk_lsh,Ay,temp6,92);
  average_cpp(cctkGH,cctk_lsh,Az,temp7,93);
  
  // compute alpha psi^6 A^x and store at temp8
#pragma omp parallel for
  for(int k=0;k<kend;k++) for(int j=0;j<jend;j++) for(int i=0;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	temp8[index] = (temp2[index]*Ax[index] + temp3[index]*temp6[index] + temp4[index]*temp7[index])*temp5[index];
      }  

  // Compute gf at points i+1/2,j,k+1/2
  // temp2: gupxy
  // temp3: gupyy
  // temp4: gupyz
  // temp5: alpha psi^2
  // temp6: A_x
  // temp7: A_z
  average_cpp(cctkGH,cctk_lsh,gupxy,temp2,13);
  average_cpp(cctkGH,cctk_lsh,gupyy,temp3,13);
  average_cpp(cctkGH,cctk_lsh,gupyz,temp4,13);
  average_cpp(cctkGH,cctk_lsh,temp1,temp5,13);
  average_cpp(cctkGH,cctk_lsh,Ax,temp6,18);
  average_cpp(cctkGH,cctk_lsh,Az,temp7,83);

  // compute alpha psi^6 A^y and store at temp9
#pragma omp parallel for
  for(int k=0;k<kend;k++) for(int j=0;j<jend;j++) for(int i=0;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	temp9[index] = (temp2[index]*temp6[index] + temp3[index]*Ay[index] + temp4[index]*temp7[index])*temp5[index];
      }

  // Compute gf at points i+1/2,j+1/2,k
  // temp2: gupxz
  // temp3: gupyz
  // temp4: gupzz
  // temp5: alpha psi^2
  // temp6: A_x
  // temp7: A_y
  average_cpp(cctkGH,cctk_lsh,gupxz,temp2,12);
  average_cpp(cctkGH,cctk_lsh,gupyz,temp3,12);
  average_cpp(cctkGH,cctk_lsh,gupzz,temp4,12);
  average_cpp(cctkGH,cctk_lsh,temp1,temp5,12);
  average_cpp(cctkGH,cctk_lsh,Ax,temp6,17);
  average_cpp(cctkGH,cctk_lsh,Ay,temp7,27);

  // compute alpha psi^6 A^z and store at temp10
#pragma omp parallel for
  for(int k=0;k<kend;k++) for(int j=0;j<jend;j++) for(int i=0;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	temp10[index] = (temp2[index]*temp6[index] + temp3[index]*temp7[index] + temp4[index]*Az[index])*temp5[index];
      }

  // Compute beta^i at i+1/2,j+1/2,k+1/2 and store in temp2-temp4
  average_cpp(cctkGH,cctk_lsh,shiftx,temp2,123);
  average_cpp(cctkGH,cctk_lsh,shifty,temp3,123);
  average_cpp(cctkGH,cctk_lsh,shiftz,temp4,123);

  //store alpha in temp5
  average_cpp(cctkGH,cctk_lsh,lapm1,temp5,123);

  // Now compute psi6phi_rhs 
#pragma omp parallel for
  for(int k=2;k<kend-2;k++) for(int j=2;j<jend-2;j++) for(int i=2;i<iend-2;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	int index1x = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
	int index1y = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
	int index1z = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
	psi6phi_rhs[index] = dXm1*(temp8[index]-temp8[index1x]) 
	  + dYm1*(temp9[index]-temp9[index1y]) 
	  + dZm1*(temp10[index]-temp10[index1z]);

	//psi6phi_rhs[index] += 1.e-5*sin(3.0*r[index]-3.0*time);
	//psi6phi_rhs[index] += -2.0*1.e-5*sin(2.0*r[index]-2.0*time)*(0.5-erf((r[index]*20-100)/40)*0.5);

	// use upwind scheme to handle partial_j (psi^6 Phi beta^j)
	if (temp2[index] < 0) { 
	  int indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
	  int indexm2 = CCTK_GFINDEX3D(cctkGH,i-2,j,k);
	  psi6phi_rhs[index] += 0.5*dXm1*( temp2[indexm2]*psi6phi[indexm2] 
					   - 4.0*temp2[indexm1]*psi6phi[indexm1] 
					   + 3.0*temp2[index]*psi6phi[index]);
	} else { 
	  int index1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
	  int index2 = CCTK_GFINDEX3D(cctkGH,i+2,j,k);
	  psi6phi_rhs[index] += 0.5*dXm1*( - temp2[index2]*psi6phi[index2]
					   + 4.0*temp2[index1]*psi6phi[index1]
					   - 3.0*temp2[index]*psi6phi[index]);
	}
	if (temp3[index] < 0) {
	  int indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
	  int indexm2 = CCTK_GFINDEX3D(cctkGH,i,j-2,k);
	  psi6phi_rhs[index] += 0.5*dYm1*(  temp3[indexm2]*psi6phi[indexm2]
					    - 4.0*temp3[indexm1]*psi6phi[indexm1]
					    + 3.0*temp3[index]*psi6phi[index]);
	} else {
	  int index1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
	  int index2 = CCTK_GFINDEX3D(cctkGH,i,j+2,k);
	  psi6phi_rhs[index] += 0.5*dYm1*( - temp3[index2]*psi6phi[index2]
					   + 4.0*temp3[index1]*psi6phi[index1]
					   - 3.0*temp3[index]*psi6phi[index]);
	}
	if (temp4[index] < 0) {
	  int indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
	  int indexm2 = CCTK_GFINDEX3D(cctkGH,i,j,k-2);
	  psi6phi_rhs[index] += 0.5*dZm1*(  temp4[indexm2]*psi6phi[indexm2]
					    - 4.0*temp4[indexm1]*psi6phi[indexm1]
					    + 3.0*temp4[index]*psi6phi[index]);
	} else {
	  int index1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
	  int index2 = CCTK_GFINDEX3D(cctkGH,i,j,k+2);
	  psi6phi_rhs[index] += 0.5*dZm1*( - temp4[index2]*psi6phi[index2]
					   + 4.0*temp4[index1]*psi6phi[index1]
					   - 3.0*temp4[index]*psi6phi[index]);
	}
	if (damp_lorenz>0.0){
	  //subtract lambda * alpha psi^6 Phi     
	  psi6phi_rhs[index] += -damp_lorenz*(temp5[index]+1.0)*psi6phi[index]; 
	}
      }

  // Store alpha psi^-6 to temp5
#pragma omp parallel for
  for(int k=0;k<kend;k++) for(int j=0;j<jend;j++) for(int i=0;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	temp5[index] = (1.0+lapm1[index])/temp11[index];
      }

  // Compute gf at i+1/2,j+1/2,k+1/2
  // temp1: alpha psi^-6
  // temp2-temp4: shiftx, shifty, shiftz (has been computed)
  // temp6-temp8: A_x, A_y, A_z
  average_cpp(cctkGH,cctk_lsh,temp5,temp1,123);
  //average_cpp(cctkGH,cctk_lsh,shiftx,temp2,123);
  //average_cpp(cctkGH,cctk_lsh,shifty,temp3,123);
  //average_cpp(cctkGH,cctk_lsh,shiftz,temp4,123);
  average_cpp(cctkGH,cctk_lsh,Ax,temp6,1);
  average_cpp(cctkGH,cctk_lsh,Ay,temp7,2);
  average_cpp(cctkGH,cctk_lsh,Az,temp8,3);

  // Compute alpha Phi - beta^i A_i at points i+1/2,j+1/2,k+1/2 and store at temp9
#pragma omp parallel for
  for(int k=0;k<kend;k++) for(int j=0;j<jend;j++) for(int i=0;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	temp9[index] = temp1[index]*psi6phi[index] - temp2[index]*temp6[index] 
	  - temp3[index]*temp7[index] - temp4[index]*temp8[index];
      }

  // Now add -grad(alpha Phi - beta^i A_i) to A_i_rhs 
#pragma omp parallel for
  for(int k=1;k<kend;k++) for(int j=1;j<jend;j++) for(int i=1;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	int indexm1x = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
	int indexm1y = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
	int indexm1z = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
	Ax_rhs[index] += dXm1*(temp9[indexm1x]-temp9[index]);
	Ay_rhs[index] += dYm1*(temp9[indexm1y]-temp9[index]);
	Az_rhs[index] += dZm1*(temp9[indexm1z]-temp9[index]);
      }
  
}

extern "C" void CCTK_FCALL CCTK_FNAME(mhd_A_psi6phi_rhs_cpp)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *psi6phi_rhs, double *Ax_rhs, double *Ay_rhs, double *Az_rhs,
   double *psi6phi, double *Ax, double *Ay, double *Az,
   double *lapm1, double *shiftx, double *shifty, double *shiftz, double *phi,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *temp1, double *temp2, double *temp3, double *temp4, double *temp5,
   double *temp6, double *temp7, double *temp8, double *temp9, double *temp10,
   double *temp11, double *dx,double *dy,double *dz,
   double &x_posn,double &y_posn,double &z_posn,double &time,
   double *r,double *damp_lorenz)
{
  mhd_A_psi6phi_rhs_cpp(*cctkGH,cctk_lsh, nghostzones, *Symmetry,
	                psi6phi_rhs, Ax_rhs,Ay_rhs,Az_rhs,
			psi6phi,Ax,Ay,Az,lapm1,shiftx,shifty,shiftz,phi,
			gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, 
	    	        temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,
			temp9,temp10,temp11, *dx,*dy,*dz,x_posn,y_posn,z_posn,time,r,*damp_lorenz);
}

// Calculate the average of a gridfunction f; output to gridfunction favg; 
// m = 1   (average over x: i,j,k -> i+1/2,j,k);
// m = 2   (average over y: i,j,k -> i,j+1/2,k);
// m = 3   (average over z: i,j,k -> i,j,k+1/2);
// m = 12  (average over x and y: i,j,k -> i+1/2,j+1/2,k);
// m = 13  (average over x and z: i,j,k -> i+1/2,j,k+1/2);
// m = 23  (average over y and z: i,j,k -> i,j+1/2,k+1/2);
// m = 123 (average over x, y and z: i,j,k -> i+1/2,j+1/2,k+1/2);
// m = 92  (average over x and y: i,j,k -> i-1/2,j+1/2,k);
// m = 93  (average over x and z: i,j,k -> i-1/2,j,k+1/2);
// m = 18  (average over x and y: i,j,k -> i+1/2,j-1/2,k);
// m = 83  (average over y and z: i,j,k -> i,j-1/2,k+1/2);
// m = 17  (average over x and z: i,j,k -> i+1/2,j,k-1/2);
// m = 27  (average over y and z: i,j,k -> i,j+1/2,k-1/2);

//extern "C" void CCTK_FCALL CCTK_FNAME(average_cpp)
//  (const cGH **cctkGH,int *ext,double *f,double *favg, int *m);

void average_cpp(const cGH *cctkGH,int *ext,double *f,double *favg,
                 int m ) {
  
  using namespace std;

  // average over x: i,j,k -> i+1/2,j,k
  if (m==1) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) { 
	  int ip1 = i+1;
	  if (ip1==ext[0]) ip1 = i;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index1 = CCTK_GFINDEX3D(cctkGH,ip1,j,k);
	  favg[index] = 0.5*(f[index]+f[index1]);
	}
    return;
  }
  
  // average over y: i,j,k -> i,j+1/2,k
  if (m==2) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int jp1 = j+1;
	  if (jp1==ext[1]) jp1 = j;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index1 = CCTK_GFINDEX3D(cctkGH,i,jp1,k);
	  favg[index] = 0.5*(f[index]+f[index1]);
	}
    return;
  }

  // average over z: i,j,k -> i,j,k+1/2
  if (m==3) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int kp1 = k+1;
	  if (kp1==ext[2]) kp1 = k;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index1 = CCTK_GFINDEX3D(cctkGH,i,j,kp1);
	  favg[index] = 0.5*(f[index]+f[index1]);
	}
    return;
  }

  // average over x and y: i,j,k -> i+1/2,j+1/2,k 
  if (m==12) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int ip1 = i+1;
	  if (ip1==ext[0]) ip1 = i;
	  int jp1 = j+1;
	  if (jp1==ext[1]) jp1 = j;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,ip1,j,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,jp1,k);
	  int index11 = CCTK_GFINDEX3D(cctkGH,ip1,jp1,k);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  // average over x and z: i,j,k -> i+1/2,j,k+1/2
  if (m==13) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int ip1 = i+1;
	  if (ip1==ext[0]) ip1 = i;
	  int kp1 = k+1;
	  if (kp1==ext[2]) kp1 = k;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,ip1,j,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,j,kp1);
	  int index11 = CCTK_GFINDEX3D(cctkGH,ip1,j,kp1);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  // average over y and z: i,j,k -> i,j+1/2,k+1/2
  if (m==23) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int jp1 = j+1;
	  if (jp1==ext[1]) jp1 = j;
	  int kp1 = k+1;
	  if (kp1==ext[2]) kp1 = k;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,i,jp1,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,j,kp1);
	  int index11 = CCTK_GFINDEX3D(cctkGH,i,jp1,kp1);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  // average over x, y and z: i,j,k -> i+1/2,j+1/2,k+1/2 
  if (m==123) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int ip1 = i+1;
	  if (ip1==ext[0]) ip1 = i;
	  int jp1 = j+1;
	  if (jp1==ext[1]) jp1 = j;
	  int kp1 = k+1;
	  if (kp1==ext[2]) kp1 = k;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index100 = CCTK_GFINDEX3D(cctkGH,ip1,j,k);
	  int index110 = CCTK_GFINDEX3D(cctkGH,ip1,jp1,k);
	  int index101 = CCTK_GFINDEX3D(cctkGH,ip1,j,kp1);
	  int index010 = CCTK_GFINDEX3D(cctkGH,i,jp1,k);
	  int index001 = CCTK_GFINDEX3D(cctkGH,i,j,kp1);
	  int index011 = CCTK_GFINDEX3D(cctkGH,i,jp1,kp1);
	  int index111 = CCTK_GFINDEX3D(cctkGH,ip1,jp1,kp1);
	  favg[index] = 0.125*(f[index]+f[index100]+f[index110]+f[index101]+f[index010]+f[index001]+f[index011]+f[index111]);
	}
    return;
  }

  //average over x and y: i,j,k -> i-1/2,j+1/2,k
  if (m==92) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int im1 = i-1;
	  if (im1==-1) im1 = 0;
	  int jp1 = j+1;
	  if (jp1==ext[1]) jp1 = j;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,im1,j,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,jp1,k);
	  int index11 = CCTK_GFINDEX3D(cctkGH,im1,jp1,k);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;     
  }

  //average over x and z: i,j,k -> i-1/2,j,k+1/2
  if (m==93) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int im1 = i-1;
	  if (im1==-1) im1 = 0;
	  int kp1 = k+1;
	  if (kp1==ext[2]) kp1 = k;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,im1,j,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,j,kp1);
	  int index11 = CCTK_GFINDEX3D(cctkGH,im1,j,kp1);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;     
  }

  // average over x and y: i,j,k -> i+1/2,j-1/2,k
  if (m==18) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int ip1 = i+1;
	  if (ip1==ext[0]) ip1 = i;
	  int jm1 = j-1;
	  if (jm1==-1) jm1 = 0;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,ip1,j,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,jm1,k);
	  int index11 = CCTK_GFINDEX3D(cctkGH,ip1,jm1,k);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  // average over y and z: i,j,k -> i,j-1/2,k+1/2
  if (m==83) {
#pragma omp parallel for 
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int jm1 = j-1;
	  if (jm1==-1) jm1 = 0;
	  int kp1 = k+1;
	  if (kp1==ext[2]) kp1 = k;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,i,jm1,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,j,kp1);
	  int index11 = CCTK_GFINDEX3D(cctkGH,i,jm1,kp1);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  // average over x and z: i,j,k -> i+1/2,j,k-1/2
  if (m==17) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int ip1 = i+1; 
	  if (ip1==ext[0]) ip1 = i;
	  int km1 = k-1;
	  if (km1==-1) km1 = 0;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,ip1,j,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,j,km1);
	  int index11 = CCTK_GFINDEX3D(cctkGH,ip1,j,km1);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  // average over y and z: i,j,k -> i,j+1/2,k-1/2
  if (m==27) {
#pragma omp parallel for
    for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
	  int jp1 = j+1;
	  if (jp1==ext[1]) jp1 = j;
	  int km1 = k-1;
	  if (km1==-1) km1 = 0;
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int index10 = CCTK_GFINDEX3D(cctkGH,i,jp1,k);
	  int index01 = CCTK_GFINDEX3D(cctkGH,i,j,km1);
	  int index11 = CCTK_GFINDEX3D(cctkGH,i,jp1,km1);
	  favg[index] = 0.25*(f[index]+f[index10]+f[index01]+f[index11]);
	}
    return;
  }

  printf("Error: m=%d,  averaging scheme not supported.",m);
  exit(1);

}

//extern "C" void CCTK_FCALL CCTK_FNAME(average_cpp)
//  (const cGH **cctkGH,int *ext,double *f,double *favg, int *m) {
//  average_cpp(*cctkGH,ext,f,favg,*m);
//}

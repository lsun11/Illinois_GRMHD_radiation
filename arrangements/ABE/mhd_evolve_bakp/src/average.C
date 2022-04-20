#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <math.h>
#include "cctk.h"

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

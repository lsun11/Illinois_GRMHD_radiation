#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <strings.h>
#include "nonKeq0_soln_psi_gauge.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "nr.h"
#define NR_END 1
#define FREE_ARG char*

//===============================================
// *Produces a data file that can be used to plot
//  the value of the field variables versus the
//  radial coordinate or areal radius
//===============================================

using namespace std;

double output_1var_1radius(double r, int N, double r_EH, double r_outr_bdry, Vec_DP * TA, Vec_DP *varIA, Vec_DP *varOA, bool & legal);
double fourth_order_smooth_extrap(double *r_pts,double *y,double r_output);
void polint_statpunc(double xa[], double ya[], int n, double x, double *y, double *dy);
double *vector_statpunc(long nl, long nh);
void free_vector_statpunc(double *v, long nl, long nh);
void nrerror_statpunc(char error_text[]);
double trans(double r, double r_EH, double r_outr_bdry, int & t);    // translate r into s
double vecMul(Vec_DP a, Vec_DP b);  //Multiply vectors a and b
void makeT(double s,Vec_DP * TA);      // find values of Cheb. polys. at point s 

extern "C" void CCTK_FCALL CCTK_FNAME(stationary_puncture_id_nonKeq0_soln_set_psi_gauge)
  (const cGH **cctkGH, int *nghostzones,int *cctk_lsh,
   double *x, double *y, double *z,
   double *PhysicalRadius, double *r,
   double *lapse,double *shiftx,double *shifty,double *shiftz,double *betar,double *psi);

void stationary_puncture_id_nonKeq0_soln_set_psi_gauge(const cGH *cctkGH, int *nghostzones,int *cctk_lsh,
					 double *x, double *y, double *z,
					 double *PhysicalRadius, double *r,
					 double *lapse,double *shiftx,double *shifty,double *shiftz,double *betar,double *psi) {
  DECLARE_CCTK_PARAMETERS;

  int colpts;                // holds number of collocation points used
  int N;                     // highest order Chebychev polynomial used
  double r_EH;                 // radial coordinate of event horizon
  double r_outr_bdry;                 // radial coordinate of shell
  Vec_DP * cFacCoOA;         // holds values of Chebychev coefficients
  Vec_DP * lapsCoOA;         // for conformal factor, lapse and shift
  Vec_DP * shifCoOA;         // for outer solutions (O suffix before A)
  Vec_DP * cFacCoIA;         // and inner solutions (I suffix before A)
  Vec_DP * lapsCoIA;         // 
  Vec_DP * shifCoIA;         //
  Vec_DP * TA;               // holds values of Cheb. polynomials at point

  ifstream file;                                     // Set up input stream
  char filename[30];
  sprintf(filename,"prclSoln26.txt");
  printf("Opening file %s to read initial data...\n",filename);
  file.open(filename);
  if(!file) {
    cerr << "\a Can't open " << filename << " for input." << endl;
    exit(1);
  }

  bool start = false;                                // Ignore lines that
  char e = 'e';                                      // begin with #
  while(!start){                                     //
    file.get(e);                                     //
    if(e == '#'){                                    // 
      while(e != '\n'){                              // 
        file.get(e);                                 //
      }                                              //
    }else start=true;                                //
  }                                                  //

  file >> colpts;                                    // Read in colpts, r_EH, r_outr_bdry
  N=colpts-1;                                        //
  file >> r_EH;                                        //
  file >> r_outr_bdry;                                        //
                                                     //

  cFacCoOA          = new Vec_DP(colpts);            // Initialize vectors
  lapsCoOA          = new Vec_DP(colpts);            //
  shifCoOA          = new Vec_DP(colpts);            //
  cFacCoIA          = new Vec_DP(colpts);            //
  lapsCoIA          = new Vec_DP(colpts);            //
  shifCoIA          = new Vec_DP(colpts);            //
  TA                = new Vec_DP(colpts);            //
  Vec_DP &cFacCoO   = *cFacCoOA;                     //
  Vec_DP &lapsCoO   = *lapsCoOA;                     //
  Vec_DP &shifCoO   = *shifCoOA;                     // 
  Vec_DP &cFacCoI   = *cFacCoIA;                     //
  Vec_DP &lapsCoI   = *lapsCoIA;                     //
  Vec_DP &shifCoI   = *shifCoIA;                     // 
  Vec_DP &T         = *TA;

  for(int n=0;n<colpts;n++){                         // Read in vector data
    file >> cFacCoI[n];                              //
    file >> lapsCoI[n];                              //
    file >> shifCoI[n];                              //
    file >> cFacCoO[n];                              //
    file >> lapsCoO[n];                              //
    file >> shifCoO[n];                              //
  }                                                  //


  double rEH  = 0.8304040892;           // radial coordinate of event horizon
  double rMin = 0.0;                     // Minimum r tested
  //double rMin = .5;                     // Minimum r tested
  double rMax = 10;                     // Maximum r tested
  int res = 1001;                       // Number of evenly-spaced r tested
  double step = (rMax-rMin)/(res-1);    // Interval between tested r
  
  //Vec_DP r(res);                        // Vector holding r values tested
  Vec_BOOL legal(res);                  // r > r_{event horizon} ?  
  //Vec_DP confFact(res);                 // Vectors holding the values
  // Vec_DP lapse(res);                    // of the field variables
  //Vec_DP shift(res);                    // at tested r

  //Here we construct the coefficients for smooth junk extrapolation.
  // We use extrapolation function of the form: y = a r^4 + b r^2 + d
  //  That way, the derivative at r=0 is zero.
  //  Also, the connection to actual data is twice differentiable smooth
  double r_extrap[3];
  double confFact_y[3],shift_y[3],lapse_y[3];
  double blah_y[3];

  int n=8;
  double *r_extrap2=vector_statpunc(1,n);
  double *confFact_y2=vector_statpunc(1,n);
  double *shift_y2=vector_statpunc(1,n);
  double *lapse_y2=vector_statpunc(1,n);
  
  double dr = 0.02;
  /*
  for (int i=0;i<3;i++) {
    r_extrap[i] = rEH + i*dr;
    confFact_y[i] = output_1var_1radius(r_extrap[i],N,r_EH,r_outr_bdry,TA,cFacCoIA,cFacCoOA,legal[i]);
    lapse_y[i]    = output_1var_1radius(r_extrap[i],N,r_EH,r_outr_bdry,TA,lapsCoIA,lapsCoOA,legal[i]);
    shift_y[i]    = output_1var_1radius(r_extrap[i],N,r_EH,r_outr_bdry,TA,shifCoIA,shifCoOA,legal[i]);
  }
  */

  for (int i=1;i<=n;i++) {
    r_extrap2[i] = rEH + i*dr;
    confFact_y2[i] = output_1var_1radius(r_extrap2[i],N,r_EH,r_outr_bdry,TA,cFacCoIA,cFacCoOA,legal[i]);
    lapse_y2[i]    = output_1var_1radius(r_extrap2[i],N,r_EH,r_outr_bdry,TA,lapsCoIA,lapsCoOA,legal[i]);
    shift_y2[i]    = output_1var_1radius(r_extrap2[i],N,r_EH,r_outr_bdry,TA,shifCoIA,shifCoOA,legal[i]);
  }

  printf("Next, we sum over spectral coefficients to get psi, lapse, and shift (beta^r) at each gridpoint\n");
  printf("   Note that we fill the excised BH with smooth junk data during this procedure.\n");
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double r_phys = PhysicalRadius[index]; 
    double r_fish = r[index];

    psi[index]   = output_1var_1radius(r_phys,N,r_EH,r_outr_bdry,TA,cFacCoIA,cFacCoOA,legal[i]);
    lapse[index] = output_1var_1radius(r_phys,N,r_EH,r_outr_bdry,TA,lapsCoIA,lapsCoOA,legal[i]) - 1.0;
    betar[index] = output_1var_1radius(r_phys,N,r_EH,r_outr_bdry,TA,shifCoIA,shifCoOA,legal[i]);
    if(!legal[i]) {
      double dy;
      polint_statpunc(r_extrap2,confFact_y2,n,r_phys,&psi[index],&dy);
      polint_statpunc(r_extrap2,lapse_y2,n,r_phys,&lapse[index],&dy);
      polint_statpunc(r_extrap2,shift_y2,n,r_phys,&betar[index],&dy);

      if(psi[index]<0.0) psi[index] = fabs(psi[index]);

      lapse[index] = fabs(lapse[index]);
      lapse[index] = lapse[index] - 1.0;
      
      //printf("hi r=%e: %e %e %e\n",r_phys,psi[index],lapse[index],betar[index]);

      /*
      psi[index]   = fourth_order_smooth_extrap(r_extrap,confFact_y,r_phys);
      lapse[index] = fourth_order_smooth_extrap(r_extrap,lapse_y,r_phys) - 1.0;
      betar[index] = fourth_order_smooth_extrap(r_extrap,shift_y,r_phys);
      */
    }

    shiftx[index] = betar[index]*x[index]/r_fish;
    shifty[index] = betar[index]*y[index]/r_fish;
    shiftz[index] = betar[index]*z[index]/r_fish;

    if(k%2==0 && j==0 && i==0) printf("Reading input file: %d%% finished...\n",(int)((double)k/cctk_lsh[2]*100.0));
  }

  printf("Finished setting shift, lapse, psi!\n");
  file.close();
  delete cFacCoOA;
  delete lapsCoOA;
  delete shifCoOA;
  delete cFacCoIA;
  delete lapsCoIA;
  delete shifCoIA;
  delete TA;
}
double fourth_order_smooth_extrap(double *r_pts,double *y,double r_output) {
  double y1,y2,y3;
  double r1,r2,r3;
  y1=y[0]; y2=y[1]; y3=y[2];
  r1=r_pts[0]; r2=r_pts[1]; r3=r_pts[2];
  double a = ( (y3-y2)/(r3*r3-r2*r2) -  (y1-y2)/(r1*r1-r2*r2) ) / ( (r3*r3*r3*r3-r2*r2*r2*r2)/(r3*r3-r2*r2) - (r1*r1*r1*r1-r2*r2*r2*r2)/(r1*r1-r2*r2) );
  double b = (y1-y2)/(r1*r1 - r2*r2) - a * (r1*r1*r1*r1-r2*r2*r2*r2)/(r1*r1-r2*r2);
  double d = 0.5 * (y2 + y3 - a * (r2*r2*r2*r2 + r3*r3*r3*r3) - b *(r2*r2 + r3*r3));

  //printf("hey! %e %e %e\n",a,b,d);

  return a*r_output*r_output*r_output*r_output + b*r_output*r_output + d;
}

extern "C" void CCTK_FCALL CCTK_FNAME(stationary_puncture_id_nonKeq0_soln_set_psi_gauge)
  (const cGH **cctkGH, int *nghostzones,int *cctk_lsh,
   double *x, double *y, double *z,
   double *PhysicalRadius, double *r,
   double *lapse,double *shiftx,double *shifty,double *shiftz,double *betar,double *psi) {
  stationary_puncture_id_nonKeq0_soln_set_psi_gauge(*cctkGH,nghostzones,cctk_lsh,
						    x,y,z,
						    PhysicalRadius,r,
						    lapse, shiftx, shifty, shiftz, betar, psi);
}

void polint_statpunc(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;

  dif=fabs(x-xa[1]);
  c=vector_statpunc(1,n);
  d=vector_statpunc(1,n);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) nrerror_statpunc("Error in routine polint_statpunc");
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free_vector_statpunc(d,1,n);
  free_vector_statpunc(c,1,n);
}

double *vector_statpunc(long nl, long nh)
  /* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror_statpunc("allocation failure in vector_statpunc()");
  return v-nl+NR_END;
}

void nrerror_statpunc(char error_text[])
  /* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error in stationary puncture...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

void free_vector_statpunc(double *v, long nl, long nh)
  /* free a double vector allocated with vector_statpunc() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

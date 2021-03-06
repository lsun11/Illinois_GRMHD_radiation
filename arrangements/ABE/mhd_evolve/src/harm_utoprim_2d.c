/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_2d.c: 
---------------

    Uses the 2D method: 
       -- solves for two independent variables (W,v^2) via a 2D
          Newton-Raphson method 
       -- can be used (in principle) with a general equation of state. 

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

#include "harm_u2p_util.h"

#define NEWT_DIM 2

// Declarations: 
static FTYPE vsq_calc(FTYPE W,FTYPE &Bsq,FTYPE &QdotBsq,FTYPE &Qtsq,FTYPE &Qdotn,FTYPE &D);
static int Utoprim_new_body(FTYPE U[], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[], int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE &T_fluid, FTYPE *rho_tab, FTYPE *P_tab,  FTYPE *eps_tab,  FTYPE *k_tab,  FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics);

//static int Utoprim_new_body(FTYPE U[], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[]);

static int general_newton_raphson( FTYPE x[], int n, void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int,FTYPE &, FTYPE &, FTYPE &, FTYPE &,FTYPE &,FTYPE &, FTYPE &, FTYPE & , int &, int &, FTYPE &,FTYPE *, FTYPE *, FTYPE *,FTYPE *, FTYPE *, int &, int &),FTYPE &Bsq,FTYPE &QdotBsq,FTYPE &Qtsq,FTYPE &Qdotn,FTYPE &D,FTYPE &T_fluid, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *gamma_tab, FTYPE *eps_tab, FTYPE *k_tab, int &enable_OS_collapse, int &compute_microphysics);
static void func_vsq( FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *f, FTYPE *df, int n,FTYPE &Bsq,FTYPE &QdotBsq,FTYPE &Qtsq,FTYPE &Qdotn,FTYPE &D, FTYPE &T_fluid, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *gamma_tab, FTYPE *eps_tab, FTYPE *k_tab, int &enable_OS_collapse, int &compute_microphysics);
static FTYPE x1_of_x0(FTYPE x0, FTYPE &Bsq, FTYPE &QdotBsq, FTYPE &Qtsq, FTYPE &Qdotn, FTYPE &D ) ;
//static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq, FTYPE &D, FTYPE &P_cold, FTYPE &eps_cold) ;
static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq, FTYPE &T_fluid, FTYPE &D, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *eps_tab, FTYPE *k_tab, FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics);

static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq, int &compute_microphysics);
//static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq, FTYPE &D, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma,FTYPE *rho_tab, FTYPE *P_tab, FTYPE *gamma_tab, FTYPE *eps_tab);
static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq, FTYPE &D, FTYPE &P_cold, FTYPE &eps_cold, FTYPE &T_fluid, int &neos, int &ergo_star, FTYPE &ergo_sigma,FTYPE *rho_tab, FTYPE *P_tab, FTYPE *eps_tab, FTYPE *k_tab, FTYPE *gamma_tab, int &enable_OS_collapse,int &compute_microphysics);

static FTYPE Gamma_cold (FTYPE rho0, int &neos, FTYPE *rho_tab, FTYPE *gamma_tab, FTYPE *eps_tab, FTYPE eps);


extern void compute_pcold_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
			      int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab,
			      double *eps_tab, double *k_tab,
				      double *gamma_tab,int &enable_OS_collapse);

/**********************************************************************/
/******************************************************************

  Utoprim_2d():
  
  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  Note that Greek
     indices run 0,1,2,3 and Latin indices run 1,2,3 (spatial only).


              /  rho u^t            \
         U =  |  T^t_t + rho u^t    |  sqrt(-det(g_{\mu\nu}))
              |  T^t_i              |
              \  B^i                /

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

int Utoprim_2d(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
	       FTYPE gdet, FTYPE prim[NPR], int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE &T_fluid, FTYPE *rho_tab, FTYPE *P_tab,  FTYPE *eps_tab,  FTYPE *k_tab,  FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics)
{

  FTYPE U_tmp[NPR], U_tmp2[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE alpha;


  if( U[0] <= 0. ) { 
    return(-100);
  }

  /* First update the primitive B-fields */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] / gdet ;

  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-gcon[0][0]);
  
  /* Transform the CONSERVED variables into the new system */
  U_tmp[RHO] = alpha * U[RHO] / gdet;
  U_tmp[UU]  = alpha * (U[UU] - U[RHO])  / gdet ;
  for( i = UTCON1; i <= UTCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }

  /* Transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  if(U_tmp[RHO] > 1.0e40) { printf("inside Utoprim_2d, U_tmp[RHO] = alpha * U[RHO] / gdet is too large to be correct. U_tmp[RHO], alpha, U[RHO], gdet are %e, %e, %e, %e \n", U_tmp[RHO], alpha, U[RHO], gdet);}


  ret = Utoprim_new_body(U_tmp, gcov, gcon, gdet, prim_tmp, neos, ergo_star, ergo_sigma, T_fluid, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse, compute_microphysics);
     

  /* Transform new primitive variables back if there was no problem : */ 
  if( ret == 0 || ret == 5 || ret==101 ) {
    for( i = 0; i < BCON1; i++ ) {
      prim[i] = prim_tmp[i];
    }
  }

  return( ret ) ;

  //printf("after harm_utoprim_2d.c, ret value is %d", ret);
}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu  |
             \  alpha B^i        /



               /    rho        \
	prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;

STATISTICS: Lots of unphysical roots!
 # FAILURES REASON
    127     -100 
  41570     5   
    392     101 
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], 
			    FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[NPR], int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE &T_fluid, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *eps_tab,  FTYPE *k_tab,  FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics)
{

  FTYPE x_2d[NEWT_DIM];
  FTYPE QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qsq,Qtcon[NDIM];
  FTYPE rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq,tmpdiff ;
  int i,j, n, retval, i_increase ;

  //TESTING vv
  double P_cold, eps, eps_cold;
  //TESTING ^^


  n = NEWT_DIM ;

  // Assume ok initially:
  retval = 0;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  
  raise_g(Qcov,gcon,Qcon) ;



  FTYPE Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  QdotB = 0. ;
  for(i=0;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  FTYPE QdotBsq = QdotB*QdotB ;

  ncov_calc(gcon,ncov) ;
  raise_g(ncov,gcon,ncon);

  FTYPE Qdotn = Qcon[0]*ncov[0] ;

  Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;

  FTYPE Qtsq = Qsq + Qdotn*Qdotn ;

  FTYPE D = U[RHO] ;

  /* calculate W from last timestep and use for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++){
      utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1];
      if (utsq > UTSQ_TOO_BIG) {printf ("i=,j=,gcov[i][j],prim[UTCON1+i-1]=,prim[UTCON1+j-1]=, %d, %d, %e, %e, %e\n", i, j, gcov[i][j], prim[UTCON1+i-1],prim[UTCON1+j-1]);}
    }

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  //if (utsq==0.0) {printf ("utsq= %e\n", utsq);}


  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	

  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0, T_fluid, u, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse, compute_microphysics) ;  

  w = rho0 + u + p ;

  W_last = w*gammasq ;


  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
	    - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  //printf("before genereal_newton_raphson, W_last is %e\n",W_last);

  if(W_last == 0.0) {printf("W_last = w*gammasq is 0, w=%e, gammasq=%e, utsq=%e, \n", w, gammasq, utsq);}
  // Calculate W and vsq: 
  x_2d[0] =  fabs( W_last );
  x_2d[1] = x1_of_x0( W_last , Bsq,QdotBsq,Qtsq,Qdotn,D) ;


  if(x_2d[0] > 1.0e60) { printf("before general N-R , x_2d[0] = fabs( W_last ) is too large to be correct. w, rho0, u, p,  and  are %e, %e, %e, %e \n", w, rho0, u,p);}
  //if(x_2d[1]==1.0) { printf("before general N-R , vsq is 1. \n");}



  //printf("before general N-R , W is %.17g\n , rho0 is %.17g\n, u is %.17g\n, p is %.17g, w is %.17g, D is %.17g, vsq is %.17g \n", x_2d[0], rho0, u , p, w, D, x_2d[1]);

  retval = general_newton_raphson( x_2d, n, func_vsq, Bsq,QdotBsq,Qtsq,Qdotn,D, T_fluid, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, gamma_tab, eps_tab, k_tab, enable_OS_collapse, compute_microphysics); 

  //printf("after general N-R , W is %.17g\n , rho0 is %.17g\n, u is %.17g\n, p is %.17g, w is %.17g,  D is %.17g, vsq is %.17g \n", x_2d[0], rho0, u , p, w, D, x_2d[1]);

  W = x_2d[0];
  vsq = x_2d[1];
	
  if(isnan(W)) { printf("after general N-R , W is nan. \n");}
  if(isnan(vsq)) { printf("after general N-R , vsq is nan. \n");}



  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }

  // Calculate v^2:
  if( vsq >= 1. ) {
    vsq = 1.-2.e-16;
    //retval = 4;
    //return(retval) ;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;

  w = W * (1. - vsq) ;
 
  // W = rho0 + rho0*eps + P = rho0 + (P-P_cold)/(Gamma_th-1) + rho0*eps_cold + P (see harm_primitives_lowlevel, serveral lines above "CALL HARM PRIMITIVES SOLVER")
  //   = rho0 + Gamma_th*P/(Gamma_th-1) + eps_cold * rho0 - P_cold/(Gamma_th-1)
  //--> (W - rho0) - (eps_cold*rho0 - P_cold/(Gamma_th-1)) = Gamma_th/(Gamma_th-1) * P, where we get P.

  p = pressure_rho0_w(rho0,T_fluid, w,P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse, compute_microphysics) ;
  //w is possibly changed in pressure_rho0_w, reassgin W:
  W = w / (1. - vsq) ;
  //  p = pressure_rho0_w(rho0,w, P_cold, eps_cold) ;
  u = w - (rho0 + p) ; // u = rho0 eps, w = rho0 h
 

  if( (rho0 <= 0.) || (u <= 0.) ) {
    // User may want to handle this case differently, e.g. do NOT return upon 
    // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:

    retval = 5;
    //return(retval) ;
  }

  prim[RHO] = rho0 ;
  prim[UU] = u ;

  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;

  //  if(isnan(prim[UTCON1])) { printf("Inside harm_utoprim_2d, prim[UTCON1] is nan. W, Bsq, QdotB are %e,%e,%e \n", W, Bsq, QdotB);}
	
  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* done! */
  return(retval) ;

}


/**********************************************************************/ 
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static FTYPE vsq_calc(FTYPE W,FTYPE &Bsq,FTYPE &QdotBsq,FTYPE &Qtsq,FTYPE &Qdotn,FTYPE &D)
{
  FTYPE Wsq,Xsq, a;
	
	Wsq = W*W ;
	Xsq = (Bsq + W) * (Bsq + W);

	a = ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
	return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );
	if (isnan(a)){printf("inside vsq_calc, vsq is nan, W=%e, Wsq=%e, Qtsq=%e, QdotBsq=%e, Bsq=%e \n", W, Wsq, Qtsq, QdotBsq, Bsq);}
}


/********************************************************************

  x1_of_x0(): 
           
    -- calculates v^2 from W  with some physical bounds checking;
    -- asumes x0 is already physical
    -- makes v^2 physical  if not;

*********************************************************************/

static FTYPE x1_of_x0(FTYPE x0, FTYPE &Bsq, FTYPE &QdotBsq, FTYPE &Qtsq, FTYPE &Qdotn, FTYPE &D ) 
{
  FTYPE x1,vsq;
  FTYPE dv = 1.e-15;
  
    vsq = fabs(vsq_calc(x0,Bsq,QdotBsq,Qtsq,Qdotn,D)) ; // guaranteed to be positive 


    if(vsq==1.0) {printf("---Inside x1_of_x0, vsq is 1. \n");}
    if(isnan(vsq)) {printf("---Inside x1_of_x0, vsq is nan, %e, %e,%e,%e,%e,%e \n", x0,Bsq,QdotBsq,Qtsq,Qdotn,D);}
  return( ( vsq >= 1. ) ? (1.0 - dv) : vsq   ); 

}

/********************************************************************

  validate_x(): 
           
    -- makes sure that x[0,1] have physical values, based upon 
       their definitions:
    
*********************************************************************/

static void validate_x(FTYPE x[2], FTYPE x0[2] ) 
{
  
  FTYPE dv = 1.e-15;

  /* Always take the absolute value of x[0] and check to see if it's too big:  */ 
  x[0] = fabs(x[0]);
  x[0] = (x[0] > W_TOO_BIG) ?  x0[0] : x[0];
  

  x[1] = (x[1] < 0.) ?   0.       : x[1];  /* if it's too small */
  x[1] = (x[1] >= 1.) ?  (1. - dv) : x[1];  /* if it's too big   */

  return;

}

/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, 
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
					   FTYPE [][NEWT_DIM], FTYPE *, 
					   FTYPE *, int,FTYPE &,FTYPE &,FTYPE &,
					   FTYPE &,FTYPE &, FTYPE &, FTYPE &,FTYPE &,
					   int &, int &, FTYPE &, FTYPE *, FTYPE *,
					   FTYPE *,FTYPE *, FTYPE *, int &, int &),FTYPE &Bsq,FTYPE &QdotBsq,FTYPE &Qtsq,FTYPE &Qdotn,FTYPE &D, FTYPE &T_fluid, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *gamma_tab, FTYPE *eps_tab, FTYPE *k_tab, int &enable_OS_collapse, int &compute_microphysics)
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM];
  FTYPE resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  vsq_old = vsq = W = W_old = 0.;
  n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

     if(x[0] > 1.0e60) { printf("before func_vsq, x[0] is too large to be correct. \n");}
     if(isnan(x[1])) { printf("before func_vsq, x[1] is nan, x[0]=%e \n", x[0]);}

     (*funcd) (x, dx, resid, jac, &f, &df, n, Bsq,QdotBsq,Qtsq,Qdotn,D, T_fluid, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, gamma_tab, eps_tab, k_tab, enable_OS_collapse, compute_microphysics);  /* returns with new dx, f, df */
      
    // if(isnan(x[0])) { printf("after func_vsq, W is nan. \n");}
    // if(x[1]==1.0) { printf("after func_vsq, vsq is 1. \n");}


    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    validate_x( x, x_old ) ;


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) 
	|| (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)

  //   printf ("inside N-R solver, f= %e, df = %e", f, df);
    /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) ||  (finite(df)==0) ) {
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
    //    printf("Newton Solver returns 1 !!!!! %d %e %e %e %e\n",n_iter,f,df,errx,MIN_NEWT_TOL);
    return(1);
  } 
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }

  return(0);

}



/**********************************************************************/
/*********************************************************************************
   func_vsq(): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W,vsq here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
*********************************************************************************/

static void func_vsq(FTYPE x[], FTYPE dx[], FTYPE resid[], 
		     FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n,
		     FTYPE &Bsq,FTYPE &QdotBsq,FTYPE &Qtsq,FTYPE &Qdotn,FTYPE &D,
		     FTYPE &T_fluid, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma,
		     FTYPE *rho_tab, FTYPE *P_tab, FTYPE *gamma_tab, FTYPE *eps_tab, FTYPE *k_tab,
		     int &enable_OS_collapse, int &compute_microphysics)
{
  
  FTYPE  W, vsq, Wsq, p_tmp, dPdvsq, dPdW, temp, detJ,tmp2,tmp3;
  FTYPE t11;
  FTYPE t16;
  FTYPE t18;
  FTYPE t2;
  FTYPE t21;
  FTYPE t23;
  FTYPE t24;
  FTYPE t25;
  FTYPE t3;
  FTYPE t35;
  FTYPE t36;
  FTYPE t4;
  FTYPE t40;
  FTYPE t9;

  W = x[0];
  vsq = x[1];

  
  Wsq = W*W;

  p_tmp  = pressure_W_vsq( W, vsq , T_fluid, D, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse, compute_microphysics);
  // W is possibly changed, reassign the values.
  x[0] = W;
  Wsq = W*W;
  dPdW   = dpdW_calc_vsq( W, vsq, compute_microphysics );
  //  dPdvsq = dpdvsq_calc( W, vsq, D, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab,gamma_tab, eps_tab);
  dPdvsq = dpdvsq_calc( W, vsq, D, P_cold, eps_cold, T_fluid, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse, compute_microphysics);


  
  if(isnan(p_tmp)) { printf("inside func_vsq, p_tmp is nan. \n");}
  if(isnan(dPdW)) { printf("inside func_vsq ,dPdW is nan. \n");}
  if(isnan(dPdvsq)) { printf("inside func_vsq, dPdvsq is nan. Here W, vsq, are %e, %e, \n", W ,vsq); }
  

  // These expressions were calculated using Mathematica, but made into efficient 
  // code using Maple.  Since we know the analytic form of the equations, we can 
  // explicitly calculate the Newton-Raphson step: 

  t2 = -0.5*Bsq+dPdvsq;
  t3 = Bsq+W;
  t4 = t3*t3;
  t9 = 1/Wsq;
  t11 = Qtsq-vsq*t4+QdotBsq*(Bsq+2.0*W)*t9;
  t16 = QdotBsq*t9;
  t18 = -Qdotn-0.5*Bsq*(1.0+vsq)+0.5*t16-W+p_tmp;
  t21 = 1/t3;
  t23 = 1/W;
  t24 = t16*t23;
  t25 = -1.0+dPdW-t24;
  t35 = t25*t3+(Bsq-2.0*dPdvsq)*(QdotBsq+vsq*Wsq*W)*t9*t23;
  t36 = 1/t35;
  dx[0] = -(t2*t11+t4*t18)*t21*t36;
  t40 = (vsq+t24)*t3;
  dx[1] = -(-t25*t11-2.0*t40*t18)*t21*t36;
  detJ = t3*t35;
  jac[0][0] = -2.0*t40;
  jac[0][1] = -t4;
  jac[1][0] = t25;
  jac[1][1] = t2;
  resid[0] = t11;
  resid[1] = t18;


  *df = -resid[0]*resid[0] - resid[1]*resid[1];

  *f = -0.5 * ( *df );

}
/********************************************************************** 
 ********************************************************************** 
   
 The following routines specify the equation of state.  All routines 
  above here should be indpendent of EOS.  If the user wishes 
  to use another equation of state, the below functions must be replaced 
  by equivalent routines based upon the new EOS. 

 **********************************************************************
**********************************************************************/

/**********************************************************************/
/********************************************************************** 
  pressure_W_vsq():  
  
        -- Gamma-law equation of state;
        -- pressure as a function of W, vsq, and D:
 **********************************************************************/
static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq, FTYPE &T_fluid, FTYPE &D, FTYPE &P_cold, FTYPE &eps_cold, int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *eps_tab, FTYPE *k_tab, FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics) 
{
  FTYPE gtmp;
  FTYPE P_total;

  gtmp = 1. - vsq;

  FTYPE rho0= D*sqrt(gtmp);
  FTYPE w = W * gtmp;
  //  FTYPE P_cold, eps_cold;
  compute_pcold_epscold_cpp(rho0, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse);
  
  if (compute_microphysics == 1){
    compute_P_T_microphys_w(P_total, w, T_fluid, P_cold, eps_cold, rho0);
    // compute_P_T_microphys_w possibly changes w
    W = w/gtmp;
    return P_total;
  }
  else{
    return(  (GAMMA_th - 1.) * ( W * gtmp  -  D*(1.0 + eps_cold) * sqrt(gtmp)) + P_cold) / GAMMA_th ; 
  }
}


/**********************************************************************/
/********************************************************************** 
  dpdW_calc_vsq(): 
 
      -- partial derivative of pressure with respect to W;
 **********************************************************************/
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq, int &compute_microphysics)
{

  //return( (GAMMA - 1.) * (1. - vsq) /  GAMMA ) ;
  //if(isnan(W)) {printf("---Inside dpdW_calc_vsq, W is nan, GAMMA_th =%e, \n", GAMMA_th);}
  
  if(vsq==1.0) {printf("---Inside dpdW_calc_vsq, vsq is 1, GAMMA_th =%e, \n", GAMMA_th);}
  if (compute_microphysics == 1){
    return (0.25*( 1. - vsq)) ;
  }else{
    return( (GAMMA_th - 1.) * (1. - vsq) /  GAMMA_th ) ;
  }
}

/**********************************************************************/
/********************************************************************** 
  dpdvsq_calc(): 
 
      -- partial derivative of pressure with respect to vsq
 **********************************************************************/
static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq, FTYPE &D, FTYPE &P_cold, FTYPE &eps_cold, FTYPE &T_fluid, int &neos, int &ergo_star, FTYPE &ergo_sigma, FTYPE *rho_tab, FTYPE *P_tab, FTYPE *eps_tab, FTYPE *k_tab, FTYPE *gamma_tab, int &enable_OS_collapse, int &compute_microphysics)
{
  FTYPE rho0, gtmp, eps, gamma_cold;
  gtmp = 1. - vsq;
  rho0 = D*sqrt(gtmp);
  gamma_cold = Gamma_cold(rho0, neos, rho_tab, gamma_tab, eps_tab, eps);

  FTYPE drho_dvsq,gammas;
  drho_dvsq = - 0.5 * D / sqrt(gtmp);
  gammas = gamma_cold*(gamma_cold-GAMMA_th)/(gamma_cold-1.0);

  //  FTYPE P_cold, eps_cold;
  compute_pcold_epscold_cpp(rho0, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse);

  // if(vsq == 1.0) {printf("---Inside dpdvsq_calc, vsq is 1, GAMMA_th =%e, W =%e, vsq =%e, D =%e, P_cold =%e, eps_cold =%e \n", GAMMA_th, W, vsq, D, P_cold, eps_cold);}

  //Here we include P=rho EOS (ergo_star==1 && rho > rho_tab[neos-1]).
  if (compute_microphysics ==1 ){
    FTYPE m_n;
    FTYPE P_total;
    m_n = 1.244e-57; //neutron mass in code unit (km)
    FTYPE w = W * gtmp;
    compute_P_T_microphys_w(P_total, w, T_fluid, P_cold, eps_cold, rho0);
    return ( 0.25*( -W + 1.5*T_fluid/m_n*drho_dvsq ) );
  }else{
    if (ergo_star == 0 || (ergo_star == 1 && rho0 <= rho_tab[neos-1])){
      return ( -(GAMMA_th - 1.) * (drho_dvsq*(1.0+eps) + W ) + (gammas * P_cold * drho_dvsq/rho0 ))/GAMMA_th;
    }
    else {
      FTYPE h_tab, Con1, Con2;
      h_tab= 1.0+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1];
      Con1 = ergo_sigma * h_tab / pow (rho_tab[neos-1], ergo_sigma)* pow(rho0, ergo_sigma);
      Con2 = ergo_sigma * (1.0 + eps_tab[neos-1]) * rho_tab[neos-1] - P_tab[neos-1];  
      return (drho_dvsq * ( Con1 - (GAMMA_th - 1.0) * ( (1.0+eps_cold) + (Con1/rho0 - Con2/rho0/rho0 ) / (ergo_sigma + 1.0) ) ) - (GAMMA_th - 1.0) * W ) / GAMMA_th;
    }
  }
}


/*********************************************************************
 Gamma_cold():                                                                                                                                                                                            
       
     -- adiabatic index for cold matter according to denisty.
                                                                                                                                                 
**********************************************************************/

static FTYPE Gamma_cold (FTYPE rho0, int &neos, FTYPE *rho_tab, FTYPE *gamma_tab, FTYPE *eps_tab, FTYPE eps)
{

  //printf("inside Gamma_cold, eps_tab=%e,%e, %e,%e", eps_tab[0],eps_tab[1],eps_tab[2], eps_tab[3]);
  if(neos==1){
    eps = eps_tab[0];
    return gamma_tab[0];
  }

  else {
    for (int nn=1; nn<neos; nn++) {
      if (rho0<=rho_tab[nn] && rho0> rho_tab[nn-1]){
	eps = eps_tab[nn-1];
	return gamma_tab[nn-1];
      }
    }
      if (rho0> rho_tab[neos-1]) {
        eps = eps_tab [neos-1];
        return gamma_tab[neos-1]; 
    }
  }
  //eps_tab = 0;
  //return 1.0;
}




/****************************************************************************** 
             END   OF   UTOPRIM_2D.C
 ******************************************************************************/




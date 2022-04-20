 /*@@
   @file      Kerr.c
   @date      August 2000
   @author    John Baker
   @desc 
      Set up initial data for a spinning black holes
   @enddesc 
   @version $Header: /cactus/CactusEinstein/IDAnalyticBH/src/Kerr.c,v 1.13 2004/05/17 12:28:56 goodale Exp $
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "IDAnalyticBH.h"

static const char *rcsid = "$Header: /cactus/CactusEinstein/IDAnalyticBH/src/Kerr.c,v 1.13 2004/05/17 12:28:56 goodale Exp $";

CCTK_FILEVERSION(CactusEinstein_IDAnalyticBH_Kerr_c)


 /*@@
   @routine    KerrID
   @date       August 2000
   @author     John Baker
   @desc 
   Create Kerr Initital data
   @enddesc 
   @calls     
   @calledby   
   @history 
   @hdate Fri Apr 26 10:04:05 2002 @hauthor Tom Goodale
   @hdesc Changed to use new StaticConformal stuff
   @endhistory 
 
 @@*/
void KerrID(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  int i, do_lapse=0, do_shift=0;
  int npoints;
  CCTK_REAL rK,R,cth,cth_2,sth,sth_2,r_2,R_2,R_3,rhoradius,rhoradius_2,xx,yy,zz;
  CCTK_REAL Sigma,sqrt_Delta,p2,lapse,beta_phi,shift_phi;
  CCTK_REAL tmp, inv_psi;
  CCTK_REAL Phi4,fourPhi3,Phi,Chi2;
  CCTK_REAL Phi_R,Phi_RR,Phi_Rq,Phi_q,Phi_qq;
  CCTK_REAL Phi4_R,Phi4_RR,Phi4_Rq,Phi4_q,Phi4_qq;
  /* CCTK_REAL Delta,Krj,gRR,gqq,gjj; */
  CCTK_REAL KRj,Kqj;
  
  CCTK_REAL dRdx,dRdy,dRdz,dqdx,dqdy,dqdz,djdx,djdy;
  CCTK_REAL d2Rdxx,d2Rdxy,d2Rdxz,d2Rdyy,d2Rdyz,d2Rdzz;
  CCTK_REAL d2qdxx,d2qdxy,d2qdxz,d2qdyy,d2qdyz,d2qdzz;
  CCTK_REAL m=mass,a=a_Kerr,a_2=a*a,m2_a2=m*m-a_2;
  int make_conformal_derivs = 0;

  CCTK_VInfo(CCTK_THORNSTRING,
             "setting up Kerr initial data");

  /* total number of points on this processor */
  npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  if (CCTK_Equals(initial_lapse,"Kerr"))
  {
    do_lapse=1;    
    CCTK_INFO("Initialise with Kerr lapse");
  }
  if (CCTK_Equals(initial_shift,"Kerr"))
  {
    do_shift=1;    
    CCTK_INFO("Initialise with Kerr shift");
  }


  /* Check if we should create and store conformal factor stuff */
  if(CCTK_EQUALS(metric_type, "static conformal"))
  {
    if      (CCTK_EQUALS(conformal_storage,"factor"))
    {
      *conformal_state = 1;
      make_conformal_derivs = 0;
    }
    else if (CCTK_EQUALS(conformal_storage,"factor+derivs"))
    {
      *conformal_state = 2;
      make_conformal_derivs = 1;
    }
    else if (CCTK_EQUALS(conformal_storage,"factor+derivs+2nd derivs"))
    {
      *conformal_state = 3;
      make_conformal_derivs = 1;
    }
    else
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
"KerrID(): impossible value for conformal_storage=\"%s\"!",
                 conformal_storage);                            /*NOTREACHED*/
    }
  }      


  /* printf("npoints: %i\n",npoints); */
  for(i = 0; i < npoints; i++)
  {
    /* if((i/1000)*1000==i||i>35930)printf("i=%i\n",i); */
    /* Define coordinate functions */
    xx=x[i];  yy=y[i];  zz=z[i];
    rhoradius_2=xx*xx+yy*yy;
    rhoradius=sqrt(rhoradius_2);
    R_2=rhoradius_2+zz*zz;
    R=sqrt(R_2);
    R_3=R*R_2;
    cth=zz/R;
    cth_2=cth*cth;
    sth_2=rhoradius_2/R_2;
    sth=rhoradius/R;

    /* Special to Kerr */
    rK=R+m+m2_a2/4/R;
    r_2=rK*rK;
    sqrt_Delta=R-m2_a2/4/R;
    /*Delta=sqrt_Delta*sqrt_Delta;*/
    Sigma=r_2+a_2*cth_2;
    beta_phi=-2*m*rK*a*sth_2/Sigma;
    p2=a_2+r_2-a*beta_phi;
    /* drdR=sqrt_Delta/R; */
    lapse=sqrt_Delta/sqrt(p2);
    shift_phi=-2*m*rK*a/p2;
    
    /* Kerr metric in quasi-isotropic coordinates */
    /* ds^2=Phi4*(dR^2+R^2*dth^2+R^2*Chi^2*sin(th)^2*dphi^2) */
    Phi4=Sigma/R_2;
    Chi2=p2/Sigma;
    /* gRR=1; */
    /* gqq=R_2; */
    /* gjj=R_2*sth_2*Chi2; */
    
    /* extrinsic curvature */
    /* tmp=p2+a*beta_phi/2; */
    /* Krj=(beta_phi*rK)/(2/lapse/p2)*(1+(2/Sigma-1/R_2)*tmp); */
    /* Kqj=(-2*m*a*sth*cth/Sigma)/lapse*((r_2+a_2)/(p2*Sigma)*tmp-.5); */
    /* transform r to R */
    /* KRj=drdR*Krj; */
    tmp=sqrt(p2)*Sigma;
    KRj=a*m*sth_2/R/tmp*(r_2-a_2+2*r_2*(r_2+a_2)/Sigma);
    Kqj=sqrt_Delta*beta_phi*a_2*sth*cth/tmp;

    /* Conformal factor derivatives */
    {
      double drdR=sqrt_Delta/R;
      double d2rdR2=m2_a2/2/R_3;
      double dSigmadR=2*rK*drdR;
      double d2SigmadR2=2*(drdR*drdR+rK*d2rdR2);
      Phi4_R=dSigmadR/R_2-2*Sigma/R_3;
      Phi4_RR=d2SigmadR2/R_2-4*dSigmadR/R_3+6*Sigma/R_2/R_2;
    }
    /* Phi4_R =2/R_3*(rK*sqrt_Delta-Sigma); */
    /* Phi4_RR=(m2_a2*rK/R+2*(Delta+Sigma)-8*rK*sqrt_Delta)/R_2/R_2; */
    Phi4_q =-2*a_2*cth*sth/R_2;
    Phi4_Rq=-2*Phi4_q/R;
    Phi4_qq=2*a_2*(-cth_2+sth_2)/R_2;

    Phi=pow(Phi4,.25);
    fourPhi3=4*Phi*Phi*Phi;
    
    Phi_R =Phi4_R/fourPhi3;
    Phi_RR=(Phi4_RR-3*Phi4_R*Phi4_R/Phi4/4)/fourPhi3;
    Phi_q =Phi4_q/fourPhi3;
    Phi_Rq=(Phi4_Rq-3*Phi4_R*Phi4_q/Phi4/4)/fourPhi3;
    Phi_qq=(Phi4_qq-3*Phi4_q*Phi4_q/Phi4/4)/fourPhi3;

    /* Now we set the cactus variables */
    psi [i] = Phi;

    if(do_lapse)
    {
      lapm1[i]=lapse-1;
    }

    /* transform to cartesian coordinates */
    dRdx=xx/R;           
    dRdy=yy/R;            
    dRdz=zz/R;
    dqdx=xx*zz/rhoradius/R_2;  
    dqdy=yy*zz/rhoradius/R_2;  
    dqdz=-rhoradius/R_2;
    djdx=-yy/rhoradius_2;       
    djdy=xx/rhoradius_2;      

    if(make_conformal_derivs)
    {
      d2Rdxx=(1-xx*xx/R_2)/R;
      d2Rdxy=-xx*yy/R_3;
      d2Rdxz=-xx*zz/R_3;
      d2Rdyy=(1-yy*yy/R_2)/R;
      d2Rdyz=-zz*yy/R_3;
      d2Rdzz=(1-zz*zz/R_2)/R;
      tmp=2/R_2+1/rhoradius_2;
      d2qdxx=dqdx*(1/xx-xx*tmp);
      d2qdxy=-dqdx*yy*tmp;
      d2qdyy=dqdy*(1/yy-yy*tmp);
      d2qdzz=-2*dqdz*zz/R_2;
      tmp=2/R_2-1/zz/zz;
      d2qdxz=-zz*dqdx*tmp;
      d2qdyz=-zz*dqdy*tmp;
    

      /* conformal factor partial derivatives */
      /* note second derivatives are not tensors */
      psix[i] = dRdx*Phi_R + dqdx*Phi_q;
      psiy[i] = dRdy*Phi_R + dqdy*Phi_q;
      psiz[i] = dRdz*Phi_R + dqdz*Phi_q;

      if(*conformal_state > 2)
      {
        psixx[i] = dRdx*dRdx*Phi_RR + d2Rdxx*Phi_R + dqdx*dqdx*Phi_qq 
            + d2qdxx*Phi_q + 2*dRdx*dqdx*Phi_Rq;
        psixy[i] = dRdy*dRdx*Phi_RR + d2Rdxy*Phi_R + dqdy*dqdx*Phi_qq 
            + d2qdxy*Phi_q + (dRdx*dqdy+dRdy*dqdx)*Phi_Rq;
        psixz[i] = dRdz*dRdx*Phi_RR + d2Rdxz*Phi_R + dqdz*dqdx*Phi_qq 
            + d2qdxz*Phi_q + (dRdx*dqdz+dRdz*dqdx)*Phi_Rq;
        psiyy[i] = dRdy*dRdy*Phi_RR + d2Rdyy*Phi_R + dqdy*dqdy*Phi_qq 
            + d2qdyy*Phi_q + 2*dRdy*dqdy*Phi_Rq;
        psiyz[i] = dRdz*dRdy*Phi_RR + d2Rdyz*Phi_R + dqdz*dqdy*Phi_qq 
            + d2qdyz*Phi_q + (dRdy*dqdz+dRdz*dqdy)*Phi_Rq;
        psizz[i] = dRdz*dRdz*Phi_RR + d2Rdzz*Phi_R + dqdz*dqdz*Phi_qq 
            + d2qdzz*Phi_q + 2*dRdz*dqdz*Phi_Rq;
      }

    /*     Cactus convention
     *     -----------------
     */
      inv_psi = 1 / psi[i];

      psix[i]  *= inv_psi;
      psiy[i]  *= inv_psi;
      psiz[i]  *= inv_psi;

      if(*conformal_state > 2)
      {
        psixx[i] *= inv_psi;
        psixy[i] *= inv_psi;
        psixz[i] *= inv_psi;
        psiyy[i] *= inv_psi;
        psiyz[i] *= inv_psi;
        psizz[i] *= inv_psi;
      }
    }

    /* metric */
    tmp=(Chi2-1)*R_2*sth_2;
    gxx[i] = 1 + djdx*djdx*tmp;
    gxy[i] = djdx*djdy*tmp;
    gxz[i] = 0;
    gyy[i] = 1 + djdy*djdy*tmp;
    gyz[i] = 0;
    gzz[i] = 1;
  

    /* extrinsic curvature */
    kxx[i] = 2*(dRdx*KRj+dqdx*Kqj)*djdx;
    kxy[i] = (dRdx*KRj+dqdx*Kqj)*djdy + (dRdy*KRj+dqdy*Kqj)*djdx;
    kxz[i] = (dRdz*KRj+dqdz*Kqj)*djdx;
    kyy[i] = 2*(dRdy*KRj+dqdy*Kqj)*djdy;
    kyz[i] = (dRdz*KRj+dqdz*Kqj)*djdy;
    kzz[i] = 0;

    /*    probable convention to define conformal extrinsic curvature:
          if (*conformal_state == CONFORMAL_METRIC)
          {
          tmp=Psi*Psi;
          kxx[i] *=tmp;
          kxy[i] *=tmp;
          kxz[i] *=tmp;
          kyy[i] *=tmp;
          kyz[i] *=tmp;
      }*/

    /* shift */
    if(do_shift)
    {
      shiftx[i]=-yy*shift_phi;
      shifty[i]=xx*shift_phi;
      shiftz[i]=0;
    }
  }

  /*     Metric depends on conformal state
   *     ---------------------------------
   */
  
  if(CCTK_EQUALS(metric_type, "physical"))
  {
    for(i = 0; i < npoints; i++)
    {
      /* if((i/1000)*1000==i||i>35930)printf("i=%i\n",i); */
      tmp=psi[i];tmp*=tmp;tmp*=tmp;
      gxx[i] *= tmp;
      gxy[i] *= tmp;
      gyy[i] *= tmp;
      gzz[i] *= tmp;
      psi[i]  = 1;
    }
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, gxz);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, gyz);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psix);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiy);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiz);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psixx);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psixy);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psixz);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiyy);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psiyz);
    IDAnalyticBH_zero_CCTK_REAL_array(npoints, psizz);
  }
  
  return;
}

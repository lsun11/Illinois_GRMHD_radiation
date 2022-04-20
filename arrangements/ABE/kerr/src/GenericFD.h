/*@@                                         
   @file      GenericFD/src/GenericFD.h
   @date      June 16 2002
   @author    S. Husa                           
   @desc

   $Id: GenericFD.h,v 1.2 2004/04/05 15:23:50 shusa Exp $                                  
   
   @enddesc                                     
   @@*/                                           

/*  Copyright 2004 Sascha Husa, Ian Hinder, Christiane Lechner

    This file is part of Kranc.

    Kranc is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Kranc is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef NOPRECOMPUTE
#define PRECOMPUTE
#endif

#if defined(FD_C0) || defined(FD_C2) || defined(FD_C4) || defined(FD_C6) || defined(FD_C2C4)
#define FD_SET_BY_USER
#endif
#ifndef FD_SET_BY_USER
#define FD_C6
#endif


#if defined(FD_C0)
#define FD_METHOD_DESC  "FD method: replace derivatives by zero"
#endif

#if defined(FD_C2)
#define FD_METHOD_DESC  "FD method: second order finite differences"
#endif

#if defined(FD_C4)
#define FD_METHOD_DESC  "FD method: fourth order finite differences"
#endif

#if defined(FD_C6)
#define FD_METHOD_DESC  "FD method: sixth order finite differences"
#endif

#if defined(FD_C2C4)
#define FD_METHOD_DESC  "FD method: weighted 2nd/4th order centered finite differences"
#endif



/* utility functions  */
#if defined(KRANC_C)
#define string(d,f) d ## f
#else
#define string(d,f) d/**/f
#endif

/* finite differencing macros */

/*                                   */
/* add method argument to shorthands */
/*                                   */

/* second derivatives */
#define D11x(gf)  string(D11,gf)
#define D22x(gf)  string(D22,gf)
#define D33x(gf)  string(D33,gf)
#define D21x(gf)  string(D21,gf)
#define D32x(gf)  string(D32,gf)
#define D31x(gf)  string(D31,gf)

/* first derivatives                */
#define D1x(gf)   string(D1,gf)
#define D2x(gf)   string(D2,gf)
#define D3x(gf)   string(D3,gf)

#ifdef PRECOMPUTE

/* second derivatives */
#define D11(gf,i,j,k)  string(D11,gf)
#define D22(gf,i,j,k)  string(D22,gf)
#define D33(gf,i,j,k)  string(D33,gf)
#define D21(gf,i,j,k)  string(D21,gf)
#define D32(gf,i,j,k)  string(D32,gf)
#define D31(gf,i,j,k)  string(D31,gf) 

/* first derivatives                */
#define D1(gf,i,j,k)   string(D1,gf)
#define D2(gf,i,j,k)   string(D2,gf) 
#define D3(gf,i,j,k)   string(D3,gf) 

#else

/* second derivatives */
#define D11(gf,i,j,k)  D11gf(gf,i,j,k)
#define D22(gf,i,j,k)  D22gf(gf,i,j,k)
#define D33(gf,i,j,k)  D33gf(gf,i,j,k)
#define D21(gf,i,j,k)  D21gf(gf,i,j,k)
#define D32(gf,i,j,k)  D32gf(gf,i,j,k)
#define D31(gf,i,j,k)  D31gf(gf,i,j,k)

/* first derivatives                */
#define D1(gf,i,j,k)   D1gf(gf, i,j,k)
#define D2(gf,i,j,k)   D2gf(gf, i,j,k)
#define D3(gf,i,j,k)   D3gf(gf, i,j,k)

#endif

#ifdef FD_C0
/* second derivatives */
#define D11gf(gf,i,j,k)  D11_c0(gf,i,j,k)
#define D22gf(gf,i,j,k)  D22_c0(gf,i,j,k)
#define D33gf(gf,i,j,k)  D33_c0(gf,i,j,k)
#define D21gf(gf,i,j,k)  D21_c0(gf,i,j,k)
#define D32gf(gf,i,j,k)  D32_c0(gf,i,j,k)
#define D31gf(gf,i,j,k)  D31_c0(gf,i,j,k)

/* first derivatives                */
#define D1gf(gf,i,j,k)   D1_c0(gf, i,j,k)
#define D2gf(gf,i,j,k)   D2_c0(gf, i,j,k)
#define D3gf(gf,i,j,k)   D3_c0(gf, i,j,k)
#endif



#ifdef FD_C2
/* second derivatives */
#define D11gf(gf,i,j,k)  D11_c2(gf,i,j,k)
#define D22gf(gf,i,j,k)  D22_c2(gf,i,j,k)
#define D33gf(gf,i,j,k)  D33_c2(gf,i,j,k)
#define D21gf(gf,i,j,k)  D21_c2(gf,i,j,k)
#define D32gf(gf,i,j,k)  D32_c2(gf,i,j,k)
#define D31gf(gf,i,j,k)  D31_c2(gf,i,j,k)

/* first derivatives                */
#define D1gf(gf,i,j,k)   D1_c2(gf, i,j,k)
#define D2gf(gf,i,j,k)   D2_c2(gf, i,j,k)
#define D3gf(gf,i,j,k)   D3_c2(gf, i,j,k)

#define D1_up_gt_gf(gf,i,j,k)  D1_up_gt_c2(gf,i,j,k)
#define D1_up_lt_gf(gf,i,j,k)  D1_up_lt_c2(gf,i,j,k)
#define D2_up_gt_gf(gf,i,j,k)  D2_up_gt_c2(gf,i,j,k)
#define D2_up_lt_gf(gf,i,j,k)  D2_up_lt_c2(gf,i,j,k)
#define D3_up_gt_gf(gf,i,j,k)  D3_up_gt_c2(gf,i,j,k)
#define D3_up_lt_gf(gf,i,j,k)  D3_up_lt_c2(gf,i,j,k)
#endif


#ifdef FD_C4
/* second derivatives */
#define D11gf(gf,i,j,k)  D11_c4(gf,i,j,k) 
#define D22gf(gf,i,j,k)  D22_c4(gf,i,j,k)
#define D33gf(gf,i,j,k)  D33_c4(gf,i,j,k)
#define D21gf(gf,i,j,k)  D21_c4(gf,i,j,k) 
#define D32gf(gf,i,j,k)  D32_c4(gf,i,j,k) 
#define D31gf(gf,i,j,k)  D31_c4(gf,i,j,k)

/* 4th order upwind differencing, as defined in Campanelli et al. gr-qc/0505055 (See below!) */

#define D1_up_gt_gf(gf,i,j,k)  D1_up_gt_c4(gf,i,j,k)
#define D1_up_lt_gf(gf,i,j,k)  D1_up_lt_c4(gf,i,j,k)
#define D2_up_gt_gf(gf,i,j,k)  D2_up_gt_c4(gf,i,j,k)
#define D2_up_lt_gf(gf,i,j,k)  D2_up_lt_c4(gf,i,j,k)
#define D3_up_gt_gf(gf,i,j,k)  D3_up_gt_c4(gf,i,j,k)
#define D3_up_lt_gf(gf,i,j,k)  D3_up_lt_c4(gf,i,j,k)

/* first derivatives                */
#define D1gf(gf,i,j,k)   D1_c4(gf, i,j,k) 
#define D2gf(gf,i,j,k)   D2_c4(gf, i,j,k)
#define D3gf(gf,i,j,k)   D3_c4(gf, i,j,k)
#endif


#ifdef FD_C6
/* second derivatives */
#define D11gf(gf,i,j,k)  D11_c6(gf,i,j,k) 
#define D22gf(gf,i,j,k)  D22_c6(gf,i,j,k)
#define D33gf(gf,i,j,k)  D33_c6(gf,i,j,k)
#define D21gf(gf,i,j,k)  D21_c6(gf,i,j,k) 
#define D32gf(gf,i,j,k)  D32_c6(gf,i,j,k) 
#define D31gf(gf,i,j,k)  D31_c6(gf,i,j,k)

/* 4th order upwind differencing, as defined in Campanelli et al. gr-qc/0505055 (See below!) */

#define D1_up_gt_gf(gf,i,j,k)  D1_up_gt_c6(gf,i,j,k)
#define D1_up_lt_gf(gf,i,j,k)  D1_up_lt_c6(gf,i,j,k)
#define D2_up_gt_gf(gf,i,j,k)  D2_up_gt_c6(gf,i,j,k)
#define D2_up_lt_gf(gf,i,j,k)  D2_up_lt_c6(gf,i,j,k)
#define D3_up_gt_gf(gf,i,j,k)  D3_up_gt_c6(gf,i,j,k)
#define D3_up_lt_gf(gf,i,j,k)  D3_up_lt_c6(gf,i,j,k)

/* first derivatives                */
#define D1gf(gf,i,j,k)   D1_c6(gf, i,j,k) 
#define D2gf(gf,i,j,k)   D2_c6(gf, i,j,k)
#define D3gf(gf,i,j,k)   D3_c6(gf, i,j,k)
#endif


#ifdef FD_C2C4
/* second derivatives */
#define D11gf(gf,i,j,k)  D11_c2c4(gf,i,j,k)
#define D22gf(gf,i,j,k)  D22_c2c4(gf,i,j,k)
#define D33gf(gf,i,j,k)  D33_c2c4(gf,i,j,k)
#define D21gf(gf,i,j,k)  D21_c2c4(gf,i,j,k)
#define D32gf(gf,i,j,k)  D32_c2c4(gf,i,j,k)
#define D31gf(gf,i,j,k)  D31_c2c4(gf,i,j,k)

/* first derivatives                */
#define D1gf(gf,i,j,k)   D1_c2c4(gf, i,j,k)
#define D2gf(gf,i,j,k)   D2_c2c4(gf, i,j,k)
#define D3gf(gf,i,j,k)   D3_c2c4(gf, i,j,k)
#endif




/*****************************************************/
/*                                                   */
/*             METHODS                               */
/*                                                   */
/*****************************************************/

/*  c0                              */

/*  set all derivatives = 0         */
/*  for debugging and benchmarking  */

/*  second derivatives              */

#define D11_c0(gf,i,j,k)   0.
#define D22_c0(gf,i,j,k)   0.
#define D33_c0(gf,i,j,k)   0.
#define D21_c0(gf,i,j,k)   0.
#define D32_c0(gf,i,j,k)   0.
#define D31_c0(gf,i,j,k)   0.

/* first derivatives                */

#define D1_c0(gf,i,j,k)    0.
#define D2_c0(gf,i,j,k)    0.
#define D3_c0(gf,i,j,k)    0.


/*  c2                  */
/*                      */
/*  2nd order centered  */
/*                      */
/* second derivatives, centered, 2nd order */

#ifndef KRANC_C
#define D11_c2(gf,i,j,k)                         \
	 ((   gf(i+1,j,k) \
	 - 2.*gf(i,  j,k) \
	 +    gf(i-1,j,k)) * dxi2 )

#define D22_c2(gf,i,j,k)                         \
	 ((   gf(i,j+1,k) \
	 - 2.*gf(i,j,  k) \
	 +    gf(i,j-1,k)) * dyi2 )

#define D33_c2(gf,i,j,k)                         \
	 ((   gf(i,j,k+1) \
	 - 2.*gf(i,j,k  ) \
	 +    gf(i,j,k-1)) * dzi2 )

#define D21_c2(gf,i,j,k)                        \
	 ((gf(i+1,j+1,k) \
	 + gf(i-1,j-1,k) \
	 - gf(i+1,j-1,k) \
	 - gf(i-1,j+1,k)) * hdxi * hdyi )

#define D32_c2(gf,i,j,k)                        \
	 ((gf(i,j+1,k+1) \
	 + gf(i,j-1,k-1) \
	 - gf(i,j+1,k-1) \
	 - gf(i,j-1,k+1)) * hdyi * hdzi )

#define D31_c2(gf,i,j,k)                        \
	 ((gf(i+1,j,k+1) \
	 + gf(i-1,j,k-1) \
	 - gf(i+1,j,k-1) \
	 - gf(i-1,j,k+1)) * hdxi * hdzi )

/* first derivatives, centered, 2nd order */

#define D1_c2(gf,i,j,k)                       \
	 ((gf(i+1,j,k) \
	 - gf(i-1,j,k)) * hdxi)

#define D2_c2(gf,i,j,k)                       \
	 ((gf(i,j+1,k) \
	 - gf(i,j-1,k)) * hdyi)

#define D3_c2(gf,i,j,k)                       \
	 ((gf(i,j,k+1) \
	 - gf(i,j,k-1)) * hdzi)
	
#else
#define D11_c2(gf,i,j,k)                         \
	 ((   gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
	 - 2.*gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
	 +    gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]) * dxi2 )

#define D22_c2(gf,i,j,k)                         \
	 ((   gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
	 - 2.*gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
	 +    gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]) * dyi2 )

#define D33_c2(gf,i,j,k)                         \
	 ((   gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
	 - 2.*gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
	 +    gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]) * dzi2 )

#define D21_c2(gf,i,j,k)                        \
	 ((gf[CCTK_GFINDEX3D(cctkGH,i+1,j+1,k)] \
	 + gf[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i+1,j-1,k)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i-1,j+1,k)]) * hdxi * hdyi )

#define D32_c2(gf,i,j,k)                        \
	 ((gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k+1)] \
	 + gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k-1)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k+1)]) * hdyi * hdzi )

#define D31_c2(gf,i,j,k)                        \
	 ((gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k+1)] \
	 + gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k-1)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k+1)]) * hdxi * hdzi )

/* first derivatives, centered, 2nd order */

#define D1_c2(gf,i,j,k)                       \
	 ((gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]) * hdxi)

#define D2_c2(gf,i,j,k)                       \
	 ((gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]) * hdyi)

#define D3_c2(gf,i,j,k)                       \
	 ((gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
	 - gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]) * hdzi)

//Standard second order upwind differencing!
#define D1_up_gt_c2(gf,i,j,k)                 \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
         + 4. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
	 - 3. * gf[CCTK_GFINDEX3D(cctkGH,i  ,j,k)]) * (hdxi))

#define D1_up_lt_c2(gf,i,j,k)                 \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)] \
         - 4. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
	 + 3. * gf[CCTK_GFINDEX3D(cctkGH,i  ,j,k)]) * (hdxi))

#define D2_up_gt_c2(gf,i,j,k)                 \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
         + 4. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
	 - 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j  ,k)]) * (hdxi))

#define D2_up_lt_c2(gf,i,j,k)                 \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)] \
         - 4. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
	 + 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j  ,k)]) * (hdxi))

#define D3_up_gt_c2(gf,i,j,k)                 \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
         + 4. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
	 - 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k  )]) * (hdxi))

#define D3_up_lt_c2(gf,i,j,k)                 \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)] \
         - 4. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
	 + 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k  )]) * (hdxi))

#endif	

	
/*  c4                  */
/*                      */
/*  4th order centered  */
/*                      */

/* second derivatives, centered, 4th order */
#ifndef KRANC_C
#define D11_c4(gf,i,j,k)                    \
	 ((-       gf(i+2,j,k) \
	   + 16. * gf(i+1,j,k) \
	   - 30. * gf(i ,j,k) \
	   + 16. * gf(i-1,j,k) \
	   -       gf(i-2,j,k)) * dxi2 / 12.)


#define D22_c4(gf,i,j,k)                    \
	 ((-       gf(i,j+2,k) \
	   + 16. * gf(i,j+1,k) \
	   - 30. * gf(i ,j,k) \
	   + 16. * gf(i,j-1,k) \
	   -       gf(i,j-2,k)) * dyi2 / 12.)


#define D33_c4(gf,i,j,k)                              \
	 ((-       gf(i,j,k+2) \
	   + 16. * gf(i,j,k+1) \
	   - 30. * gf(i ,j,k) \
	   + 16. * gf(i,j,k-1) \
	   -       gf(i,j,k-2)) * dzi2 / 12.)

#define D21_c4(gf,i,j,k) (                                  \
	 (+1.*(        gf(i-2,j-2,k) \
               -  8. * gf(i-1,j-2,k) \
               +  8. * gf(i+1,j-2,k) \
               -       gf(i+2,j-2,k)) \
          -8.*(        gf(i-2,j-1,k) \
               -  8. * gf(i-1,j-1,k) \
               +  8. * gf(i+1,j-1,k) \
               -       gf(i+2,j-1,k))\
          +8.*(        gf(i-2,j+1,k) \
               -  8. * gf(i-1,j+1,k) \
               +  8. * gf(i+1,j+1,k) \
               -       gf(i+2,j+1,k))\
          -1.*(        gf(i-2,j+2,k) \
               -  8. * gf(i-1,j+2,k) \
               +  8. * gf(i+1,j+2,k) \
               -       gf(i+2,j+2,k))) * dxi * dyi / 144.)

#define D31_c4(gf,i,j,k) (                                  \
	 (+1.*(        gf(i-2,j,k-2) \
               -  8. * gf(i-1,j,k-2) \
               +  8. * gf(i+1,j,k-2) \
               -       gf(i+2,j,k-2))\
          -8.*(        gf(i-2,j,k-1) \
               -  8. * gf(i-1,j,k-1) \
               +  8. * gf(i+1,j,k-1) \
               -       gf(i+2,j,k-1))\
          +8.*(        gf(i-2,j,k+1) \
               -  8. * gf(i-1,j,k+1) \
               +  8. * gf(i+1,j,k+1) \
               -       gf(i+2,j,k+1))\
          -1.*(        gf(i-2,j,k+2) \
               -  8. * gf(i-1,j,k+2) \
               +  8. * gf(i+1,j,k+2) \
               -       gf(i+2,j,k+2))) * dxi * dzi / 144.)

#define D32_c4(gf,i,j,k) (                                  \
	 (+1.*(        gf(i,j-2,k-2) \
               -  8. * gf(i,j-1,k-2) \
               +  8. * gf(i,j+1,k-2) \
               -       gf(i,j+2,k-2))\
          -8.*(        gf(i,j-2,k-1) \
               -  8. * gf(i,j-1,k-1) \
               +  8. * gf(i,j+1,k-1) \
               -       gf(i,j+2,k-1))\
          +8.*(        gf(i,j-2,k+1) \
               -  8. * gf(i,j-1,k+1) \
               +  8. * gf(i,j+1,k+1) \
               -       gf(i,j+2,k+1))\
          -1.*(        gf(i,j-2,k+2) \
               -  8. * gf(i,j-1,k+2) \
               +  8. * gf(i,j+1,k+2) \
               -       gf(i,j+2,k+2))) * dyi * dzi / 144.)

/*
#define D21_c4(gf,i,j,k)                                    \
	     ((-       gf(i+2,j+2,k) \
               +       gf(i+2,j-2,k) \
	       +       gf(i-2,j+2,k) \
               -       gf(i-2,j-2,k) \
               + 16. * gf(i+1,j+1,k) \
    	       - 16. * gf(i+1,j-1,k) \
	       - 16. * gf(i-1,j+1,k) \
	       + 16. * gf(i-1,j-1,k)) * dxi * dyi / 48.)

#define D31_c4(gf,i,j,k)                                    \
             ((-       gf(i+2,j,k+2) \
               +       gf(i+2,j,k-2) \
               +       gf(i-2,j,k+2) \
               -       gf(i-2,j,k-2) \
               + 16. * gf(i+1,j,k+1) \
               - 16. * gf(i+1,j,k-1) \
               - 16. * gf(i-1,j,k+1) \
               + 16. * gf(i-1,j,k-1)) * dxi * dzi / 48.)

	
#define D32_c4(gf,i,j,k)                                    \
             ((-       gf(i,j+2,k+2) \
               +       gf(i,j+2,k-2) \
               +       gf(i,j-2,k+2) \
               -       gf(i,j-2,k-2) \
               + 16. * gf(i,j+1,k+1) \
               - 16. * gf(i,j+1,k-1) \
               - 16. * gf(i,j-1,k+1) \
               + 16. * gf(i,j-1,k-1)) * dzi * dyi / 48.)
*/

#define D1_c4(gf,i,j,k)                            \
       ((-      gf(i+2,j,k) \
         + 8. * gf(i+1,j,k) \
         - 8. * gf(i-1,j,k) \
	 +      gf(i-2,j,k)) * (dxi / 12.))

#define D2_c4(gf,i,j,k)                            \
       ((-      gf(i,j+2,k) \
         + 8. * gf(i,j+1,k) \
         - 8. * gf(i,j-1,k) \
	 +      gf(i,j-2,k)) * (dyi / 12.))

#define D3_c4(gf,i,j,k)                            \
       ((-      gf(i,j,k+2) \
         + 8. * gf(i,j,k+1) \
         - 8. * gf(i,j,k-1) \
	 +      gf(i,j,k-2)) * (dzi / 12.))

#else
#define D11_c4(gf,i,j,k)                    \
	 ((-       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
	   + 16. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
	   - 30. * gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
	   + 16. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
	   -       gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]) * dxi2 / 12.)


#define D22_c4(gf,i,j,k)                    \
	 ((-       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
	   + 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
	   - 30. * gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
	   + 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
	   -       gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]) * dyi2 / 12.)


#define D33_c4(gf,i,j,k)                              \
	 ((-       gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
	   + 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
	   - 30. * gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
	   + 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
	   -       gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]) * dzi2 / 12.)

/*
#define D21_c4(gf,i,j,k)                                    \
	     ((-       gf[CCTK_GFINDEX3D(cctkGH,i+2,j+2,k)] \
               +       gf[CCTK_GFINDEX3D(cctkGH,i+2,j-2,k)] \
	       +       gf[CCTK_GFINDEX3D(cctkGH,i-2,j+2,k)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i-2,j-2,k)] \
               + 16. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j+1,k)] \
    	       - 16. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j-1,k)] \
	       - 16. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j+1,k)] \
	       + 16. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)]) * dxi * dyi / 48.)

#define D31_c4(gf,i,j,k)                                    \
             ((-       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k+2)] \
               +       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k-2)] \
               +       gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k+2)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k-2)] \
               + 16. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k+1)] \
               - 16. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k-1)] \
               - 16. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k+1)] \
               + 16. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)]) * dxi * dzi / 48.)

	
#define D32_c4(gf,i,j,k)                                    \
             ((-       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k+2)] \
               +       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k-2)] \
               +       gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k+2)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k-2)] \
               + 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k+1)] \
               - 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k-1)] \
               - 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k+1)] \
               + 16. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)]) * dzi * dyi / 48.)
*/

/* AS IN CAMP ET AL: */

#define D21_c4(gf,i,j,k) (                                  \
	 (+1.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j-2,k)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j-2,k)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j-2,k)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j-2,k)])\
          -8.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j-1,k)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j-1,k)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j-1,k)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j-1,k)])\
          +8.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j+1,k)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j+1,k)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j+1,k)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j+1,k)])\
          -1.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j+2,k)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j+2,k)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j+2,k)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j+2,k)])) * dxi * dyi / 144.)

#define D31_c4(gf,i,j,k) (                                  \
	 (+1.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k-2)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k-2)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k-2)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k-2)])\
          -8.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k-1)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k-1)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k-1)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k-1)])\
          +8.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k+1)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k+1)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k+1)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k+1)])\
          -1.*(        gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k+2)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k+2)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k+2)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k+2)])) * dxi * dzi / 144.)

#define D32_c4(gf,i,j,k) (                                  \
	 (+1.*(        gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k-2)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k-2)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k-2)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k-2)])\
          -8.*(        gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k-1)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k-1)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k-1)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k-1)])\
          +8.*(        gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k+1)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k+1)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k+1)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k+1)])\
          -1.*(        gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k+2)] \
               -  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k+2)] \
               +  8. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k+2)] \
               -       gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k+2)])) * dyi * dzi / 144.)


/* 4th order upwind differencing, as defined in Campanelli et al. gr-qc/0505055 */

#define D1_up_gt_c4(gf,i,j,k)                      \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i+3,j,k)] \
         - 6. * gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
         +18. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
         -10. * gf[CCTK_GFINDEX3D(cctkGH,i  ,j,k)] \
	 - 3. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]) * (dxi / 12.))

#define D1_up_lt_c4(gf,i,j,k)                      \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i-3,j,k)] \
         + 6. * gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)] \
         -18. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
         +10. * gf[CCTK_GFINDEX3D(cctkGH,i  ,j,k)] \
	 + 3. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]) * (dxi / 12.))

#define D2_up_gt_c4(gf,i,j,k)                      \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j+3,k)] \
         - 6. * gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
         +18. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
         -10. * gf[CCTK_GFINDEX3D(cctkGH,i,j  ,k)] \
	 - 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]) * (dyi / 12.))

#define D2_up_lt_c4(gf,i,j,k)                      \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j-3,k)] \
         + 6. * gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)] \
         -18. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
         +10. * gf[CCTK_GFINDEX3D(cctkGH,i,j  ,k)] \
	 + 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]) * (dyi / 12.))

#define D3_up_gt_c4(gf,i,j,k)                      \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j,k+3)] \
         - 6. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
         +18. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
         -10. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k  )] \
	 - 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]) * (dzi / 12.))

#define D3_up_lt_c4(gf,i,j,k)                      \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j,k-3)] \
         + 6. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)] \
         -18. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
         +10. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k  )] \
	 + 3. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]) * (dzi / 12.))

/* first derivatives, centered, 4th order */

#define D1_c4(gf,i,j,k)                            \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
         + 8. * gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
         - 8. * gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
	 +      gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]) * (dxi / 12.))

#define D2_c4(gf,i,j,k)                            \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
         + 8. * gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
         - 8. * gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
	 +      gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]) * (dyi / 12.))

#define D3_c4(gf,i,j,k)                            \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
         + 8. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
         - 8. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
	 +      gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]) * (dxi / 12.))
#endif
/* END OF FOURTH ORDER */

/* SIXTH ORDER FD STENCIL DEFINITIONS */

#ifndef KRANC_C
/* Sixth order second derivatives a la Husa et al arxiv.org 0706.0740. 
   Note that we use more efficient stencils for cross derivatives! */

#define D11_c6(gf,i,j,k)                    \
         ((+   2.* gf(i+3,j,k) \
           -  27.* gf(i+2,j,k) \
           + 270.* gf(i+1,j,k) \
           - 490.* gf(i ,j,k) \
           + 270.* gf(i-1,j,k) \
           -  27.* gf(i-2,j,k) \
           +   2.* gf(i-3,j,k)) * dxi2 / 180.)

#define D22_c6(gf,i,j,k)                    \
         ((+   2.* gf(i,j+3,k) \
           -  27.* gf(i,j+2,k) \
           + 270.* gf(i,j+1,k) \
           - 490.* gf(i,j ,k) \
           + 270.* gf(i,j-1,k) \
           -  27.* gf(i,j-2,k) \
           +   2.* gf(i,j-3,k)) * dyi2 / 180.)

#define D33_c6(gf,i,j,k)                    \
         ((+   2.* gf(i,j,k+3) \
           -  27.* gf(i,j,k+2) \
           + 270.* gf(i,j,k+1) \
           - 490.* gf(i,j,k ) \
           + 270.* gf(i,j,k-1) \
           -  27.* gf(i,j,k-2) \
           +   2.* gf(i,j,k-3)) * dzi2 / 180.)

#define D21_c6(gf,i,j,k)                    \
        ((+     D1_c6(gf,i,j+3,k) \
         - 9. * D1_c6(gf,i,j+2,k) \
         + 45.* D1_c6(gf,i,j+1,k) \
         - 45.* D1_c6(gf,i,j-1,k) \
         + 9. * D1_c6(gf,i,j-2,k) \
         -      D1_c6(gf,i,j-3,k)) * (dyi / 60.))

#define D32_c6(gf,i,j,k)                    \
        ((+     D2_c6(gf,i,j,k+3) \
         - 9. * D2_c6(gf,i,j,k+2) \
         + 45.* D2_c6(gf,i,j,k+1) \
         - 45.* D2_c6(gf,i,j,k-1) \
         + 9. * D2_c6(gf,i,j,k-2) \
         -      D2_c6(gf,i,j,k-3)) * (dzi / 60.))

#define D31_c6(gf,i,j,k)                    \
        ((+     D1_c6(gf,i,j,k+3) \
         - 9. * D1_c6(gf,i,j,k+2) \
         + 45.* D1_c6(gf,i,j,k+1) \
         - 45.* D1_c6(gf,i,j,k-1) \
         + 9. * D1_c6(gf,i,j,k-2) \
         -      D1_c6(gf,i,j,k-3)) * (dzi / 60.))

/* first derivatives, centered, 6th order */

#define D1_c6(gf,i,j,k)                            \
       ((+      gf(i+3,j,k) \
         - 9. * gf(i+2,j,k) \
         + 45.* gf(i+1,j,k) \
         - 45.* gf(i-1,j,k) \
         + 9. * gf(i-2,j,k) \
         -      gf(i-3,j,k)) * (dxi / 60.))

#define D2_c6(gf,i,j,k)                            \
       ((+      gf(i,j+3,k) \
         - 9. * gf(i,j+2,k) \
         + 45.* gf(i,j+1,k) \
         - 45.* gf(i,j-1,k) \
         + 9. * gf(i,j-2,k) \
         -      gf(i,j-3,k)) * (dyi / 60.))

#define D3_c6(gf,i,j,k)                            \
       ((+      gf(i,j,k+3) \
         - 9. * gf(i,j,k+2) \
         + 45.* gf(i,j,k+1) \
         - 45.* gf(i,j,k-1) \
         + 9. * gf(i,j,k-2) \
         -      gf(i,j,k-3)) * (dzi / 60.))

#else

/* Sixth order second derivatives a la Husa et al arxiv.org 0706.0740. 
   Note that we use more efficient stencils for cross derivatives! */
/* first derivatives, centered, 6th order */

#define D11_c6(gf,i,j,k)                    \
         ((+   2.* gf[CCTK_GFINDEX3D(cctkGH,i+3,j,k)] \
           -  27.* gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
           + 270.* gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
           - 490.* gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] \
           + 270.* gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
           -  27.* gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)] \
           +   2.* gf[CCTK_GFINDEX3D(cctkGH,i-3,j,k)]) * dxi2 / 180.)

#define D22_c6(gf,i,j,k)                    \
         ((+   2.* gf[CCTK_GFINDEX3D(cctkGH,i,j+3,k)] \
           -  27.* gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
           + 270.* gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
           - 490.* gf[CCTK_GFINDEX3D(cctkGH,i,j ,k)] \
           + 270.* gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
           -  27.* gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)] \
           +   2.* gf[CCTK_GFINDEX3D(cctkGH,i,j-3,k)]) * dyi2 / 180.)

#define D33_c6(gf,i,j,k)                    \
         ((+   2.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+3)] \
           -  27.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
           + 270.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
           - 490.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k )] \
           + 270.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
           -  27.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)] \
           +   2.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-3)]) * dzi2 / 180.)

#define D21_c6(gf,i,j,k)                    \
        ((+     D1_c6(gf,i,j+3,k) \
         - 9. * D1_c6(gf,i,j+2,k) \
         + 45.* D1_c6(gf,i,j+1,k) \
         - 45.* D1_c6(gf,i,j-1,k) \
         + 9. * D1_c6(gf,i,j-2,k) \
         -      D1_c6(gf,i,j-3,k)) * (dyi / 60.))

#define D32_c6(gf,i,j,k)                    \
        ((+     D2_c6(gf,i,j,k+3) \
         - 9. * D2_c6(gf,i,j,k+2) \
         + 45.* D2_c6(gf,i,j,k+1) \
         - 45.* D2_c6(gf,i,j,k-1) \
         + 9. * D2_c6(gf,i,j,k-2) \
         -      D2_c6(gf,i,j,k-3)) * (dzi / 60.))

#define D31_c6(gf,i,j,k)                    \
        ((+     D1_c6(gf,i,j,k+3) \
         - 9. * D1_c6(gf,i,j,k+2) \
         + 45.* D1_c6(gf,i,j,k+1) \
         - 45.* D1_c6(gf,i,j,k-1) \
         + 9. * D1_c6(gf,i,j,k-2) \
         -      D1_c6(gf,i,j,k-3)) * (dzi / 60.))

#define D1_c6(gf,i,j,k)                            \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i+3,j,k)] \
         - 9. * gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
         + 45.* gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
         - 45.* gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
         + 9. * gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)] \
         -      gf[CCTK_GFINDEX3D(cctkGH,i-3,j,k)]) * (dxi / 60.))

#define D2_c6(gf,i,j,k)                            \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j+3,k)] \
         - 9. * gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
         + 45.* gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
         - 45.* gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
         + 9. * gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)] \
         -      gf[CCTK_GFINDEX3D(cctkGH,i,j-3,k)]) * (dyi / 60.))

#define D3_c6(gf,i,j,k)                            \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j,k+3)] \
         - 9. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
         + 45.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
         - 45.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
         + 9. * gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)] \
         -      gf[CCTK_GFINDEX3D(cctkGH,i,j,k-3)]) * (dzi / 60.))

/* 6th order upwind differencing, as defined in Husa et al arxiv.org 0706.0740 */

#define D1_up_gt_c6(gf,i,j,k)                      \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i+4,j,k)] \
         +  8.* gf[CCTK_GFINDEX3D(cctkGH,i+3,j,k)] \
         - 30.* gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)] \
         + 80.* gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
         - 35.* gf[CCTK_GFINDEX3D(cctkGH,i  ,j,k)] \
         - 24.* gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
         +  2.* gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)]) * (dxi / 60.))

#define D1_up_lt_c6(gf,i,j,k)                      \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i-4,j,k)] \
         -  8.* gf[CCTK_GFINDEX3D(cctkGH,i-3,j,k)] \
         + 30.* gf[CCTK_GFINDEX3D(cctkGH,i-2,j,k)] \
         - 80.* gf[CCTK_GFINDEX3D(cctkGH,i-1,j,k)] \
         + 35.* gf[CCTK_GFINDEX3D(cctkGH,i  ,j,k)] \
         + 24.* gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
         -  2.* gf[CCTK_GFINDEX3D(cctkGH,i+2,j,k)]) * (dxi / 60.))

#define D2_up_gt_c6(gf,i,j,k)                      \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j+4,k)] \
         +  8.* gf[CCTK_GFINDEX3D(cctkGH,i,j+3,k)] \
         - 30.* gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)] \
         + 80.* gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
         - 35.* gf[CCTK_GFINDEX3D(cctkGH,i,j  ,k)] \
         - 24.* gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
         +  2.* gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)]) * (dyi / 60.))

#define D2_up_lt_c6(gf,i,j,k)                      \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j-4,k)] \
         -  8.* gf[CCTK_GFINDEX3D(cctkGH,i,j-3,k)] \
         + 30.* gf[CCTK_GFINDEX3D(cctkGH,i,j-2,k)] \
         - 80.* gf[CCTK_GFINDEX3D(cctkGH,i,j-1,k)] \
         + 35.* gf[CCTK_GFINDEX3D(cctkGH,i,j  ,k)] \
         + 24.* gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
         -  2.* gf[CCTK_GFINDEX3D(cctkGH,i,j+2,k)]) * (dyi / 60.))

#define D3_up_gt_c6(gf,i,j,k)                      \
       ((-      gf[CCTK_GFINDEX3D(cctkGH,i,j,k+4)] \
         +  8.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+3)] \
         - 30.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)] \
         + 80.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
         - 35.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k  )] \
         - 24.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
         +  2.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)]) * (dzi / 60.))

#define D3_up_lt_c6(gf,i,j,k)                      \
       ((+      gf[CCTK_GFINDEX3D(cctkGH,i,j,k-4)] \
         -  8.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-3)] \
         + 30.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-2)] \
         - 80.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k-1)] \
         + 35.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k  )] \
         + 24.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
         -  2.* gf[CCTK_GFINDEX3D(cctkGH,i,j,k+2)]) * (dzi / 60.))

#endif

/*****************************************************/
/*                                                    */
/*         DERIVED METHODS                            */
/*                                                    */
/******************************************************/


/* blend c2 and c4 */
/* second derivatives */
#define D11_c2c4(gf,i,j,k) (fdweight_c2*D11_c2(gf,i,j,k) + fdweight_c4*D11_c4(gf,i,j,k))
#define D22_c2c4(gf,i,j,k) (fdweight_c2*D22_c2(gf,i,j,k) + fdweight_c4*D22_c4(gf,i,j,k))
#define D33_c2c4(gf,i,j,k) (fdweight_c2*D33_c2(gf,i,j,k) + fdweight_c4*D33_c4(gf,i,j,k))
#define D21_c2c4(gf,i,j,k) (fdweight_c2*D21_c2(gf,i,j,k) + fdweight_c4*D21_c4(gf,i,j,k))
#define D32_c2c4(gf,i,j,k) (fdweight_c2*D32_c2(gf,i,j,k) + fdweight_c4*D32_c4(gf,i,j,k))
#define D31_c2c4(gf,i,j,k) (fdweight_c2*D31_c2(gf,i,j,k) + fdweight_c4*D31_c4(gf,i,j,k))

/* first derivatives */
#define D1_c2c4(gf,i,j,k)  (fdweight_c2*D1_c2(gf, i,j,k) + fdweight_c4*D1_c4(gf,i,j,k))
#define D2_c2c4(gf,i,j,k)  (fdweight_c2*D2_c2(gf, i,j,k) + fdweight_c4*D2_c4(gf,i,j,k))
#define D3_c2c4(gf,i,j,k)  (fdweight_c2*D3_c2(gf, i,j,k) + fdweight_c4*D3_c4(gf,i,j,k))
	

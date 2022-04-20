#include "cctk.h"
#include "include/k2g_decl.h"

#define NAPF_HELP whichvec,det,uxx,uxy,uxz,uyy,uyz,uzz,tracek
#define DECLARE_NAPF_HELP CCTK_REAL det,uxx,uxy,uxz,uyy,uyz,uzz,tracek &&\
        integer whichvec

#define RICCI_POINT ricxx_loc,ricxy_loc,ricxz_loc,ricyy_loc,ricyz_loc,riczz_loc,ricscal_loc
#define DECLARE_RICCI_POINT CCTK_REAL ricxx_loc,ricxy_loc,ricxz_loc,ricyy_loc,ricyz_loc,riczz_loc,ricscal_loc

#define PFVEC v1xL,v2xL,v3xL,v1yL,v2yL,v3yL,v1zL,v2zL,v3zL
#define DECLARE_PFVEC \
      CCTK_REAL v1xL,v2xL,v3xL &&\
      CCTK_REAL v1yL,v2yL,v3yL &&\
      CCTK_REAL v1zL,v2zL,v3zL 

#define ELECMAG elecxx,elecxy,elecxz,elecyy,elecyz,eleczz, \
      magxx,magxy,magxz,magyy,magyz,magzz
#define DECLARE_ELECMAG  \
      CCTK_REAL elecxx,elecxy,elecxz,elecyy,elecyz,eleczz &&\
      CCTK_REAL magxx,magxy,magxz,magyy,magyz,magzz

#define NEWPEN psi0re_loc,psi0im_loc,psi1re_loc,psi1im_loc,psi2re_loc,\
      psi2im_loc,psi3re_loc,psi3im_loc,psi4re_loc,psi4im_loc
#define DECLARE_NEWPEN \
      CCTK_REAL psi0re_loc,psi0im_loc,psi1re_loc,psi1im_loc,psi2re_loc &&\
      CCTK_REAL psi2im_loc,psi3re_loc,psi3im_loc,psi4re_loc,psi4im_loc

#define RIEMIJ icurvre_loc, icurvim_loc, jcurvre_loc, jcurvim_loc
#define DECLARE_RIEMIJ CCTK_REAL RIEMIJ

#define LOCALKGx kkxx,kkxy,kkxz,kkyy,kkyz,kkzz,\
                 ggxx,ggxy,ggxz,ggyy,ggyz,ggzz,\
                 xx,yy,zz
#define DECLARE_LOCALKGx CCTK_REAL LOCALKGx


#define NAPSIFRIENDS NAPF_HELP, RICCI_POINT, PFVEC, ELECMAG, NEWPEN, RIEMIJ, LOCALKGx
#define DECLARE_NAPSIFRIENDS \
        DECLARE_NAPF_HELP &&\
        DECLARE_RICCI_POINT &&\
        DECLARE_PFVEC &&\
        DECLARE_ELECMAG &&\
        DECLARE_NEWPEN &&\
        DECLARE_RIEMIJ &&\
        DECLARE_LOCALKGx


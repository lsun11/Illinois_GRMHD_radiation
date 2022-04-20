/*This code handles prolongation of Az (z-comp of vec. potential),
  where array element (i,j,k) corresponds to points on the 
  semi-staggered grid [i+1/2,j+1/2,k].
  FOR DETAILED COMMENTS ON WHAT THIS CODE DOES, SEE
  prolongate_3d_stagger011.cc , whose code is just a 
  permutation of what this code does.
*/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "operator_prototypes.hh"
#include "typeprops.hh"

//#include "weno.h"

using namespace std;

double weno_5pt(double *f, double x, double xj, double h);

namespace CarpetLib {
  

  
#define SRCIND3(i,j,k)                          \
  index3 (i, j, k,                              \
          srciext, srcjext, srckext)
#define DSTIND3(i,j,k)                          \
  index3 (i, j, k,                              \
          dstiext, dstjext, dstkext)
  
  
  
  template <typename T>
  void
  prolongate_3d_stagger110 (T const * restrict const src,
                        ivect3 const & restrict srcext,
                        T * restrict const dst,
                        ivect3 const & restrict dstext,
                        ibbox3 const & restrict srcbbox,
                        ibbox3 const & restrict dstbbox,
                        ibbox3 const & restrict regbbox)
  {
    typedef typename typeprops<T>::real RT;
    
    if (any (srcbbox.stride() <= regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: strides disagree");
    }
    
    if (any (srcbbox.stride() != reffact2 * dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: source strides are not twice the destination strides");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
    
    
    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all ((regbbox.lower() - srcbbox.lower()) % regbbox.stride() == 0));
    ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / regbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % regbbox.stride() == 0));
    ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / regbbox.stride();
    
    
    
    bvect3 const needoffsetlo = srcoff % reffact2 != 0 or regext > 1;
    bvect3 const needoffsethi = (srcoff + regext - 1) % reffact2 != 0 or regext > 1;
    ivect3 const offsetlo = either (needoffsetlo, 3, 0);
    ivect3 const offsethi = either (needoffsethi, 3, 0);
    
    
    
    if (not regbbox.expand(offsetlo, offsethi).is_contained_in(srcbbox) or
        not regbbox                           .is_contained_in(dstbbox))
    {
      CCTK_WARN (0, "Internal error: region extent is not contained in array extent");
    }
    
    if (any (srcext != srcbbox.shape() / srcbbox.stride() or
             dstext != dstbbox.shape() / dstbbox.stride()))
    {
      CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes");
    }
    
    
    
    size_t const srciext = srcext[0];
    size_t const srcjext = srcext[1];
    size_t const srckext = srcext[2];
    
    size_t const dstiext = dstext[0];
    size_t const dstjext = dstext[1];
    size_t const dstkext = dstext[2];
    
    size_t const regiext = regext[0];
    size_t const regjext = regext[1];
    size_t const regkext = regext[2];
    
    size_t const srcioff = srcoff[0];
    size_t const srcjoff = srcoff[1];
    size_t const srckoff = srcoff[2];
    
    size_t const dstioff = dstoff[0];
    size_t const dstjoff = dstoff[1];
    size_t const dstkoff = dstoff[2];
    
    
    
    size_t const fi = srcioff % 2;
    size_t const fj = srcjoff % 2;
    size_t const fk = srckoff % 2;
    
    size_t const i0 = srcioff / 2;
    size_t const j0 = srcjoff / 2;
    size_t const k0 = srckoff / 2;
    
    //RT const one = 1;
    //
    //RT const f1 =    3*one/256;
    //RT const f2 = - 25*one/256;
    //RT const f3 =  150*one/256;
    //RT const f4 =  150*one/256;
    //RT const f5 = - 25*one/256;
    //RT const f6 =    3*one/256;

    double f1di[5], f1dj[5], f1dk[5];

    // Loop over fine region
    // Label scheme: l 8 fk fj fi
    
    size_t is, js, ks;
    size_t id, jd, kd;
    size_t i, j, k;
    
    // begin k loop
    k = 0;
    ks = k0;
    kd = dstkoff;
    if (fk == 0) goto l80;
    goto l81;
    
    // begin j loop
   l80:
    j = 0;
    js = j0;
    jd = dstjoff;
    if (fj == 0) goto l800;
    goto l801;
    
    // begin i loop
   l800:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8000;
    goto l8001;
    
    // kernel
   l8000:
    //dst[DSTIND3(id,jd,kd)] = src[SRCIND3(is,js,ks)];

    // Interpolate to (is-1/4, js-1/4, ks) 
    for (int jj=0; jj<5; jj++) {
        for (int ii=0; ii<5; ii++) {
            f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks)];
        }
        f1dj[jj] = weno_5pt(f1di,-0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dj,-0.25,0.0,1.0);

//printf("l8000: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);


    i = i+1;
    id = id+1;
    if (i < regiext) goto l8001;
    goto l900;
    
    // kernel
   l8001:
    //dst[DSTIND3(id,jd,kd)] =
    //  + f1 * src[SRCIND3(is-2,js,ks)]
    //  + f2 * src[SRCIND3(is-1,js,ks)]
    //  + f3 * src[SRCIND3(is  ,js,ks)]
    //  + f4 * src[SRCIND3(is+1,js,ks)]
    //  + f5 * src[SRCIND3(is+2,js,ks)]
    //  + f6 * src[SRCIND3(is+3,js,ks)];

    // Interpolate to (is+1/4,js-1/4,ks)
    for (int jj=0; jj<5; jj++) {
        for (int ii=0; ii<5; ii++) {
            f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks)];
        }
        f1dj[jj] = weno_5pt(f1di,0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dj,-0.25,0.0,1.0);

//printf("l8001: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8000;
    goto l900;
    
    // end i loop
   l900:
    j = j+1;
    jd = jd+1;
    if (j < regjext) goto l801;
    goto l90;
    
    // begin i loop
   l801:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8010;
    goto l8011;
    
    // kernel
   l8010:
    //dst[DSTIND3(id,jd,kd)] =
    //  + f1 * src[SRCIND3(is,js-2,ks)]
    //  + f2 * src[SRCIND3(is,js-1,ks)]
    //  + f3 * src[SRCIND3(is,js  ,ks)]
    //  + f4 * src[SRCIND3(is,js+1,ks)]
    //  + f5 * src[SRCIND3(is,js+2,ks)]
    //  + f6 * src[SRCIND3(is,js+3,ks)];

    // Interpolate to (is-1/4,js+1/4,ks)
    for (int jj=0; jj<5; jj++) {
        for (int ii=0; ii<5; ii++) {
            f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks)];
        }
        f1dj[jj] = weno_5pt(f1di,-0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dj,0.25,0.0,1.0);

//printf("l8010: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    if (i < regiext) goto l8011;
    goto l901;
    
    // kernel
   l8011:
    //dst[DSTIND3(id,jd,kd)] =
    //  + f1*f1 * src[SRCIND3(is-2,js-2,ks)]
    //  + f2*f1 * src[SRCIND3(is-1,js-2,ks)]
    //  + f3*f1 * src[SRCIND3(is  ,js-2,ks)]
    //  + f4*f1 * src[SRCIND3(is+1,js-2,ks)]
    //  + f5*f1 * src[SRCIND3(is+2,js-2,ks)]
    //  + f6*f1 * src[SRCIND3(is+3,js-2,ks)]
    //  + f1*f2 * src[SRCIND3(is-2,js-1,ks)]
    //  + f2*f2 * src[SRCIND3(is-1,js-1,ks)]
    //  + f3*f2 * src[SRCIND3(is  ,js-1,ks)]
    //  + f4*f2 * src[SRCIND3(is+1,js-1,ks)]
    //  + f5*f2 * src[SRCIND3(is+2,js-1,ks)]
    //  + f6*f2 * src[SRCIND3(is+3,js-1,ks)]
    //  + f1*f3 * src[SRCIND3(is-2,js  ,ks)]
    //  + f2*f3 * src[SRCIND3(is-1,js  ,ks)]
    //  + f3*f3 * src[SRCIND3(is  ,js  ,ks)]
    //  + f4*f3 * src[SRCIND3(is+1,js  ,ks)]
    //  + f5*f3 * src[SRCIND3(is+2,js  ,ks)]
    //  + f6*f3 * src[SRCIND3(is+3,js  ,ks)]
    //  + f1*f4 * src[SRCIND3(is-2,js+1,ks)]
    //  + f2*f4 * src[SRCIND3(is-1,js+1,ks)]
    //  + f3*f4 * src[SRCIND3(is  ,js+1,ks)]
    //  + f4*f4 * src[SRCIND3(is+1,js+1,ks)]
    //  + f5*f4 * src[SRCIND3(is+2,js+1,ks)]
    //  + f6*f4 * src[SRCIND3(is+3,js+1,ks)]
    //  + f1*f5 * src[SRCIND3(is-2,js+2,ks)]
    //  + f2*f5 * src[SRCIND3(is-1,js+2,ks)]
    //  + f3*f5 * src[SRCIND3(is  ,js+2,ks)]
    //  + f4*f5 * src[SRCIND3(is+1,js+2,ks)]
    //  + f5*f5 * src[SRCIND3(is+2,js+2,ks)]
    //  + f6*f5 * src[SRCIND3(is+3,js+2,ks)]
    //  + f1*f6 * src[SRCIND3(is-2,js+3,ks)]
    //  + f2*f6 * src[SRCIND3(is-1,js+3,ks)]
    //  + f3*f6 * src[SRCIND3(is  ,js+3,ks)]
    //  + f4*f6 * src[SRCIND3(is+1,js+3,ks)]
    //  + f5*f6 * src[SRCIND3(is+2,js+3,ks)]
    //  + f6*f6 * src[SRCIND3(is+3,js+3,ks)];

    // Interpolate to (is+1/4,js+1/4,ks)
    for (int jj=0; jj<5; jj++) {
        for (int ii=0; ii<5; ii++) {
            f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks)];
        }
        f1dj[jj] = weno_5pt(f1di,0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dj,0.25,0.0,1.0);

//printf("l8011: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8010;
    goto l901;
    
    // end i loop
   l901:
    j = j+1;
    jd = jd+1;
    js = js+1;
    if (j < regjext) goto l800;
    goto l90;
    
    // end j loop
   l90:
    k = k+1;
    kd = kd+1;
    if (k < regkext) goto l81;
    goto l9;
    
    // begin j loop
   l81:
    j = 0;
    js = j0;
    jd = dstjoff;
    if (fj == 0) goto l810;
    goto l811;
    
    // begin i loop
   l810:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8100;
    goto l8101;
    
    // kernel
   l8100:
    //dst[DSTIND3(id,jd,kd)] =
    //  + f1 * src[SRCIND3(is,js,ks-2)]
    //  + f2 * src[SRCIND3(is,js,ks-1)]
    //  + f3 * src[SRCIND3(is,js,ks  )]
    //  + f4 * src[SRCIND3(is,js,ks+1)]
    //  + f5 * src[SRCIND3(is,js,ks+2)]
    //  + f6 * src[SRCIND3(is,js,ks+3)];

    // Interpolate to (is-1/4,js-1/4,ks+1/2)
    for (int kk=0; kk<5; kk++) {
        for(int jj=0; jj<5; jj++) {
           for (int ii=0; ii<5; ii++) {
               f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks+kk-2)];
           }
           f1dj[jj] = weno_5pt(f1di,-0.25,0.0,1.0);
        }
        f1dk[kk] = weno_5pt(f1dj,-0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dk,0.5,0.0,1.0);

//printf("l8100: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    if (i < regiext) goto l8101;
    goto l910;
    
    // kernel
   l8101:
    //dst[DSTIND3(id,jd,kd)] =
    //  + f1*f1 * src[SRCIND3(is-2,js,ks-2)]
    //  + f2*f1 * src[SRCIND3(is-1,js,ks-2)]
    //  + f3*f1 * src[SRCIND3(is  ,js,ks-2)]
    //  + f4*f1 * src[SRCIND3(is+1,js,ks-2)]
    //  + f5*f1 * src[SRCIND3(is+2,js,ks-2)]
    //  + f6*f1 * src[SRCIND3(is+3,js,ks-2)]
    //  + f1*f2 * src[SRCIND3(is-2,js,ks-1)]
    //  + f2*f2 * src[SRCIND3(is-1,js,ks-1)]
    //  + f3*f2 * src[SRCIND3(is  ,js,ks-1)]
    //  + f4*f2 * src[SRCIND3(is+1,js,ks-1)]
    //  + f5*f2 * src[SRCIND3(is+2,js,ks-1)]
    //  + f6*f2 * src[SRCIND3(is+3,js,ks-1)]
    //  + f1*f3 * src[SRCIND3(is-2,js,ks  )]
    //  + f2*f3 * src[SRCIND3(is-1,js,ks  )]
    //  + f3*f3 * src[SRCIND3(is  ,js,ks  )]
    //  + f4*f3 * src[SRCIND3(is+1,js,ks  )]
    //  + f5*f3 * src[SRCIND3(is+2,js,ks  )]
    //  + f6*f3 * src[SRCIND3(is+3,js,ks  )]
    //  + f1*f4 * src[SRCIND3(is-2,js,ks+1)]
    //  + f2*f4 * src[SRCIND3(is-1,js,ks+1)]
    //  + f3*f4 * src[SRCIND3(is  ,js,ks+1)]
    //  + f4*f4 * src[SRCIND3(is+1,js,ks+1)]
    //  + f5*f4 * src[SRCIND3(is+2,js,ks+1)]
    //  + f6*f4 * src[SRCIND3(is+3,js,ks+1)]
    //  + f1*f5 * src[SRCIND3(is-2,js,ks+2)]
    //  + f2*f5 * src[SRCIND3(is-1,js,ks+2)]
    //  + f3*f5 * src[SRCIND3(is  ,js,ks+2)]
    //  + f4*f5 * src[SRCIND3(is+1,js,ks+2)]
    //  + f5*f5 * src[SRCIND3(is+2,js,ks+2)]
    //  + f6*f5 * src[SRCIND3(is+3,js,ks+2)]
    //  + f1*f6 * src[SRCIND3(is-2,js,ks+3)]
    //  + f2*f6 * src[SRCIND3(is-1,js,ks+3)]
    //  + f3*f6 * src[SRCIND3(is  ,js,ks+3)]
    //  + f4*f6 * src[SRCIND3(is+1,js,ks+3)]
    //  + f5*f6 * src[SRCIND3(is+2,js,ks+3)]
    //  + f6*f6 * src[SRCIND3(is+3,js,ks+3)];

    // Interpolate to (is+1/4,js-1/4,ks+1/2)
    for (int kk=0; kk<5; kk++) {
        for(int jj=0; jj<5; jj++) {
           for (int ii=0; ii<5; ii++) {
               f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks+kk-2)];
           }
           f1dj[jj] = weno_5pt(f1di,0.25,0.0,1.0);
        }
        f1dk[kk] = weno_5pt(f1dj,-0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dk,0.5,0.0,1.0);

//printf("l8101: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8100;
    goto l910;
    
    // end i loop
   l910:
    j = j+1;
    jd = jd+1;
    if (j < regjext) goto l811;
    goto l91;
    
    // begin i loop
   l811:
    i = 0;
    is = i0;
    id = dstioff;
    if (fi == 0) goto l8110;
    goto l8111;
    
    // kernel
   l8110:
    //dst[DSTIND3(id,jd,kd)] =
    //  + f1*f1 * src[SRCIND3(is,js-2,ks-2)]
    //  + f2*f1 * src[SRCIND3(is,js-1,ks-2)]
    //  + f3*f1 * src[SRCIND3(is,js  ,ks-2)]
    //  + f4*f1 * src[SRCIND3(is,js+1,ks-2)]
    //  + f5*f1 * src[SRCIND3(is,js+2,ks-2)]
    //  + f6*f1 * src[SRCIND3(is,js+3,ks-2)]
    //  + f1*f2 * src[SRCIND3(is,js-2,ks-1)]
    //  + f2*f2 * src[SRCIND3(is,js-1,ks-1)]
    //  + f3*f2 * src[SRCIND3(is,js  ,ks-1)]
    //  + f4*f2 * src[SRCIND3(is,js+1,ks-1)]
    //  + f5*f2 * src[SRCIND3(is,js+2,ks-1)]
    //  + f6*f2 * src[SRCIND3(is,js+3,ks-1)]
    //  + f1*f3 * src[SRCIND3(is,js-2,ks  )]
    //  + f2*f3 * src[SRCIND3(is,js-1,ks  )]
    //  + f3*f3 * src[SRCIND3(is,js  ,ks  )]
    //  + f4*f3 * src[SRCIND3(is,js+1,ks  )]
    //  + f5*f3 * src[SRCIND3(is,js+2,ks  )]
    //  + f6*f3 * src[SRCIND3(is,js+3,ks  )]
    //  + f1*f4 * src[SRCIND3(is,js-2,ks+1)]
    //  + f2*f4 * src[SRCIND3(is,js-1,ks+1)]
    //  + f3*f4 * src[SRCIND3(is,js  ,ks+1)]
    //  + f4*f4 * src[SRCIND3(is,js+1,ks+1)]
    //  + f5*f4 * src[SRCIND3(is,js+2,ks+1)]
    //  + f6*f4 * src[SRCIND3(is,js+3,ks+1)]
    //  + f1*f5 * src[SRCIND3(is,js-2,ks+2)]
    //  + f2*f5 * src[SRCIND3(is,js-1,ks+2)]
    //  + f3*f5 * src[SRCIND3(is,js  ,ks+2)]
    //  + f4*f5 * src[SRCIND3(is,js+1,ks+2)]
    //  + f5*f5 * src[SRCIND3(is,js+2,ks+2)]
    //  + f6*f5 * src[SRCIND3(is,js+3,ks+2)]
    //  + f1*f6 * src[SRCIND3(is,js-2,ks+3)]
    //  + f2*f6 * src[SRCIND3(is,js-1,ks+3)]
    //  + f3*f6 * src[SRCIND3(is,js  ,ks+3)]
    //  + f4*f6 * src[SRCIND3(is,js+1,ks+3)]
    //  + f5*f6 * src[SRCIND3(is,js+2,ks+3)]
    //  + f6*f6 * src[SRCIND3(is,js+3,ks+3)];

    // Interpolate to (is-1/4,js+1/4,ks+1/2)
    for (int kk=0; kk<5; kk++) {
        for(int jj=0; jj<5; jj++) {
           for (int ii=0; ii<5; ii++) {
               f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks+kk-2)];
           }
           f1dj[jj] = weno_5pt(f1di,-0.25,0.0,1.0);
        }
        f1dk[kk] = weno_5pt(f1dj,0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dk,0.5,0.0,1.0);

//printf("l8110: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    if (i < regiext) goto l8111;
    goto l911;
    
    // kernel
   l8111:
    //{
    //  T const res1 =
    //    + f1*f1*f1 * src[SRCIND3(is-2,js-2,ks-2)]
    //    + f2*f1*f1 * src[SRCIND3(is-1,js-2,ks-2)]
    //    + f3*f1*f1 * src[SRCIND3(is  ,js-2,ks-2)]
    //    + f4*f1*f1 * src[SRCIND3(is+1,js-2,ks-2)]
    //    + f5*f1*f1 * src[SRCIND3(is+2,js-2,ks-2)]
    //    + f6*f1*f1 * src[SRCIND3(is+3,js-2,ks-2)]
    //    + f1*f2*f1 * src[SRCIND3(is-2,js-1,ks-2)]
    //    + f2*f2*f1 * src[SRCIND3(is-1,js-1,ks-2)]
    //    + f3*f2*f1 * src[SRCIND3(is  ,js-1,ks-2)]
    //    + f4*f2*f1 * src[SRCIND3(is+1,js-1,ks-2)]
    //    + f5*f2*f1 * src[SRCIND3(is+2,js-1,ks-2)]
    //    + f6*f2*f1 * src[SRCIND3(is+3,js-1,ks-2)]
    //    + f1*f3*f1 * src[SRCIND3(is-2,js  ,ks-2)]
    //    + f2*f3*f1 * src[SRCIND3(is-1,js  ,ks-2)]
    //    + f3*f3*f1 * src[SRCIND3(is  ,js  ,ks-2)]
    //    + f4*f3*f1 * src[SRCIND3(is+1,js  ,ks-2)]
    //    + f5*f3*f1 * src[SRCIND3(is+2,js  ,ks-2)]
    //    + f6*f3*f1 * src[SRCIND3(is+3,js  ,ks-2)]
    //    + f1*f4*f1 * src[SRCIND3(is-2,js+1,ks-2)]
    //    + f2*f4*f1 * src[SRCIND3(is-1,js+1,ks-2)]
    //    + f3*f4*f1 * src[SRCIND3(is  ,js+1,ks-2)]
    //    + f4*f4*f1 * src[SRCIND3(is+1,js+1,ks-2)]
    //    + f5*f4*f1 * src[SRCIND3(is+2,js+1,ks-2)]
    //    + f6*f4*f1 * src[SRCIND3(is+3,js+1,ks-2)]
    //    + f1*f5*f1 * src[SRCIND3(is-2,js+2,ks-2)]
    //    + f2*f5*f1 * src[SRCIND3(is-1,js+2,ks-2)]
    //    + f3*f5*f1 * src[SRCIND3(is  ,js+2,ks-2)]
    //    + f4*f5*f1 * src[SRCIND3(is+1,js+2,ks-2)]
    //    + f5*f5*f1 * src[SRCIND3(is+2,js+2,ks-2)]
    //    + f6*f5*f1 * src[SRCIND3(is+3,js+2,ks-2)]
    //    + f1*f6*f1 * src[SRCIND3(is-2,js+3,ks-2)]
    //    + f2*f6*f1 * src[SRCIND3(is-1,js+3,ks-2)]
    //    + f3*f6*f1 * src[SRCIND3(is  ,js+3,ks-2)]
    //    + f4*f6*f1 * src[SRCIND3(is+1,js+3,ks-2)]
    //    + f5*f6*f1 * src[SRCIND3(is+2,js+3,ks-2)]
    //    + f6*f6*f1 * src[SRCIND3(is+3,js+3,ks-2)];
    //  T const res2 =
    //    + f1*f1*f2 * src[SRCIND3(is-2,js-2,ks-1)]
    //    + f2*f1*f2 * src[SRCIND3(is-1,js-2,ks-1)]
    //    + f3*f1*f2 * src[SRCIND3(is  ,js-2,ks-1)]
    //    + f4*f1*f2 * src[SRCIND3(is+1,js-2,ks-1)]
    //    + f5*f1*f2 * src[SRCIND3(is+2,js-2,ks-1)]
    //    + f6*f1*f2 * src[SRCIND3(is+3,js-2,ks-1)]
    //    + f1*f2*f2 * src[SRCIND3(is-2,js-1,ks-1)]
    //    + f2*f2*f2 * src[SRCIND3(is-1,js-1,ks-1)]
    //    + f3*f2*f2 * src[SRCIND3(is  ,js-1,ks-1)]
    //    + f4*f2*f2 * src[SRCIND3(is+1,js-1,ks-1)]
    //    + f5*f2*f2 * src[SRCIND3(is+2,js-1,ks-1)]
    //    + f6*f2*f2 * src[SRCIND3(is+3,js-1,ks-1)]
    //    + f1*f3*f2 * src[SRCIND3(is-2,js  ,ks-1)]
    //    + f2*f3*f2 * src[SRCIND3(is-1,js  ,ks-1)]
    //    + f3*f3*f2 * src[SRCIND3(is  ,js  ,ks-1)]
    //    + f4*f3*f2 * src[SRCIND3(is+1,js  ,ks-1)]
    //    + f5*f3*f2 * src[SRCIND3(is+2,js  ,ks-1)]
    //    + f6*f3*f2 * src[SRCIND3(is+3,js  ,ks-1)]
    //    + f1*f4*f2 * src[SRCIND3(is-2,js+1,ks-1)]
    //    + f2*f4*f2 * src[SRCIND3(is-1,js+1,ks-1)]
    //    + f3*f4*f2 * src[SRCIND3(is  ,js+1,ks-1)]
    //    + f4*f4*f2 * src[SRCIND3(is+1,js+1,ks-1)]
    //    + f5*f4*f2 * src[SRCIND3(is+2,js+1,ks-1)]
    //    + f6*f4*f2 * src[SRCIND3(is+3,js+1,ks-1)]
    //    + f1*f5*f2 * src[SRCIND3(is-2,js+2,ks-1)]
    //    + f2*f5*f2 * src[SRCIND3(is-1,js+2,ks-1)]
    //    + f3*f5*f2 * src[SRCIND3(is  ,js+2,ks-1)]
    //    + f4*f5*f2 * src[SRCIND3(is+1,js+2,ks-1)]
    //    + f5*f5*f2 * src[SRCIND3(is+2,js+2,ks-1)]
    //    + f6*f5*f2 * src[SRCIND3(is+3,js+2,ks-1)]
    //    + f1*f6*f2 * src[SRCIND3(is-2,js+3,ks-1)]
    //    + f2*f6*f2 * src[SRCIND3(is-1,js+3,ks-1)]
    //    + f3*f6*f2 * src[SRCIND3(is  ,js+3,ks-1)]
    //    + f4*f6*f2 * src[SRCIND3(is+1,js+3,ks-1)]
    //    + f5*f6*f2 * src[SRCIND3(is+2,js+3,ks-1)]
    //    + f6*f6*f2 * src[SRCIND3(is+3,js+3,ks-1)];
    //  T const res3 =
    //    + f1*f1*f3 * src[SRCIND3(is-2,js-2,ks  )]
    //    + f2*f1*f3 * src[SRCIND3(is-1,js-2,ks  )]
    //    + f3*f1*f3 * src[SRCIND3(is  ,js-2,ks  )]
    //    + f4*f1*f3 * src[SRCIND3(is+1,js-2,ks  )]
    //    + f5*f1*f3 * src[SRCIND3(is+2,js-2,ks  )]
    //    + f6*f1*f3 * src[SRCIND3(is+3,js-2,ks  )]
    //    + f1*f2*f3 * src[SRCIND3(is-2,js-1,ks  )]
    //    + f2*f2*f3 * src[SRCIND3(is-1,js-1,ks  )]
    //    + f3*f2*f3 * src[SRCIND3(is  ,js-1,ks  )]
    //    + f4*f2*f3 * src[SRCIND3(is+1,js-1,ks  )]
    //    + f5*f2*f3 * src[SRCIND3(is+2,js-1,ks  )]
    //    + f6*f2*f3 * src[SRCIND3(is+3,js-1,ks  )]
    //    + f1*f3*f3 * src[SRCIND3(is-2,js  ,ks  )]
    //    + f2*f3*f3 * src[SRCIND3(is-1,js  ,ks  )]
    //    + f3*f3*f3 * src[SRCIND3(is  ,js  ,ks  )]
    //    + f4*f3*f3 * src[SRCIND3(is+1,js  ,ks  )]
    //    + f5*f3*f3 * src[SRCIND3(is+2,js  ,ks  )]
    //    + f6*f3*f3 * src[SRCIND3(is+3,js  ,ks  )]
    //    + f1*f4*f3 * src[SRCIND3(is-2,js+1,ks  )]
    //    + f2*f4*f3 * src[SRCIND3(is-1,js+1,ks  )]
    //    + f3*f4*f3 * src[SRCIND3(is  ,js+1,ks  )]
    //    + f4*f4*f3 * src[SRCIND3(is+1,js+1,ks  )]
    //    + f5*f4*f3 * src[SRCIND3(is+2,js+1,ks  )]
    //    + f6*f4*f3 * src[SRCIND3(is+3,js+1,ks  )]
    //    + f1*f5*f3 * src[SRCIND3(is-2,js+2,ks  )]
    //    + f2*f5*f3 * src[SRCIND3(is-1,js+2,ks  )]
    //    + f3*f5*f3 * src[SRCIND3(is  ,js+2,ks  )]
    //    + f4*f5*f3 * src[SRCIND3(is+1,js+2,ks  )]
    //    + f5*f5*f3 * src[SRCIND3(is+2,js+2,ks  )]
    //    + f6*f5*f3 * src[SRCIND3(is+3,js+2,ks  )]
    //    + f1*f6*f3 * src[SRCIND3(is-2,js+3,ks  )]
    //    + f2*f6*f3 * src[SRCIND3(is-1,js+3,ks  )]
    //    + f3*f6*f3 * src[SRCIND3(is  ,js+3,ks  )]
    //    + f4*f6*f3 * src[SRCIND3(is+1,js+3,ks  )]
    //    + f5*f6*f3 * src[SRCIND3(is+2,js+3,ks  )]
    //    + f6*f6*f3 * src[SRCIND3(is+3,js+3,ks  )];
    //  T const res4 =
    //    + f1*f1*f4 * src[SRCIND3(is-2,js-2,ks+1)]
    //    + f2*f1*f4 * src[SRCIND3(is-1,js-2,ks+1)]
    //    + f3*f1*f4 * src[SRCIND3(is  ,js-2,ks+1)]
    //    + f4*f1*f4 * src[SRCIND3(is+1,js-2,ks+1)]
    //    + f5*f1*f4 * src[SRCIND3(is+2,js-2,ks+1)]
    //    + f6*f1*f4 * src[SRCIND3(is+3,js-2,ks+1)]
    //    + f1*f2*f4 * src[SRCIND3(is-2,js-1,ks+1)]
    //    + f2*f2*f4 * src[SRCIND3(is-1,js-1,ks+1)]
    //    + f3*f2*f4 * src[SRCIND3(is  ,js-1,ks+1)]
    //    + f4*f2*f4 * src[SRCIND3(is+1,js-1,ks+1)]
    //    + f5*f2*f4 * src[SRCIND3(is+2,js-1,ks+1)]
    //    + f6*f2*f4 * src[SRCIND3(is+3,js-1,ks+1)]
    //    + f1*f3*f4 * src[SRCIND3(is-2,js  ,ks+1)]
    //    + f2*f3*f4 * src[SRCIND3(is-1,js  ,ks+1)]
    //    + f3*f3*f4 * src[SRCIND3(is  ,js  ,ks+1)]
    //    + f4*f3*f4 * src[SRCIND3(is+1,js  ,ks+1)]
    //    + f5*f3*f4 * src[SRCIND3(is+2,js  ,ks+1)]
    //    + f6*f3*f4 * src[SRCIND3(is+3,js  ,ks+1)]
    //    + f1*f4*f4 * src[SRCIND3(is-2,js+1,ks+1)]
    //    + f2*f4*f4 * src[SRCIND3(is-1,js+1,ks+1)]
    //    + f3*f4*f4 * src[SRCIND3(is  ,js+1,ks+1)]
    //    + f4*f4*f4 * src[SRCIND3(is+1,js+1,ks+1)]
    //    + f5*f4*f4 * src[SRCIND3(is+2,js+1,ks+1)]
    //    + f6*f4*f4 * src[SRCIND3(is+3,js+1,ks+1)]
    //    + f1*f5*f4 * src[SRCIND3(is-2,js+2,ks+1)]
    //    + f2*f5*f4 * src[SRCIND3(is-1,js+2,ks+1)]
    //    + f3*f5*f4 * src[SRCIND3(is  ,js+2,ks+1)]
    //    + f4*f5*f4 * src[SRCIND3(is+1,js+2,ks+1)]
    //    + f5*f5*f4 * src[SRCIND3(is+2,js+2,ks+1)]
    //    + f6*f5*f4 * src[SRCIND3(is+3,js+2,ks+1)]
    //    + f1*f6*f4 * src[SRCIND3(is-2,js+3,ks+1)]
    //    + f2*f6*f4 * src[SRCIND3(is-1,js+3,ks+1)]
    //    + f3*f6*f4 * src[SRCIND3(is  ,js+3,ks+1)]
    //    + f4*f6*f4 * src[SRCIND3(is+1,js+3,ks+1)]
    //    + f5*f6*f4 * src[SRCIND3(is+2,js+3,ks+1)]
    //    + f6*f6*f4 * src[SRCIND3(is+3,js+3,ks+1)];
    //  T const res5 =
    //    + f1*f1*f5 * src[SRCIND3(is-2,js-2,ks+2)]
    //    + f2*f1*f5 * src[SRCIND3(is-1,js-2,ks+2)]
    //    + f3*f1*f5 * src[SRCIND3(is  ,js-2,ks+2)]
    //    + f4*f1*f5 * src[SRCIND3(is+1,js-2,ks+2)]
    //    + f5*f1*f5 * src[SRCIND3(is+2,js-2,ks+2)]
    //    + f6*f1*f5 * src[SRCIND3(is+3,js-2,ks+2)]
    //    + f1*f2*f5 * src[SRCIND3(is-2,js-1,ks+2)]
    //    + f2*f2*f5 * src[SRCIND3(is-1,js-1,ks+2)]
    //    + f3*f2*f5 * src[SRCIND3(is  ,js-1,ks+2)]
    //    + f4*f2*f5 * src[SRCIND3(is+1,js-1,ks+2)]
    //    + f5*f2*f5 * src[SRCIND3(is+2,js-1,ks+2)]
    //    + f6*f2*f5 * src[SRCIND3(is+3,js-1,ks+2)]
    //    + f1*f3*f5 * src[SRCIND3(is-2,js  ,ks+2)]
    //    + f2*f3*f5 * src[SRCIND3(is-1,js  ,ks+2)]
    //    + f3*f3*f5 * src[SRCIND3(is  ,js  ,ks+2)]
    //    + f4*f3*f5 * src[SRCIND3(is+1,js  ,ks+2)]
    //    + f5*f3*f5 * src[SRCIND3(is+2,js  ,ks+2)]
    //    + f6*f3*f5 * src[SRCIND3(is+3,js  ,ks+2)]
    //    + f1*f4*f5 * src[SRCIND3(is-2,js+1,ks+2)]
    //    + f2*f4*f5 * src[SRCIND3(is-1,js+1,ks+2)]
    //    + f3*f4*f5 * src[SRCIND3(is  ,js+1,ks+2)]
    //    + f4*f4*f5 * src[SRCIND3(is+1,js+1,ks+2)]
    //    + f5*f4*f5 * src[SRCIND3(is+2,js+1,ks+2)]
    //    + f6*f4*f5 * src[SRCIND3(is+3,js+1,ks+2)]
    //    + f1*f5*f5 * src[SRCIND3(is-2,js+2,ks+2)]
    //    + f2*f5*f5 * src[SRCIND3(is-1,js+2,ks+2)]
    //    + f3*f5*f5 * src[SRCIND3(is  ,js+2,ks+2)]
    //    + f4*f5*f5 * src[SRCIND3(is+1,js+2,ks+2)]
    //    + f5*f5*f5 * src[SRCIND3(is+2,js+2,ks+2)]
    //    + f6*f5*f5 * src[SRCIND3(is+3,js+2,ks+2)]
    //    + f1*f6*f5 * src[SRCIND3(is-2,js+3,ks+2)]
    //    + f2*f6*f5 * src[SRCIND3(is-1,js+3,ks+2)]
    //    + f3*f6*f5 * src[SRCIND3(is  ,js+3,ks+2)]
    //    + f4*f6*f5 * src[SRCIND3(is+1,js+3,ks+2)]
    //    + f5*f6*f5 * src[SRCIND3(is+2,js+3,ks+2)]
    //    + f6*f6*f5 * src[SRCIND3(is+3,js+3,ks+2)];
    //  T const res6 =
    //    + f1*f1*f6 * src[SRCIND3(is-2,js-2,ks+3)]
    //    + f2*f1*f6 * src[SRCIND3(is-1,js-2,ks+3)]
    //    + f3*f1*f6 * src[SRCIND3(is  ,js-2,ks+3)]
    //    + f4*f1*f6 * src[SRCIND3(is+1,js-2,ks+3)]
    //    + f5*f1*f6 * src[SRCIND3(is+2,js-2,ks+3)]
    //    + f6*f1*f6 * src[SRCIND3(is+3,js-2,ks+3)]
    //    + f1*f2*f6 * src[SRCIND3(is-2,js-1,ks+3)]
    //    + f2*f2*f6 * src[SRCIND3(is-1,js-1,ks+3)]
    //    + f3*f2*f6 * src[SRCIND3(is  ,js-1,ks+3)]
    //    + f4*f2*f6 * src[SRCIND3(is+1,js-1,ks+3)]
    //    + f5*f2*f6 * src[SRCIND3(is+2,js-1,ks+3)]
    //    + f6*f2*f6 * src[SRCIND3(is+3,js-1,ks+3)]
    //    + f1*f3*f6 * src[SRCIND3(is-2,js  ,ks+3)]
    //    + f2*f3*f6 * src[SRCIND3(is-1,js  ,ks+3)]
    //    + f3*f3*f6 * src[SRCIND3(is  ,js  ,ks+3)]
    //    + f4*f3*f6 * src[SRCIND3(is+1,js  ,ks+3)]
    //    + f5*f3*f6 * src[SRCIND3(is+2,js  ,ks+3)]
    //    + f6*f3*f6 * src[SRCIND3(is+3,js  ,ks+3)]
    //    + f1*f4*f6 * src[SRCIND3(is-2,js+1,ks+3)]
    //    + f2*f4*f6 * src[SRCIND3(is-1,js+1,ks+3)]
    //    + f3*f4*f6 * src[SRCIND3(is  ,js+1,ks+3)]
    //    + f4*f4*f6 * src[SRCIND3(is+1,js+1,ks+3)]
    //    + f5*f4*f6 * src[SRCIND3(is+2,js+1,ks+3)]
    //    + f6*f4*f6 * src[SRCIND3(is+3,js+1,ks+3)]
    //    + f1*f5*f6 * src[SRCIND3(is-2,js+2,ks+3)]
    //    + f2*f5*f6 * src[SRCIND3(is-1,js+2,ks+3)]
    //    + f3*f5*f6 * src[SRCIND3(is  ,js+2,ks+3)]
    //    + f4*f5*f6 * src[SRCIND3(is+1,js+2,ks+3)]
    //    + f5*f5*f6 * src[SRCIND3(is+2,js+2,ks+3)]
    //    + f6*f5*f6 * src[SRCIND3(is+3,js+2,ks+3)]
    //    + f1*f6*f6 * src[SRCIND3(is-2,js+3,ks+3)]
    //    + f2*f6*f6 * src[SRCIND3(is-1,js+3,ks+3)]
    //    + f3*f6*f6 * src[SRCIND3(is  ,js+3,ks+3)]
    //    + f4*f6*f6 * src[SRCIND3(is+1,js+3,ks+3)]
    //    + f5*f6*f6 * src[SRCIND3(is+2,js+3,ks+3)]
    //    + f6*f6*f6 * src[SRCIND3(is+3,js+3,ks+3)];
    //  dst[DSTIND3(id,jd,kd)] = res1 + res2 + res3 + res4 + res5 + res6;
    //}

    // Interpolate to (is+1/4,js+1/4,ks+1/2)
    for (int kk=0; kk<5; kk++) {
        for(int jj=0; jj<5; jj++) {
           for (int ii=0; ii<5; ii++) {
               f1di[ii] = src[SRCIND3(is+ii-2,js+jj-2,ks+kk-2)];
           }
           f1dj[jj] = weno_5pt(f1di,0.25,0.0,1.0);
        }
        f1dk[kk] = weno_5pt(f1dj,0.25,0.0,1.0);
    }
    dst[DSTIND3(id,jd,kd)] = weno_5pt(f1dk,0.5,0.0,1.0);

//printf("l8111: %d %d %d %d %d %d %e %e\n",id,jd,kd,is,js,ks,dst[DSTIND3(id,jd,kd)],src[SRCIND3(is,js,ks)]);

    i = i+1;
    id = id+1;
    is = is+1;
    if (i < regiext) goto l8110;
    goto l911;
    
    // end i loop
   l911:
    j = j+1;
    jd = jd+1;
    js = js+1;
    if (j < regjext) goto l810;
    goto l91;
    
    // end j loop
   l91:
    k = k+1;
    kd = kd+1;
    ks = ks+1;
    if (k < regkext) goto l80;
    goto l9;
    
    // end k loop
   l9:;
    
  }

#ifdef HAVE_CCTK_COMPLEX8
  template <>
  void
  prolongate_3d_stagger110 (CCTK_COMPLEX8 const * restrict const src,
                        ivect3 const & restrict srcext,
                        CCTK_COMPLEX8 * restrict const dst,
                        ivect3 const & restrict dstext,
                        ibbox3 const & restrict srcbbox,
                        ibbox3 const & restrict dstbbox,
                        ibbox3 const & restrict regbbox)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "STAGGER110 for complex numbers is not supported");
  }
#endif

#ifdef HAVE_CCTK_COMPLEX16
  template <>
  void
  prolongate_3d_stagger110 (CCTK_COMPLEX16 const * restrict const src,
                        ivect3 const & restrict srcext,
                        CCTK_COMPLEX16 * restrict const dst,
                        ivect3 const & restrict dstext,
                        ibbox3 const & restrict srcbbox,
                        ibbox3 const & restrict dstbbox,
                        ibbox3 const & restrict regbbox)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "STAGGER110 for complex numbers is not supported");
  }
#endif

#ifdef HAVE_CCTK_COMPLEX32
  template <>
  void
  prolongate_3d_stagger110 (CCTK_COMPLEX32 const * restrict const src,
                        ivect3 const & restrict srcext,
                        CCTK_COMPLEX32 * restrict const dst,
                        ivect3 const & restrict dstext,
                        ibbox3 const & restrict srcbbox,
                        ibbox3 const & restrict dstbbox,
                        ibbox3 const & restrict regbbox)
  {
    CCTK_WARN (CCTK_WARN_ABORT, "STAGGER110 for complex numbers is not supported");
  }
#endif

  
  
  
#define INSTANTIATE(T)                                          \
  template                                                      \
  void                                                          \
  prolongate_3d_stagger110 (T const * restrict const src,           \
                        ivect3 const & restrict srcext,         \
                        T * restrict const dst,                 \
                        ivect3 const & restrict dstext,         \
                        ibbox3 const & restrict srcbbox,        \
                        ibbox3 const & restrict dstbbox,        \
                        ibbox3 const & restrict regbbox);
#include "instantiate"
#undef INSTANTIATE
  
  
  
} // namespace CarpetLib

/*
  Restriction (fine -> coarse) for Ax (x-comp of vec. potential).

  ========
  Grid key:
  ========
  Interpolating fine -> coarse here:
  coarse(1) = fine(2,1,0)
  coarse(1.5) = fine(2,1,0)
  Staggered direction (y and z):
     * | *   * | *   * | *
interpolate to 2k+1/2, use points 2k-1 2k 2k+1

jf   0   1   2   3   4   5
jc     0       1       2       3       4

yc 0   .5     1.5

  | = coarse gridpoints
  * = fine gridpoints

  Unstaggered direction (x):
   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |
j  0   1   2   3   4   5   6
i  0       1       2       3
  | = coarse gridpoints
  | and * = fine gridpoints

  See below for further technical details.
*/

#include <algorithm>
#include <cassert>
#include <cmath>

#include <cctk.h>
#include <cctk_Parameters.h>

#include "operator_prototypes.hh"
#include "typeprops.hh"

#include "order_restrict.h"

using namespace std;

namespace CarpetLib {


  
#define SRCIND3(i,j,k)                                  \
  index3 (srcioff + (i), srcjoff + (j), srckoff + (k),  \
          srciext, srcjext, srckext)
#define DSTIND3(i,j,k)                                  \
  index3 (dstioff + (i), dstjoff + (j), dstkoff + (k),  \
          dstiext, dstjext, dstkext)
  
  
  
  template <typename T>
  void
  restrict_3d_stagger111 (T const * restrict const src,
			  ivect3 const & restrict srcext,
			  T * restrict const dst,
			  ivect3 const & restrict dstext,
			  ibbox3 const & restrict srcbbox,
			  ibbox3 const & restrict dstbbox,
			  ibbox3 const & restrict regbbox)
  {
    DECLARE_CCTK_PARAMETERS;
    
    typedef typename typeprops<T>::real RT;
    
    
    
    if (any (srcbbox.stride() >= regbbox.stride() or
             dstbbox.stride() != regbbox.stride()))
      {
	CCTK_WARN (0, "Internal error: strides disagree");
      }
    
    if (any (reffact2 * srcbbox.stride() != dstbbox.stride())) {
      CCTK_WARN (0, "Internal error: destination strides are not twice the source strides");
    }
    
    // This could be handled, but is likely to point to an error
    // elsewhere
    if (regbbox.empty()) {
      CCTK_WARN (0, "Internal error: region extent is empty");
    }
    
    if (not regbbox.expanded_for(srcbbox).is_contained_in(srcbbox) or
        not regbbox.is_contained_in(dstbbox))
      {
	CCTK_WARN (0, "Internal error: region extent is not contained in array extent");
      }
    
    if (any (srcext != srcbbox.shape() / srcbbox.stride() or
             dstext != dstbbox.shape() / dstbbox.stride()))
      {
	CCTK_WARN (0, "Internal error: array sizes don't agree with bounding boxes");
      }
    
    
    
    //ivect3 const regext = regbbox.shape() / regbbox.stride();
    //assert (all (srcbbox.stride() % 2 == 0));
    //assert (all ((regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) %
    //             srcbbox.stride() == 0));
    //ivect3 const srcoff =
    //  (regbbox.lower() - srcbbox.lower() - srcbbox.stride() / 2) /
    //  srcbbox.stride();
    //assert (all ((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
    //ivect3 const dstoff =
    //  (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();

    ivect3 const regext = regbbox.shape() / regbbox.stride();
    assert (all ((regbbox.lower() - srcbbox.lower()) % srcbbox.stride() == 0));
    ivect3 const srcoff = (regbbox.lower() - srcbbox.lower()) / srcbbox.stride();
    assert (all ((regbbox.lower() - dstbbox.lower()) % dstbbox.stride() == 0));
    ivect3 const dstoff = (regbbox.lower() - dstbbox.lower()) / dstbbox.stride();
    
    
    
    int const srciext = srcext[0];
    int const srcjext = srcext[1];
    int const srckext = srcext[2];
    
    int const dstiext = dstext[0];
    int const dstjext = dstext[1];
    int const dstkext = dstext[2];
    
    int const regiext = regext[0];
    int const regjext = regext[1];
    int const regkext = regext[2];
    
    int const srcioff = srcoff[0];
    int const srcjoff = srcoff[1];
    int const srckoff = srcoff[2];
    
    int const dstioff = dstoff[0];
    int const dstjoff = dstoff[1];
    int const dstkoff = dstoff[2];
    
#include "coeffs_restrict.h"

    // Loop over coarse region
#pragma omp parallel for
    for (int k=0; k<regkext; ++k) {
      for (int j=0; j<regjext; ++j) {
        for (int i=0; i<regiext; ++i) {

	  dst[DSTIND3(i, j, k)] = typeprops<T>::fromreal (0);

	  if(ORDER==2) {
	    for(int ii=1;ii<=ORDER+1;ii++) for(int jj=1;jj<=ORDER+1;jj++) for(int kk=1;kk<=ORDER+1;kk++) {
		  dst[DSTIND3(i, j, k)] += coeff[ii]*coeff[jj]*coeff[kk] * src [SRCIND3(2*i-1+(ii-1), 2*j-1+(jj-1), 2*k-1+(kk-1))];
		}
	  }

	  if(ORDER==3 || ORDER==4) {
	    for(int ii=1;ii<=ORDER+1;ii++) for(int jj=1;jj<=ORDER+1;jj++) for(int kk=1;kk<=ORDER+1;kk++) {
		  dst[DSTIND3(i, j, k)] += coeff[ii]*coeff[jj]*coeff[kk]*src [SRCIND3(2*i-2+(ii-1)+(4-ORDER), 2*j-2+(jj-1)+(4-ORDER), 2*k-2+(kk-1)+(4-ORDER))];
		}
	  }

        }
      }
    }
    
  }
  
#define INSTANTIATE(T)						\
  template							\
  void								\
  restrict_3d_stagger111 (T const * restrict const src,		\
			  ivect3 const & restrict srcext,	\
			  T * restrict const dst,		\
			  ivect3 const & restrict dstext,	\
			  ibbox3 const & restrict srcbbox,	\
			  ibbox3 const & restrict dstbbox,	\
			  ibbox3 const & restrict regbbox);

#include "instantiate"
#undef INSTANTIATE

} // namespace CarpetLib

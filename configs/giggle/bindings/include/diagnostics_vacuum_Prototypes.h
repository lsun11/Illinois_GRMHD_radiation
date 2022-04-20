/*@@
   @header  diagnostics_vacuum_Prototypes.h
   @author  Automatically generated by CreateFunctionBindings.pl
   @desc
            Prototypes for overloaded functions used by this thorn
   @enddesc
  @@*/


#ifndef _DIAGNOSTICS_VACUUM_PROTOTYPES_H_
#define _DIAGNOSTICS_VACUUM_PROTOTYPES_H_  1

#ifdef CCODE
#ifdef __cplusplus
extern "C" {
#endif

CCTK_INT HorizonCentroid(const CCTK_INT horizon_number ,
 CCTK_REAL* centroid_x ,
 CCTK_REAL* centroid_y ,
 CCTK_REAL* centroid_z);

CCTK_INT HorizonLocalCoordinateOrigin(const CCTK_INT horizon_number ,
 CCTK_REAL* origin_x ,
 CCTK_REAL* origin_y ,
 CCTK_REAL* origin_z);

CCTK_INT HorizonRadiusInDirection(const CCTK_INT horizon_number ,
 const CCTK_INT N_points ,
 const CCTK_REAL* x ,
 const CCTK_REAL* y ,
 const CCTK_REAL* z ,
 CCTK_REAL* radius);

CCTK_INT HorizonWasFound(const CCTK_INT horizon_number);

#ifdef __cplusplus
}
#endif
#endif /* CCODE */

#ifdef FCODE
#ifdef F90CODE
#define DECLARE_DIAGNOSTICS_VACUUM_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  interface &&\
     CCTK_INT function HorizonCentroid (horizon_number, centroid_x, centroid_y, centroid_z) &&\
       implicit none &&\
       CCTK_INT horizon_number &&\
       CCTK_REAL centroid_x &&\
       CCTK_REAL centroid_y &&\
       CCTK_REAL centroid_z &&\
     end function HorizonCentroid &&\
  end interface &&\
  interface &&\
     CCTK_INT function HorizonLocalCoordinateOrigin (horizon_number, origin_x, origin_y, origin_z) &&\
       implicit none &&\
       CCTK_INT horizon_number &&\
       CCTK_REAL origin_x &&\
       CCTK_REAL origin_y &&\
       CCTK_REAL origin_z &&\
     end function HorizonLocalCoordinateOrigin &&\
  end interface &&\
  interface &&\
     CCTK_INT function HorizonRadiusInDirection (horizon_number, N_points, x, y, z, radius) &&\
       implicit none &&\
       CCTK_INT horizon_number &&\
       CCTK_INT N_points &&\
       CCTK_REAL x(*) &&\
       CCTK_REAL y(*) &&\
       CCTK_REAL z(*) &&\
       CCTK_REAL radius(*) &&\
     end function HorizonRadiusInDirection &&\
  end interface &&\
  interface &&\
     CCTK_INT function HorizonWasFound (horizon_number) &&\
       implicit none &&\
       CCTK_INT horizon_number &&\
     end function HorizonWasFound &&\
  end interface &&\

#else /* ! F90CODE */

#define DECLARE_DIAGNOSTICS_VACUUM_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  external HorizonCentroid &&\
  CCTK_INT HorizonCentroid &&\
  external HorizonLocalCoordinateOrigin &&\
  CCTK_INT HorizonLocalCoordinateOrigin &&\
  external HorizonRadiusInDirection &&\
  CCTK_INT HorizonRadiusInDirection &&\
  external HorizonWasFound &&\
  CCTK_INT HorizonWasFound &&\

#endif /* ! F90CODE */
#endif /* FCODE */

#endif


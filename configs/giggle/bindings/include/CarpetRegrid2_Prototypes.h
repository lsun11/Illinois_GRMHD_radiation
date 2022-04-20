/*@@
   @header  CarpetRegrid2_Prototypes.h
   @author  Automatically generated by CreateFunctionBindings.pl
   @desc
            Prototypes for overloaded functions used by this thorn
   @enddesc
  @@*/


#ifndef _CARPETREGRID2_PROTOTYPES_H_
#define _CARPETREGRID2_PROTOTYPES_H_  1

#ifdef CCODE
#ifdef __cplusplus
extern "C" {
#endif

CCTK_INT Carpet_Regrid(const CCTK_POINTER_TO_CONST cctkGH ,
 const CCTK_POINTER superregss ,
 const CCTK_POINTER regsss ,
 const CCTK_INT force);
CCTK_INT CarpetRegrid2_Regrid(const CCTK_POINTER_TO_CONST cctkGH ,
 const CCTK_POINTER superregss ,
 const CCTK_POINTER regsss ,
 const CCTK_INT force);

CCTK_INT Carpet_RegridMaps(const CCTK_POINTER_TO_CONST cctkGH ,
 const CCTK_POINTER superregsss ,
 const CCTK_POINTER regssss ,
 const CCTK_INT force);
CCTK_INT CarpetRegrid2_RegridMaps(const CCTK_POINTER_TO_CONST cctkGH ,
 const CCTK_POINTER superregsss ,
 const CCTK_POINTER regssss ,
 const CCTK_INT force);

CCTK_INT ConvertFromPhysicalBoundary(const CCTK_INT size ,
 const CCTK_REAL* physical_min ,
 const CCTK_REAL* physical_max ,
 CCTK_REAL* interior_min ,
 CCTK_REAL* interior_max ,
 CCTK_REAL* exterior_min ,
 CCTK_REAL* exterior_max ,
 const CCTK_REAL* spacing);

CCTK_INT GetBoundarySpecification(const CCTK_INT size ,
 CCTK_INT* nboundaryzones ,
 CCTK_INT* is_internal ,
 CCTK_INT* is_staggered ,
 CCTK_INT* shiftout);

CCTK_INT GetDomainSpecification(const CCTK_INT size ,
 CCTK_REAL* physical_min ,
 CCTK_REAL* physical_max ,
 CCTK_REAL* interior_min ,
 CCTK_REAL* interior_max ,
 CCTK_REAL* exterior_min ,
 CCTK_REAL* exterior_max ,
 CCTK_REAL* spacing);

CCTK_INT MultiPatch_ConvertFromPhysicalBoundary(const CCTK_INT map ,
 const CCTK_INT size ,
 const CCTK_REAL* physical_min ,
 const CCTK_REAL* physical_max ,
 CCTK_REAL* interior_min ,
 CCTK_REAL* interior_max ,
 CCTK_REAL* exterior_min ,
 CCTK_REAL* exterior_max ,
 const CCTK_REAL* spacing);

CCTK_INT MultiPatch_GetBoundarySpecification(const CCTK_INT map ,
 const CCTK_INT size ,
 CCTK_INT* nboundaryzones ,
 CCTK_INT* is_internal ,
 CCTK_INT* is_staggered ,
 CCTK_INT* shiftout);

CCTK_INT MultiPatch_GetDomainSpecification(const CCTK_INT map ,
 const CCTK_INT size ,
 CCTK_REAL* physical_min ,
 CCTK_REAL* physical_max ,
 CCTK_REAL* interior_min ,
 CCTK_REAL* interior_max ,
 CCTK_REAL* exterior_min ,
 CCTK_REAL* exterior_max ,
 CCTK_REAL* spacing);

#ifdef __cplusplus
}
#endif
#endif /* CCODE */

#ifdef FCODE
#ifdef F90CODE
#define DECLARE_CARPETREGRID2_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  interface &&\
     CCTK_INT function Carpet_Regrid (cctkGH, superregss, regsss, force) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
       CCTK_POINTER superregss &&\
       CCTK_POINTER regsss &&\
       CCTK_INT force &&\
     end function Carpet_Regrid &&\
  end interface &&\
  interface &&\
     CCTK_INT function Carpet_RegridMaps (cctkGH, superregsss, regssss, force) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
       CCTK_POINTER superregsss &&\
       CCTK_POINTER regssss &&\
       CCTK_INT force &&\
     end function Carpet_RegridMaps &&\
  end interface &&\
  interface &&\
     CCTK_INT function ConvertFromPhysicalBoundary (size, physical_min, physical_max, interior_min, interior_max, exterior_min, exterior_max, spacing) &&\
       implicit none &&\
       CCTK_INT size &&\
       CCTK_REAL physical_min(*) &&\
       CCTK_REAL physical_max(*) &&\
       CCTK_REAL interior_min(*) &&\
       CCTK_REAL interior_max(*) &&\
       CCTK_REAL exterior_min(*) &&\
       CCTK_REAL exterior_max(*) &&\
       CCTK_REAL spacing(*) &&\
     end function ConvertFromPhysicalBoundary &&\
  end interface &&\
  interface &&\
     CCTK_INT function GetBoundarySpecification (size, nboundaryzones, is_internal, is_staggered, shiftout) &&\
       implicit none &&\
       CCTK_INT size &&\
       CCTK_INT nboundaryzones(*) &&\
       CCTK_INT is_internal(*) &&\
       CCTK_INT is_staggered(*) &&\
       CCTK_INT shiftout(*) &&\
     end function GetBoundarySpecification &&\
  end interface &&\
  interface &&\
     CCTK_INT function GetDomainSpecification (size, physical_min, physical_max, interior_min, interior_max, exterior_min, exterior_max, spacing) &&\
       implicit none &&\
       CCTK_INT size &&\
       CCTK_REAL physical_min(*) &&\
       CCTK_REAL physical_max(*) &&\
       CCTK_REAL interior_min(*) &&\
       CCTK_REAL interior_max(*) &&\
       CCTK_REAL exterior_min(*) &&\
       CCTK_REAL exterior_max(*) &&\
       CCTK_REAL spacing(*) &&\
     end function GetDomainSpecification &&\
  end interface &&\
  interface &&\
     CCTK_INT function MultiPatch_ConvertFromPhysicalBoundary (map, size, physical_min, physical_max, interior_min, interior_max, exterior_min, exterior_max, spacing) &&\
       implicit none &&\
       CCTK_INT map &&\
       CCTK_INT size &&\
       CCTK_REAL physical_min(*) &&\
       CCTK_REAL physical_max(*) &&\
       CCTK_REAL interior_min(*) &&\
       CCTK_REAL interior_max(*) &&\
       CCTK_REAL exterior_min(*) &&\
       CCTK_REAL exterior_max(*) &&\
       CCTK_REAL spacing(*) &&\
     end function MultiPatch_ConvertFromPhysicalBoundary &&\
  end interface &&\
  interface &&\
     CCTK_INT function MultiPatch_GetBoundarySpecification (map, size, nboundaryzones, is_internal, is_staggered, shiftout) &&\
       implicit none &&\
       CCTK_INT map &&\
       CCTK_INT size &&\
       CCTK_INT nboundaryzones(*) &&\
       CCTK_INT is_internal(*) &&\
       CCTK_INT is_staggered(*) &&\
       CCTK_INT shiftout(*) &&\
     end function MultiPatch_GetBoundarySpecification &&\
  end interface &&\
  interface &&\
     CCTK_INT function MultiPatch_GetDomainSpecification (map, size, physical_min, physical_max, interior_min, interior_max, exterior_min, exterior_max, spacing) &&\
       implicit none &&\
       CCTK_INT map &&\
       CCTK_INT size &&\
       CCTK_REAL physical_min(*) &&\
       CCTK_REAL physical_max(*) &&\
       CCTK_REAL interior_min(*) &&\
       CCTK_REAL interior_max(*) &&\
       CCTK_REAL exterior_min(*) &&\
       CCTK_REAL exterior_max(*) &&\
       CCTK_REAL spacing(*) &&\
     end function MultiPatch_GetDomainSpecification &&\
  end interface &&\

#else /* ! F90CODE */

#define DECLARE_CARPETREGRID2_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  external Carpet_Regrid &&\
  CCTK_INT Carpet_Regrid &&\
  external Carpet_RegridMaps &&\
  CCTK_INT Carpet_RegridMaps &&\
  external ConvertFromPhysicalBoundary &&\
  CCTK_INT ConvertFromPhysicalBoundary &&\
  external GetBoundarySpecification &&\
  CCTK_INT GetBoundarySpecification &&\
  external GetDomainSpecification &&\
  CCTK_INT GetDomainSpecification &&\
  external MultiPatch_ConvertFromPhysicalBoundary &&\
  CCTK_INT MultiPatch_ConvertFromPhysicalBoundary &&\
  external MultiPatch_GetBoundarySpecification &&\
  CCTK_INT MultiPatch_GetBoundarySpecification &&\
  external MultiPatch_GetDomainSpecification &&\
  CCTK_INT MultiPatch_GetDomainSpecification &&\

#endif /* ! F90CODE */
#endif /* FCODE */

#endif

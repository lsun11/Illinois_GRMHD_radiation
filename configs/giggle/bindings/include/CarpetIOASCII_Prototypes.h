/*@@
   @header  CarpetIOASCII_Prototypes.h
   @author  Automatically generated by CreateFunctionBindings.pl
   @desc
            Prototypes for overloaded functions used by this thorn
   @enddesc
  @@*/


#ifndef _CARPETIOASCII_PROTOTYPES_H_
#define _CARPETIOASCII_PROTOTYPES_H_  1

#ifdef CCODE
#ifdef __cplusplus
extern "C" {
#endif

CCTK_INT IO_TruncateOutputFiles(const CCTK_POINTER_TO_CONST GH);

CCTK_INT SymmetryTableHandleForGI(const CCTK_POINTER_TO_CONST cctkGH ,
 const CCTK_INT group_index);

CCTK_INT SymmetryTableHandleForGrid(const CCTK_POINTER_TO_CONST cctkGH);

CCTK_POINTER_TO_CONST UniqueBuildID(const CCTK_POINTER_TO_CONST cctkGH);

CCTK_POINTER_TO_CONST UniqueRunID(const CCTK_POINTER_TO_CONST cctkGH);

CCTK_POINTER_TO_CONST UniqueSimulationID(const CCTK_POINTER_TO_CONST cctkGH);

#ifdef __cplusplus
}
#endif
#endif /* CCODE */

#ifdef FCODE
#ifdef F90CODE
#define DECLARE_CARPETIOASCII_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  interface &&\
     CCTK_INT function IO_TruncateOutputFiles (GH) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST GH &&\
     end function IO_TruncateOutputFiles &&\
  end interface &&\
  interface &&\
     CCTK_INT function SymmetryTableHandleForGI (cctkGH, group_index) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
       CCTK_INT group_index &&\
     end function SymmetryTableHandleForGI &&\
  end interface &&\
  interface &&\
     CCTK_INT function SymmetryTableHandleForGrid (cctkGH) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
     end function SymmetryTableHandleForGrid &&\
  end interface &&\
  interface &&\
     CCTK_POINTER_TO_CONST function UniqueBuildID (cctkGH) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
     end function UniqueBuildID &&\
  end interface &&\
  interface &&\
     CCTK_POINTER_TO_CONST function UniqueRunID (cctkGH) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
     end function UniqueRunID &&\
  end interface &&\
  interface &&\
     CCTK_POINTER_TO_CONST function UniqueSimulationID (cctkGH) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST cctkGH &&\
     end function UniqueSimulationID &&\
  end interface &&\

#else /* ! F90CODE */

#define DECLARE_CARPETIOASCII_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  external IO_TruncateOutputFiles &&\
  CCTK_INT IO_TruncateOutputFiles &&\
  external SymmetryTableHandleForGI &&\
  CCTK_INT SymmetryTableHandleForGI &&\
  external SymmetryTableHandleForGrid &&\
  CCTK_INT SymmetryTableHandleForGrid &&\
  external UniqueBuildID &&\
  CCTK_POINTER_TO_CONST UniqueBuildID &&\
  external UniqueRunID &&\
  CCTK_POINTER_TO_CONST UniqueRunID &&\
  external UniqueSimulationID &&\
  CCTK_POINTER_TO_CONST UniqueSimulationID &&\

#endif /* ! F90CODE */
#endif /* FCODE */

#endif

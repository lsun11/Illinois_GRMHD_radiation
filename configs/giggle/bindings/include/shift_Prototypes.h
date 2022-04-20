/*@@
   @header  shift_Prototypes.h
   @author  Automatically generated by CreateFunctionBindings.pl
   @desc
            Prototypes for overloaded functions used by this thorn
   @enddesc
  @@*/


#ifndef _SHIFT_PROTOTYPES_H_
#define _SHIFT_PROTOTYPES_H_  1

#ifdef CCODE
#ifdef __cplusplus
extern "C" {
#endif

CCTK_INT Boundary_SelectGroupForBC(const CCTK_POINTER_TO_CONST GH ,
 const CCTK_INT faces ,
 const CCTK_INT boundary_width ,
 const CCTK_INT table_handle ,
 CCTK_STRING var_name, CCTK_STRING bc_name);

CCTK_INT HorizonRadiusInDirection(const CCTK_INT horizon_number ,
 const CCTK_INT N_points ,
 const CCTK_REAL* x ,
 const CCTK_REAL* y ,
 const CCTK_REAL* z ,
 CCTK_REAL* radius);

CCTK_INT MoLRegisterConstrained(const CCTK_INT ConstrainedIndex);

CCTK_INT MoLRegisterEvolvedGroup(const CCTK_INT EvolvedIndex ,
 const CCTK_INT RHSIndex);

#ifdef __cplusplus
}
#endif
#endif /* CCODE */

#ifdef FCODE
#ifdef F90CODE
#define DECLARE_SHIFT_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  interface &&\
     CCTK_INT function Boundary_SelectGroupForBC (GH, faces, boundary_width, table_handle, var_name, bc_name) &&\
       implicit none &&\
       CCTK_POINTER_TO_CONST GH &&\
       CCTK_INT faces &&\
       CCTK_INT boundary_width &&\
       CCTK_INT table_handle &&\
       character(*) var_name &&\
       character(*) bc_name &&\
     end function Boundary_SelectGroupForBC &&\
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
     CCTK_INT function MoLRegisterConstrained (ConstrainedIndex) &&\
       implicit none &&\
       CCTK_INT ConstrainedIndex &&\
     end function MoLRegisterConstrained &&\
  end interface &&\
  interface &&\
     CCTK_INT function MoLRegisterEvolvedGroup (EvolvedIndex, RHSIndex) &&\
       implicit none &&\
       CCTK_INT EvolvedIndex &&\
       CCTK_INT RHSIndex &&\
     end function MoLRegisterEvolvedGroup &&\
  end interface &&\

#else /* ! F90CODE */

#define DECLARE_SHIFT_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \
  external Boundary_SelectGroupForBC &&\
  CCTK_INT Boundary_SelectGroupForBC &&\
  external HorizonRadiusInDirection &&\
  CCTK_INT HorizonRadiusInDirection &&\
  external MoLRegisterConstrained &&\
  CCTK_INT MoLRegisterConstrained &&\
  external MoLRegisterEvolvedGroup &&\
  CCTK_INT MoLRegisterEvolvedGroup &&\

#endif /* ! F90CODE */
#endif /* FCODE */

#endif

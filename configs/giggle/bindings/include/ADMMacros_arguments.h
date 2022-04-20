/*@@
   @header  ADMMacros_arguments.h
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Defines macros to declare/define/pass function arguments
            in calls from C to Fortran for thorn ADMMacros
   @enddesc
 @@*/


#ifdef FCODE
#define DECLARE_ADMMACROS_PRIVATE_FARGUMENTS \


#define ADMMACROS_PRIVATE_FARGUMENTS \


#endif /* FCODE */

#ifdef CCODE
#define DECLARE_ADMMACROS_PRIVATE_CARGUMENTS \


#define DECLARE_ADMMACROS_PRIVATE_C2F \


#define INITIALISE_ADMMACROS_PRIVATE_C2F \


#define ADMMACROS_PRIVATE_C2F_PROTO \


#define PASS_ADMMACROS_PRIVATE_C2F(GH) \


#endif /* CCODE */

#ifdef FCODE
#define DECLARE_ADMMACROS_PROTECTED_FARGUMENTS \


#define ADMMACROS_PROTECTED_FARGUMENTS \


#endif /* FCODE */

#ifdef CCODE
#define DECLARE_ADMMACROS_PROTECTED_CARGUMENTS \


#define DECLARE_ADMMACROS_PROTECTED_C2F \


#define INITIALISE_ADMMACROS_PROTECTED_C2F \


#define ADMMACROS_PROTECTED_C2F_PROTO \


#define PASS_ADMMACROS_PROTECTED_C2F(GH) \


#endif /* CCODE */

#ifdef FCODE
#define DECLARE_ADMMACROS_PUBLIC_FARGUMENTS \
CCTK_INT local_spatial_order&&\


#define ADMMACROS_PUBLIC_FARGUMENTS \
local_spatial_order

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_ADMMACROS_PUBLIC_CARGUMENTS \
CCTK_INT * CCTK_RESTRICT local_spatial_order = (cctki_dummy_int = &local_spatial_order - &local_spatial_order, (CCTK_INT *) CCTKi_VarDataPtr(cctkGH, 0, "ADMMACROS", "local_spatial_order"));

#define DECLARE_ADMMACROS_PUBLIC_C2F \
static int CCTKARGNUM_local_spatial_order = -1; \
static int CCTKGROUPNUM_local_spatial_order = -1;

#define INITIALISE_ADMMACROS_PUBLIC_C2F \
if(CCTKARGNUM_local_spatial_order == -1) CCTKARGNUM_local_spatial_order = CCTK_VarIndex("ADMMacros::local_spatial_order"); \
if(CCTKGROUPNUM_local_spatial_order == -1) CCTKGROUPNUM_local_spatial_order = CCTK_GroupIndex("ADMMacros::local_spatial_order");

#define ADMMACROS_PUBLIC_C2F_PROTO \
CCTK_INT *

#define PASS_ADMMACROS_PUBLIC_C2F(GH) \
(CCTK_INT *)(PASS_REFERENCE(local_spatial_order, 0))

#endif /* CCODE */

#ifdef FCODE
#define ADMMACROS_FARGUMENTS _CCTK_FARGUMENTS, ADMMACROS_PUBLIC_FARGUMENTS

#define DECLARE_ADMMACROS_FARGUMENTS _DECLARE_CCTK_FARGUMENTS DECLARE_ADMMACROS_PUBLIC_FARGUMENTS

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_ADMMACROS_CARGUMENTS _DECLARE_CCTK_CARGUMENTS DECLARE_ADMMACROS_PUBLIC_CARGUMENTS

#define ADMMACROS_C2F_PROTO _CCTK_C2F_PROTO, ADMMACROS_PUBLIC_C2F_PROTO

#define PASS_ADMMACROS_C2F(GH) _PASS_CCTK_C2F(GH), PASS_ADMMACROS_PUBLIC_C2F(GH)

#define DECLARE_ADMMACROS_C2F _DECLARE_CCTK_C2F DECLARE_ADMMACROS_PUBLIC_C2F

#define INITIALISE_ADMMACROS_C2F _INITIALISE_CCTK_C2F INITIALISE_ADMMACROS_PUBLIC_C2F

#define ADMMACROS_CARGUMENTS cGH *cctkGH

#endif /* CCODE */
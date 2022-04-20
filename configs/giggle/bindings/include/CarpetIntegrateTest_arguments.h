/*@@
   @header  CarpetIntegrateTest_arguments.h
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Defines macros to declare/define/pass function arguments
            in calls from C to Fortran for thorn CarpetIntegrateTest
   @enddesc
 @@*/


#ifdef FCODE
#define DECLARE_CARPETINTEGRATETEST_PRIVATE_FARGUMENTS \
INTEGER X0integrand&&\
INTEGER X1integrand&&\
INTEGER X2integrand&&\
CCTK_REAL integrand(X0integrand,X1integrand,X2integrand)&&\
CCTK_REAL integrand_p(X0integrand,X1integrand,X2integrand)&&\
CCTK_REAL integrand_p_p(X0integrand,X1integrand,X2integrand)&&\


#define CARPETINTEGRATETEST_PRIVATE_FARGUMENTS \
X0integrand,X1integrand,X2integrand,integrand,integrand_p,integrand_p_p

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_CARPETINTEGRATETEST_PRIVATE_CARGUMENTS \
CCTK_REAL * CCTK_RESTRICT integrand = (cctki_dummy_int = &integrand - &integrand, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "CARPETINTEGRATETEST", "integrand")); \
CCTK_REAL * CCTK_RESTRICT integrand_p = (cctki_dummy_int = &integrand_p - &integrand_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 1, "CARPETINTEGRATETEST", "integrand")); \
CCTK_REAL * CCTK_RESTRICT integrand_p_p = (cctki_dummy_int = &integrand_p_p - &integrand_p_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 2, "CARPETINTEGRATETEST", "integrand"));

#define DECLARE_CARPETINTEGRATETEST_PRIVATE_C2F \
static int CCTKARGNUM_integrand = -1; \
static int CCTKGROUPNUM_integrand = -1;

#define INITIALISE_CARPETINTEGRATETEST_PRIVATE_C2F \
if(CCTKARGNUM_integrand == -1) CCTKARGNUM_integrand = CCTK_VarIndex("CarpetIntegrateTest::integrand"); \
if(CCTKGROUPNUM_integrand == -1) CCTKGROUPNUM_integrand = CCTK_GroupIndex("CarpetIntegrateTest::integrand");

#define CARPETINTEGRATETEST_PRIVATE_C2F_PROTO \
const int *,const int *,const int *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *

#define PASS_CARPETINTEGRATETEST_PRIVATE_C2F(GH) \
PASS_GROUPSIZE(integrand, 0),\
PASS_GROUPSIZE(integrand, 1),\
PASS_GROUPSIZE(integrand, 2),\
(CCTK_REAL *)(PASS_REFERENCE(integrand, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(integrand, 1)),\
(CCTK_REAL *)(PASS_REFERENCE(integrand, 2))

#endif /* CCODE */

#ifdef FCODE
#define DECLARE_CARPETINTEGRATETEST_PROTECTED_FARGUMENTS \


#define CARPETINTEGRATETEST_PROTECTED_FARGUMENTS \


#endif /* FCODE */

#ifdef CCODE
#define DECLARE_CARPETINTEGRATETEST_PROTECTED_CARGUMENTS \


#define DECLARE_CARPETINTEGRATETEST_PROTECTED_C2F \


#define INITIALISE_CARPETINTEGRATETEST_PROTECTED_C2F \


#define CARPETINTEGRATETEST_PROTECTED_C2F_PROTO \


#define PASS_CARPETINTEGRATETEST_PROTECTED_C2F(GH) \


#endif /* CCODE */

#ifdef FCODE
#define DECLARE_CARPETINTEGRATETEST_PUBLIC_FARGUMENTS \
INTEGER X0coordinates&&\
INTEGER X1coordinates&&\
INTEGER X2coordinates&&\
CCTK_REAL coarse_dx&&\
CCTK_REAL coarse_dy&&\
CCTK_REAL coarse_dz&&\
CCTK_REAL r(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL x(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL y(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL z(X0coordinates,X1coordinates,X2coordinates)&&\


#define CARPETINTEGRATETEST_PUBLIC_FARGUMENTS \
X0coordinates,X1coordinates,X2coordinates,coarse_dx,coarse_dy,coarse_dz,r,x,y,z

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_CARPETINTEGRATETEST_PUBLIC_CARGUMENTS \
CCTK_REAL * CCTK_RESTRICT coarse_dx = (cctki_dummy_int = &coarse_dx - &coarse_dx, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dx")); \
CCTK_REAL * CCTK_RESTRICT coarse_dy = (cctki_dummy_int = &coarse_dy - &coarse_dy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dy")); \
CCTK_REAL * CCTK_RESTRICT coarse_dz = (cctki_dummy_int = &coarse_dz - &coarse_dz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dz")); \
CCTK_REAL * CCTK_RESTRICT r = (cctki_dummy_int = &r - &r, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "r")); \
CCTK_REAL * CCTK_RESTRICT x = (cctki_dummy_int = &x - &x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "x")); \
CCTK_REAL * CCTK_RESTRICT y = (cctki_dummy_int = &y - &y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "y")); \
CCTK_REAL * CCTK_RESTRICT z = (cctki_dummy_int = &z - &z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "z"));

#define DECLARE_CARPETINTEGRATETEST_PUBLIC_C2F \
static int CCTKARGNUM_coarse_dx = -1; \
static int CCTKGROUPNUM_gridspacings = -1; \
static int CCTKARGNUM_coarse_dy = -1; \
static int CCTKARGNUM_coarse_dz = -1; \
static int CCTKARGNUM_r = -1; \
static int CCTKGROUPNUM_coordinates = -1; \
static int CCTKARGNUM_x = -1; \
static int CCTKARGNUM_y = -1; \
static int CCTKARGNUM_z = -1;

#define INITIALISE_CARPETINTEGRATETEST_PUBLIC_C2F \
if(CCTKARGNUM_coarse_dx == -1) CCTKARGNUM_coarse_dx = CCTK_VarIndex("GRID::coarse_dx"); \
if(CCTKGROUPNUM_gridspacings == -1) CCTKGROUPNUM_gridspacings = CCTK_GroupIndex("GRID::gridspacings"); \
if(CCTKARGNUM_coarse_dy == -1) CCTKARGNUM_coarse_dy = CCTK_VarIndex("GRID::coarse_dy"); \
if(CCTKARGNUM_coarse_dz == -1) CCTKARGNUM_coarse_dz = CCTK_VarIndex("GRID::coarse_dz"); \
if(CCTKARGNUM_r == -1) CCTKARGNUM_r = CCTK_VarIndex("GRID::r"); \
if(CCTKGROUPNUM_coordinates == -1) CCTKGROUPNUM_coordinates = CCTK_GroupIndex("GRID::coordinates"); \
if(CCTKARGNUM_x == -1) CCTKARGNUM_x = CCTK_VarIndex("GRID::x"); \
if(CCTKARGNUM_y == -1) CCTKARGNUM_y = CCTK_VarIndex("GRID::y"); \
if(CCTKARGNUM_z == -1) CCTKARGNUM_z = CCTK_VarIndex("GRID::z");

#define CARPETINTEGRATETEST_PUBLIC_C2F_PROTO \
const int *,const int *,const int *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *

#define PASS_CARPETINTEGRATETEST_PUBLIC_C2F(GH) \
PASS_GROUPSIZE(coordinates, 0),\
PASS_GROUPSIZE(coordinates, 1),\
PASS_GROUPSIZE(coordinates, 2),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dx, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(r, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(z, 0))

#endif /* CCODE */

#ifdef FCODE
#define CARPETINTEGRATETEST_FARGUMENTS _CCTK_FARGUMENTS, CARPETINTEGRATETEST_PRIVATE_FARGUMENTS, CARPETINTEGRATETEST_PUBLIC_FARGUMENTS

#define DECLARE_CARPETINTEGRATETEST_FARGUMENTS _DECLARE_CCTK_FARGUMENTS DECLARE_CARPETINTEGRATETEST_PRIVATE_FARGUMENTS DECLARE_CARPETINTEGRATETEST_PUBLIC_FARGUMENTS

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_CARPETINTEGRATETEST_CARGUMENTS _DECLARE_CCTK_CARGUMENTS DECLARE_CARPETINTEGRATETEST_PRIVATE_CARGUMENTS DECLARE_CARPETINTEGRATETEST_PUBLIC_CARGUMENTS

#define CARPETINTEGRATETEST_C2F_PROTO _CCTK_C2F_PROTO, CARPETINTEGRATETEST_PRIVATE_C2F_PROTO, CARPETINTEGRATETEST_PUBLIC_C2F_PROTO

#define PASS_CARPETINTEGRATETEST_C2F(GH) _PASS_CCTK_C2F(GH), PASS_CARPETINTEGRATETEST_PRIVATE_C2F(GH), PASS_CARPETINTEGRATETEST_PUBLIC_C2F(GH)

#define DECLARE_CARPETINTEGRATETEST_C2F _DECLARE_CCTK_C2F DECLARE_CARPETINTEGRATETEST_PRIVATE_C2F DECLARE_CARPETINTEGRATETEST_PUBLIC_C2F

#define INITIALISE_CARPETINTEGRATETEST_C2F _INITIALISE_CCTK_C2F INITIALISE_CARPETINTEGRATETEST_PRIVATE_C2F INITIALISE_CARPETINTEGRATETEST_PUBLIC_C2F

#define CARPETINTEGRATETEST_CARGUMENTS cGH *cctkGH

#endif /* CCODE */

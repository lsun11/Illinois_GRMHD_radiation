/*@@
   @header  IDWaveMoL_arguments.h
   @author  Automatically generated by GridFuncStuff.pl
   @desc
            Defines macros to declare/define/pass function arguments
            in calls from C to Fortran for thorn IDWaveMoL
   @enddesc
 @@*/


#ifdef FCODE
#define DECLARE_IDWAVEMOL_PRIVATE_FARGUMENTS \


#define IDWAVEMOL_PRIVATE_FARGUMENTS \


#endif /* FCODE */

#ifdef CCODE
#define DECLARE_IDWAVEMOL_PRIVATE_CARGUMENTS \


#define DECLARE_IDWAVEMOL_PRIVATE_C2F \


#define INITIALISE_IDWAVEMOL_PRIVATE_C2F \


#define IDWAVEMOL_PRIVATE_C2F_PROTO \


#define PASS_IDWAVEMOL_PRIVATE_C2F(GH) \


#endif /* CCODE */

#ifdef FCODE
#define DECLARE_IDWAVEMOL_PROTECTED_FARGUMENTS \


#define IDWAVEMOL_PROTECTED_FARGUMENTS \


#endif /* FCODE */

#ifdef CCODE
#define DECLARE_IDWAVEMOL_PROTECTED_CARGUMENTS \


#define DECLARE_IDWAVEMOL_PROTECTED_C2F \


#define INITIALISE_IDWAVEMOL_PROTECTED_C2F \


#define IDWAVEMOL_PROTECTED_C2F_PROTO \


#define PASS_IDWAVEMOL_PROTECTED_C2F(GH) \


#endif /* CCODE */

#ifdef FCODE
#define DECLARE_IDWAVEMOL_PUBLIC_FARGUMENTS \
INTEGER X0coordinates&&\
INTEGER X0energy&&\
INTEGER X0scalarevolvemol_scalar&&\
INTEGER X0scalarevolvemol_vector&&\
INTEGER X0scalarrhsmol_scalar&&\
INTEGER X0scalarrhsmol_vector&&\
INTEGER X1coordinates&&\
INTEGER X1energy&&\
INTEGER X1scalarevolvemol_scalar&&\
INTEGER X1scalarevolvemol_vector&&\
INTEGER X1scalarrhsmol_scalar&&\
INTEGER X1scalarrhsmol_vector&&\
INTEGER X2coordinates&&\
INTEGER X2energy&&\
INTEGER X2scalarevolvemol_scalar&&\
INTEGER X2scalarevolvemol_vector&&\
INTEGER X2scalarrhsmol_scalar&&\
INTEGER X2scalarrhsmol_vector&&\
CCTK_REAL coarse_dx&&\
CCTK_REAL coarse_dy&&\
CCTK_REAL coarse_dz&&\
CCTK_REAL energy(X0energy,X1energy,X2energy)&&\
CCTK_REAL phi(X0scalarevolvemol_scalar,X1scalarevolvemol_scalar,X2scalarevolvemol_scalar)&&\
CCTK_REAL phi_p(X0scalarevolvemol_scalar,X1scalarevolvemol_scalar,X2scalarevolvemol_scalar)&&\
CCTK_REAL phi_p_p(X0scalarevolvemol_scalar,X1scalarevolvemol_scalar,X2scalarevolvemol_scalar)&&\
CCTK_REAL phirhs(X0scalarrhsmol_scalar,X1scalarrhsmol_scalar,X2scalarrhsmol_scalar)&&\
CCTK_REAL phit(X0scalarevolvemol_scalar,X1scalarevolvemol_scalar,X2scalarevolvemol_scalar)&&\
CCTK_REAL phit_p(X0scalarevolvemol_scalar,X1scalarevolvemol_scalar,X2scalarevolvemol_scalar)&&\
CCTK_REAL phit_p_p(X0scalarevolvemol_scalar,X1scalarevolvemol_scalar,X2scalarevolvemol_scalar)&&\
CCTK_REAL phitrhs(X0scalarrhsmol_scalar,X1scalarrhsmol_scalar,X2scalarrhsmol_scalar)&&\
CCTK_REAL phix(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phix_p(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phix_p_p(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phixrhs(X0scalarrhsmol_vector,X1scalarrhsmol_vector,X2scalarrhsmol_vector)&&\
CCTK_REAL phiy(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phiy_p(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phiy_p_p(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phiyrhs(X0scalarrhsmol_vector,X1scalarrhsmol_vector,X2scalarrhsmol_vector)&&\
CCTK_REAL phiz(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phiz_p(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phiz_p_p(X0scalarevolvemol_vector,X1scalarevolvemol_vector,X2scalarevolvemol_vector)&&\
CCTK_REAL phizrhs(X0scalarrhsmol_vector,X1scalarrhsmol_vector,X2scalarrhsmol_vector)&&\
CCTK_REAL r(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL x(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL y(X0coordinates,X1coordinates,X2coordinates)&&\
CCTK_REAL z(X0coordinates,X1coordinates,X2coordinates)&&\


#define IDWAVEMOL_PUBLIC_FARGUMENTS \
X0coordinates,X0energy,X0scalarevolvemol_scalar,X0scalarevolvemol_vector,X0scalarrhsmol_scalar,X0scalarrhsmol_vector,X1coordinates,X1energy,X1scalarevolvemol_scalar,X1scalarevolvemol_vector,X1scalarrhsmol_scalar,X1scalarrhsmol_vector,X2coordinates,X2energy,X2scalarevolvemol_scalar,X2scalarevolvemol_vector,X2scalarrhsmol_scalar,X2scalarrhsmol_vector,coarse_dx,coarse_dy,coarse_dz,energy,phi,phi_p,phi_p_p,phirhs,phit,phit_p,phit_p_p,phitrhs,phix,phix_p,phix_p_p,phixrhs,phiy,phiy_p,phiy_p_p,phiyrhs,phiz,phiz_p,phiz_p_p,phizrhs,r,x,y,z

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_IDWAVEMOL_PUBLIC_CARGUMENTS \
CCTK_REAL * CCTK_RESTRICT coarse_dx = (cctki_dummy_int = &coarse_dx - &coarse_dx, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dx")); \
CCTK_REAL * CCTK_RESTRICT coarse_dy = (cctki_dummy_int = &coarse_dy - &coarse_dy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dy")); \
CCTK_REAL * CCTK_RESTRICT coarse_dz = (cctki_dummy_int = &coarse_dz - &coarse_dz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "coarse_dz")); \
CCTK_REAL * CCTK_RESTRICT energy = (cctki_dummy_int = &energy - &energy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "energy")); \
CCTK_REAL * CCTK_RESTRICT phi = (cctki_dummy_int = &phi - &phi, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phi")); \
CCTK_REAL * CCTK_RESTRICT phi_p = (cctki_dummy_int = &phi_p - &phi_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 1, "WAVEMOL", "phi")); \
CCTK_REAL * CCTK_RESTRICT phi_p_p = (cctki_dummy_int = &phi_p_p - &phi_p_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 2, "WAVEMOL", "phi")); \
CCTK_REAL * CCTK_RESTRICT phirhs = (cctki_dummy_int = &phirhs - &phirhs, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phirhs")); \
CCTK_REAL * CCTK_RESTRICT phit = (cctki_dummy_int = &phit - &phit, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phit")); \
CCTK_REAL * CCTK_RESTRICT phit_p = (cctki_dummy_int = &phit_p - &phit_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 1, "WAVEMOL", "phit")); \
CCTK_REAL * CCTK_RESTRICT phit_p_p = (cctki_dummy_int = &phit_p_p - &phit_p_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 2, "WAVEMOL", "phit")); \
CCTK_REAL * CCTK_RESTRICT phitrhs = (cctki_dummy_int = &phitrhs - &phitrhs, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phitrhs")); \
CCTK_REAL * CCTK_RESTRICT phix = (cctki_dummy_int = &phix - &phix, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phix")); \
CCTK_REAL * CCTK_RESTRICT phix_p = (cctki_dummy_int = &phix_p - &phix_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 1, "WAVEMOL", "phix")); \
CCTK_REAL * CCTK_RESTRICT phix_p_p = (cctki_dummy_int = &phix_p_p - &phix_p_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 2, "WAVEMOL", "phix")); \
CCTK_REAL * CCTK_RESTRICT phixrhs = (cctki_dummy_int = &phixrhs - &phixrhs, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phixrhs")); \
CCTK_REAL * CCTK_RESTRICT phiy = (cctki_dummy_int = &phiy - &phiy, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phiy")); \
CCTK_REAL * CCTK_RESTRICT phiy_p = (cctki_dummy_int = &phiy_p - &phiy_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 1, "WAVEMOL", "phiy")); \
CCTK_REAL * CCTK_RESTRICT phiy_p_p = (cctki_dummy_int = &phiy_p_p - &phiy_p_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 2, "WAVEMOL", "phiy")); \
CCTK_REAL * CCTK_RESTRICT phiyrhs = (cctki_dummy_int = &phiyrhs - &phiyrhs, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phiyrhs")); \
CCTK_REAL * CCTK_RESTRICT phiz = (cctki_dummy_int = &phiz - &phiz, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phiz")); \
CCTK_REAL * CCTK_RESTRICT phiz_p = (cctki_dummy_int = &phiz_p - &phiz_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 1, "WAVEMOL", "phiz")); \
CCTK_REAL * CCTK_RESTRICT phiz_p_p = (cctki_dummy_int = &phiz_p_p - &phiz_p_p, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 2, "WAVEMOL", "phiz")); \
CCTK_REAL * CCTK_RESTRICT phizrhs = (cctki_dummy_int = &phizrhs - &phizrhs, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "WAVEMOL", "phizrhs")); \
CCTK_REAL * CCTK_RESTRICT r = (cctki_dummy_int = &r - &r, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "r")); \
CCTK_REAL * CCTK_RESTRICT x = (cctki_dummy_int = &x - &x, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "x")); \
CCTK_REAL * CCTK_RESTRICT y = (cctki_dummy_int = &y - &y, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "y")); \
CCTK_REAL * CCTK_RESTRICT z = (cctki_dummy_int = &z - &z, (CCTK_REAL *) CCTKi_VarDataPtr(cctkGH, 0, "GRID", "z"));

#define DECLARE_IDWAVEMOL_PUBLIC_C2F \
static int CCTKARGNUM_coarse_dx = -1; \
static int CCTKGROUPNUM_gridspacings = -1; \
static int CCTKARGNUM_coarse_dy = -1; \
static int CCTKARGNUM_coarse_dz = -1; \
static int CCTKARGNUM_energy = -1; \
static int CCTKGROUPNUM_energy = -1; \
static int CCTKARGNUM_phi = -1; \
static int CCTKGROUPNUM_scalarevolvemol_scalar = -1; \
static int CCTKARGNUM_phirhs = -1; \
static int CCTKGROUPNUM_scalarrhsmol_scalar = -1; \
static int CCTKARGNUM_phit = -1; \
static int CCTKARGNUM_phitrhs = -1; \
static int CCTKARGNUM_phix = -1; \
static int CCTKGROUPNUM_scalarevolvemol_vector = -1; \
static int CCTKARGNUM_phixrhs = -1; \
static int CCTKGROUPNUM_scalarrhsmol_vector = -1; \
static int CCTKARGNUM_phiy = -1; \
static int CCTKARGNUM_phiyrhs = -1; \
static int CCTKARGNUM_phiz = -1; \
static int CCTKARGNUM_phizrhs = -1; \
static int CCTKARGNUM_r = -1; \
static int CCTKGROUPNUM_coordinates = -1; \
static int CCTKARGNUM_x = -1; \
static int CCTKARGNUM_y = -1; \
static int CCTKARGNUM_z = -1;

#define INITIALISE_IDWAVEMOL_PUBLIC_C2F \
if(CCTKARGNUM_coarse_dx == -1) CCTKARGNUM_coarse_dx = CCTK_VarIndex("GRID::coarse_dx"); \
if(CCTKGROUPNUM_gridspacings == -1) CCTKGROUPNUM_gridspacings = CCTK_GroupIndex("GRID::gridspacings"); \
if(CCTKARGNUM_coarse_dy == -1) CCTKARGNUM_coarse_dy = CCTK_VarIndex("GRID::coarse_dy"); \
if(CCTKARGNUM_coarse_dz == -1) CCTKARGNUM_coarse_dz = CCTK_VarIndex("GRID::coarse_dz"); \
if(CCTKARGNUM_energy == -1) CCTKARGNUM_energy = CCTK_VarIndex("WAVEMOL::energy"); \
if(CCTKGROUPNUM_energy == -1) CCTKGROUPNUM_energy = CCTK_GroupIndex("WAVEMOL::energy"); \
if(CCTKARGNUM_phi == -1) CCTKARGNUM_phi = CCTK_VarIndex("WAVEMOL::phi"); \
if(CCTKGROUPNUM_scalarevolvemol_scalar == -1) CCTKGROUPNUM_scalarevolvemol_scalar = CCTK_GroupIndex("WAVEMOL::scalarevolvemol_scalar"); \
if(CCTKARGNUM_phirhs == -1) CCTKARGNUM_phirhs = CCTK_VarIndex("WAVEMOL::phirhs"); \
if(CCTKGROUPNUM_scalarrhsmol_scalar == -1) CCTKGROUPNUM_scalarrhsmol_scalar = CCTK_GroupIndex("WAVEMOL::scalarrhsmol_scalar"); \
if(CCTKARGNUM_phit == -1) CCTKARGNUM_phit = CCTK_VarIndex("WAVEMOL::phit"); \
if(CCTKARGNUM_phitrhs == -1) CCTKARGNUM_phitrhs = CCTK_VarIndex("WAVEMOL::phitrhs"); \
if(CCTKARGNUM_phix == -1) CCTKARGNUM_phix = CCTK_VarIndex("WAVEMOL::phix"); \
if(CCTKGROUPNUM_scalarevolvemol_vector == -1) CCTKGROUPNUM_scalarevolvemol_vector = CCTK_GroupIndex("WAVEMOL::scalarevolvemol_vector"); \
if(CCTKARGNUM_phixrhs == -1) CCTKARGNUM_phixrhs = CCTK_VarIndex("WAVEMOL::phixrhs"); \
if(CCTKGROUPNUM_scalarrhsmol_vector == -1) CCTKGROUPNUM_scalarrhsmol_vector = CCTK_GroupIndex("WAVEMOL::scalarrhsmol_vector"); \
if(CCTKARGNUM_phiy == -1) CCTKARGNUM_phiy = CCTK_VarIndex("WAVEMOL::phiy"); \
if(CCTKARGNUM_phiyrhs == -1) CCTKARGNUM_phiyrhs = CCTK_VarIndex("WAVEMOL::phiyrhs"); \
if(CCTKARGNUM_phiz == -1) CCTKARGNUM_phiz = CCTK_VarIndex("WAVEMOL::phiz"); \
if(CCTKARGNUM_phizrhs == -1) CCTKARGNUM_phizrhs = CCTK_VarIndex("WAVEMOL::phizrhs"); \
if(CCTKARGNUM_r == -1) CCTKARGNUM_r = CCTK_VarIndex("GRID::r"); \
if(CCTKGROUPNUM_coordinates == -1) CCTKGROUPNUM_coordinates = CCTK_GroupIndex("GRID::coordinates"); \
if(CCTKARGNUM_x == -1) CCTKARGNUM_x = CCTK_VarIndex("GRID::x"); \
if(CCTKARGNUM_y == -1) CCTKARGNUM_y = CCTK_VarIndex("GRID::y"); \
if(CCTKARGNUM_z == -1) CCTKARGNUM_z = CCTK_VarIndex("GRID::z");

#define IDWAVEMOL_PUBLIC_C2F_PROTO \
const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,const int *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *,CCTK_REAL *

#define PASS_IDWAVEMOL_PUBLIC_C2F(GH) \
PASS_GROUPSIZE(coordinates, 0),\
PASS_GROUPSIZE(energy, 0),\
PASS_GROUPSIZE(scalarevolvemol_scalar, 0),\
PASS_GROUPSIZE(scalarevolvemol_vector, 0),\
PASS_GROUPSIZE(scalarrhsmol_scalar, 0),\
PASS_GROUPSIZE(scalarrhsmol_vector, 0),\
PASS_GROUPSIZE(coordinates, 1),\
PASS_GROUPSIZE(energy, 1),\
PASS_GROUPSIZE(scalarevolvemol_scalar, 1),\
PASS_GROUPSIZE(scalarevolvemol_vector, 1),\
PASS_GROUPSIZE(scalarrhsmol_scalar, 1),\
PASS_GROUPSIZE(scalarrhsmol_vector, 1),\
PASS_GROUPSIZE(coordinates, 2),\
PASS_GROUPSIZE(energy, 2),\
PASS_GROUPSIZE(scalarevolvemol_scalar, 2),\
PASS_GROUPSIZE(scalarevolvemol_vector, 2),\
PASS_GROUPSIZE(scalarrhsmol_scalar, 2),\
PASS_GROUPSIZE(scalarrhsmol_vector, 2),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dx, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(coarse_dz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(energy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phi, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phi, 1)),\
(CCTK_REAL *)(PASS_REFERENCE(phi, 2)),\
(CCTK_REAL *)(PASS_REFERENCE(phirhs, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phit, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phit, 1)),\
(CCTK_REAL *)(PASS_REFERENCE(phit, 2)),\
(CCTK_REAL *)(PASS_REFERENCE(phitrhs, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phix, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phix, 1)),\
(CCTK_REAL *)(PASS_REFERENCE(phix, 2)),\
(CCTK_REAL *)(PASS_REFERENCE(phixrhs, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phiy, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phiy, 1)),\
(CCTK_REAL *)(PASS_REFERENCE(phiy, 2)),\
(CCTK_REAL *)(PASS_REFERENCE(phiyrhs, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phiz, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(phiz, 1)),\
(CCTK_REAL *)(PASS_REFERENCE(phiz, 2)),\
(CCTK_REAL *)(PASS_REFERENCE(phizrhs, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(r, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(x, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(y, 0)),\
(CCTK_REAL *)(PASS_REFERENCE(z, 0))

#endif /* CCODE */

#ifdef FCODE
#define IDWAVEMOL_FARGUMENTS _CCTK_FARGUMENTS, IDWAVEMOL_PUBLIC_FARGUMENTS

#define DECLARE_IDWAVEMOL_FARGUMENTS _DECLARE_CCTK_FARGUMENTS DECLARE_IDWAVEMOL_PUBLIC_FARGUMENTS

#endif /* FCODE */

#ifdef CCODE
#define DECLARE_IDWAVEMOL_CARGUMENTS _DECLARE_CCTK_CARGUMENTS DECLARE_IDWAVEMOL_PUBLIC_CARGUMENTS

#define IDWAVEMOL_C2F_PROTO _CCTK_C2F_PROTO, IDWAVEMOL_PUBLIC_C2F_PROTO

#define PASS_IDWAVEMOL_C2F(GH) _PASS_CCTK_C2F(GH), PASS_IDWAVEMOL_PUBLIC_C2F(GH)

#define DECLARE_IDWAVEMOL_C2F _DECLARE_CCTK_C2F DECLARE_IDWAVEMOL_PUBLIC_C2F

#define INITIALISE_IDWAVEMOL_C2F _INITIALISE_CCTK_C2F INITIALISE_IDWAVEMOL_PUBLIC_C2F

#define IDWAVEMOL_CARGUMENTS cGH *cctkGH

#endif /* CCODE */

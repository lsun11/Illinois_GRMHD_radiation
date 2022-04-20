//----------------------------------------------------------------------------------
// Advection routines for hydro variables, valid for arbitrary numbers of ghostzones
//----------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define SCALAR 1
#define VECTORX 2
#define VECTORY 3
#define VECTORZ 4
#define BVECTORX 5
#define BVECTORY 6
#define BVECTORZ 7

//FORTRAN-callable headers:
extern "C" void CCTK_FCALL CCTK_FNAME(set_U_syms)
  (int *Symmetry,int *Utype,double *Sym_Bz,int *syms);

extern "C" void set_U_syms(int Symmetry,int Utype,double Sym_Bz,int *syms) {
  //int NO_SYMM = 0; // <- Unused, so commented out.
  int EQUATORIAL = 1;
  int OCTANT = 2;
  int AXISYM = 4;

  syms[0] = 1;
  syms[1] = 1;
  syms[2] = 1;

  if(Symmetry==EQUATORIAL) {
    if(Utype==SCALAR) {
      // do nothing; syms[0]=syms[1]=syms[2]=1 is fine.
    } else if(Utype==VECTORX) {
      // do nothing; syms[0]=syms[1]=syms[2]=1 is fine.
    } else if(Utype==VECTORY) {
      // do nothing; syms[0]=syms[1]=syms[2]=1 is fine.
    } else if(Utype==VECTORZ) {
      syms[2] = -1;
    } else if(Utype==BVECTORX) {
      syms[2] = -Sym_Bz;
    } else if(Utype==BVECTORY) {
      syms[2] = -Sym_Bz;
    } else if(Utype==BVECTORZ) {
      syms[2] = Sym_Bz;
    } else printf("Unsupported symmetry type.  Expected a value between 1-7!\n");
  } else if(Symmetry==AXISYM) {
    if(Utype==SCALAR) {
      // do nothing; syms[0]=syms[1]=syms[2]=1 is fine.
    } else if(Utype==VECTORX) {
      syms[0] = -1;
    } else if(Utype==VECTORY) {
      syms[0] = -1;
    } else if(Utype==VECTORZ) {
      syms[2] = -1;
    } else if(Utype==BVECTORX) {
      syms[0] = -1;
      syms[2] = -Sym_Bz;
    } else if(Utype==BVECTORY) {
      syms[0] = -1;
      syms[2] = -Sym_Bz;
    } else if(Utype==BVECTORZ) {
      syms[2] = Sym_Bz;
    } else printf("Unsupported symmetry type.  Expected a value between 1-7!\n");
  } else if(Symmetry==OCTANT) {
    if(Utype==SCALAR) {
      // do nothing; syms[0]=syms[1]=syms[2]=1 is fine.
    } else if(Utype==VECTORX) {
      syms[0] = -1;
    } else if(Utype==VECTORY) {
      syms[1] = -1;
    } else if(Utype==VECTORZ) {
      syms[2] = -1;
    } else if(Utype==BVECTORX) {
      syms[0] = -1;
      syms[2] = -Sym_Bz;
    } else if(Utype==BVECTORY) {
      syms[1] = -1;
      syms[2] = -Sym_Bz;
    } else if(Utype==BVECTORZ) {
      syms[2] = Sym_Bz;
    } else printf("Unsupported symmetry type.  Expected a value between 1-7!\n");
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(set_U_syms)
  (int *Symmetry,int *Utype,double *Sym_Bz,int *syms)
{
  set_U_syms(*Symmetry,*Utype,*Sym_Bz,syms);
}

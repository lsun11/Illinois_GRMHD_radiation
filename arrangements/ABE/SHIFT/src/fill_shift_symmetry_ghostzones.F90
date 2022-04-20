!-------------------------------------------------------
! Driver routine for updating shift symmetry ghostzones
!-------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine fill_shift_symmetry_gz(ext,X,Y,Z,Symmetry,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt)
  implicit none

  integer, dimension(3)                  :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,shiftx,shifty,shiftz,shiftxt,shiftyt,shiftzt
  integer                                :: Symmetry
  integer, parameter                     :: AXISYM = 4
  !
  gridfunc_syms(1)=-1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,shiftx,gridfunc_syms)
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,shiftxt,gridfunc_syms)
  gridfunc_syms(1)= 1;gridfunc_syms(2)=-1;gridfunc_syms(3)= 1
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,shifty,gridfunc_syms)
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,shiftyt,gridfunc_syms)
  gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)=-1
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,shiftz,gridfunc_syms)
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,shiftzt,gridfunc_syms)

end subroutine fill_shift_symmetry_gz

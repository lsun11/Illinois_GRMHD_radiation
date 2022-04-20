!-------------------------------------------------------
! Driver routine for updating lapse symmetry ghostzones
!-------------------------------------------------------
subroutine fill_lapse_symmetry_gz(ext,X,Y,Z,Symmetry,lapm1)
  implicit none
  integer, dimension(3)                  :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,lapm1
  integer                                :: Symmetry
  integer, parameter                     :: AXISYM = 4
  !
  gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,lapm1,gridfunc_syms)
end subroutine fill_lapse_symmetry_gz
subroutine fill_lapsederivs_symmetry_gz(ext,X,Y,Z,Symmetry,lapsex,lapsey,lapsez)
  implicit none
  integer, dimension(3)                  :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,lapsex,lapsey,lapsez
  integer                                :: Symmetry
  integer, parameter                     :: AXISYM = 4
  !
  gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,lapsex,gridfunc_syms)
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,lapsey,gridfunc_syms)
  call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,lapsez,gridfunc_syms)
end subroutine fill_lapsederivs_symmetry_gz

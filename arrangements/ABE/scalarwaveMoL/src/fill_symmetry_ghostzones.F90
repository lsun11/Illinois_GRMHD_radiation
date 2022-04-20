!------------------------------------------------------------------------------------------------------
! Routines for filling symmetry ghostzones, including master function: scalarwave_do_symmetry_gz_one_gridfunction
!------------------------------------------------------------------------------------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine fill_scalarwave_symmetry_gz(ext,X,Y,Z,Symmetry,phi,phidot)
  implicit none

  integer, dimension(3)                  :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,phi,phidot
  integer                                :: Symmetry
  integer, parameter                     :: AXISYM = 4
  !
  gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
  call scalarwave_do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,phi,gridfunc_syms)
  gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
  call scalarwave_do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,phidot,gridfunc_syms)

end subroutine fill_scalarwave_symmetry_gz

subroutine scalarwave_do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gridfunc,gridfunc_syms)
  implicit none

  integer, dimension(3)                    :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3))   :: X,Y,Z,gridfunc
  real*8                                   :: dX,dY,dZ
  integer                                  :: Symmetry
  integer :: i

  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  if(Symmetry.eq.1) then
     if(abs(z(1,1,2)).lt.dz*0.1) then
        gridfunc(:,:,1) = gridfunc_syms(3)*gridfunc(:,:,3)
     else if(abs(z(1,1,3)).lt.dz*0.1) then
        do i=1,2
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,6-i)
        end do
     else if(abs(z(1,1,4)).lt.dz*0.1) then
        do i=1,3
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,8-i)
        end do
     else if(abs(z(1,1,5)).lt.dz*0.1) then
        do i=1,4
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,10-i)
        end do
     else if(abs(z(1,1,6)).lt.dz*0.1) then
        do i=1,5
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,12-i)
        end do
        !=============================================
     else if(abs(z(1,1,2)-dZ*0.5D0).lt.dz*0.1) then
        gridfunc(:,:,1) = gridfunc_syms(3)*gridfunc(:,:,2)
     else if(abs(z(1,1,3)-dZ*0.5D0).lt.dz*0.1) then
        do i=1,2
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,5-i)
        end do
     else if(abs(z(1,1,4)-dZ*0.5D0).lt.dz*0.1) then
        do i=1,3
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,7-i)
        end do
     else if(abs(z(1,1,5)-dZ*0.5D0).lt.dz*0.1) then
        do i=1,4
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,9-i)
        end do
     else if(abs(z(1,1,6)-dZ*0.5D0).lt.dz*0.1) then
        do i=1,5
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,11-i)
        end do
     else
        ! Uncomment the following only on 1 processor:
        !write(*,*) "This grid type is not supported!  Z-axis coordinates = ",ext,gridfunc_syms,z(1,1,:)
        !stop
     end if
  else 
     write(*,*) "Sorry, Symmetry!=1 not currently supported!"
     stop
  end if
end subroutine scalarwave_do_symmetry_gz_one_gridfunction

subroutine scalarwave_do_symmetry_gz_one_gridfunction_stagger(ext,X,Y,Z,Symmetry,gridfunc,gridfunc_syms, & 
                              stagger_x,stagger_y,stagger_z)
  implicit none

  integer, dimension(3)                    :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3))   :: X,Y,Z,gridfunc
  real*8                                   :: dX,dY,dZ
  integer				   :: stagger_x,stagger_y,stagger_z
  real*8                                   :: x_offset,y_offset,z_offset
  integer                                  :: Symmetry
  integer :: i

  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)

  if(Symmetry.eq.1) then
     z_offset = dZ*0.5d0*stagger_z
     if(abs(z(1,1,2)+z_offset).lt.dz*0.1) then
        gridfunc(:,:,1) = gridfunc_syms(3)*gridfunc(:,:,3)
     else if(abs(z(1,1,3)+z_offset).lt.dz*0.1) then
        do i=1,2
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,6-i)
        end do
     else if(abs(z(1,1,4)+z_offset).lt.dz*0.1) then
        do i=1,3
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,8-i)
        end do
     else if(abs(z(1,1,5)+z_offset).lt.dz*0.1) then
        do i=1,4
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,10-i)
        end do
     else if(abs(z(1,1,6)+z_offset).lt.dz*0.1) then
        do i=1,5
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,12-i)
        end do
        !=============================================
     else if(abs(z(1,1,2)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        gridfunc(:,:,1) = gridfunc_syms(3)*gridfunc(:,:,2)
     else if(abs(z(1,1,3)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,2
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,5-i)
        end do
     else if(abs(z(1,1,4)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,3
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,7-i)
        end do
     else if(abs(z(1,1,5)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,4
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,9-i)
        end do
     else if(abs(z(1,1,6)+z_offset-dZ*0.5D0).lt.dz*0.1) then
        do i=1,5
           gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,11-i)
        end do
     else
        ! Uncomment the following only on 1 processor:
        !write(*,*) "This grid type is not supported!  Z-axis coordinates = ",ext,gridfunc_syms,z(1,1,:)
        !stop
     end if
  else 
     write(*,*) "Sorry, Symmetry!=1 not currently supported!"
     stop
  end if
end subroutine scalarwave_do_symmetry_gz_one_gridfunction_stagger

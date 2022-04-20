!------------------------------------------------------------------------------------------------------
! Routines for filling symmetry ghostzones, including master function: do_symmetry_gz_one_gridfunction
!------------------------------------------------------------------------------------------------------
subroutine fill_bssn_symmetry_gz_gupij(ext,X,Y,Z,Symmetry,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz)
  implicit none
  integer, dimension(3)                  :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  integer                                :: Symmetry
  integer, parameter                     :: AXISYM = 4
  !
  if(Symmetry==1) then
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gupxx,gridfunc_syms)
     gridfunc_syms(1)=-1;gridfunc_syms(2)=-1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gupxy,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)=-1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gupxz,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gupyy,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)=-1;gridfunc_syms(3)=-1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gupyz,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gupzz,gridfunc_syms)
  end if
end subroutine fill_bssn_symmetry_gz_gupij
subroutine do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gridfunc,gridfunc_syms)
  implicit none
  integer, dimension(3)                    :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3))   :: X,Y,Z,gridfunc
  real*8                                   :: dX,dY,dZ
  integer                                  :: Symmetry
  integer :: i,j
  dX = X(2,1,1) - X(1,1,1)
  dY = Y(1,2,1) - Y(1,1,1)
  dZ = Z(1,1,2) - Z(1,1,1)
  if (Symmetry .ne. 0) then
     do j=2,ext(3)
        !Vertex centered (z=0,dz,2*dz,etc)
        if(abs(z(1,1,j)).lt.dz*0.1) then
           do i=1,j-1
              gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,2*j-i)
           end do
           !Cell centered (z=0.5*dz,1.5*dz,2.5*dz,etc)
        else if(abs(z(1,1,j)-dZ*0.5D0).lt.dz*0.1) then
           do i=1,j-1
              gridfunc(:,:,i) = gridfunc_syms(3)*gridfunc(:,:,2*j-i-1)
           end do
        end if
     end do
  end if
  !=============================================
end subroutine do_symmetry_gz_one_gridfunction
subroutine fill_bssn_symmetry_gz_bssn_vars(ext,X,Y,Z,Symmetry,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,&
  &Ayy,Ayz,Azz)
  implicit none
  integer, dimension(3)                  :: ext,gridfunc_syms
  real*8,dimension(ext(1),ext(2),ext(3)) :: X,Y,Z,phi,chi,trK,Gammax,Gammay,Gammaz,gxx,gxy,gxz,gyy,gyz,gzz,Axx,Axy,Axz,Ayy,Ayz,Azz
  integer                                :: Symmetry
  integer, parameter                     :: AXISYM = 4
  !
  if(Symmetry==1) then
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,phi,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,chi,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,trK,gridfunc_syms)
     gridfunc_syms(1)=-1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Gammax,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)=-1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Gammay,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)=-1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Gammaz,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gxx,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Axx,gridfunc_syms)
     gridfunc_syms(1)=-1;gridfunc_syms(2)=-1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gxy,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Axy,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)=-1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gxz,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Axz,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gyy,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Ayy,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)=-1;gridfunc_syms(3)=-1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gyz,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Ayz,gridfunc_syms)
     gridfunc_syms(1)= 1;gridfunc_syms(2)= 1;gridfunc_syms(3)= 1
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,gzz,gridfunc_syms)
     call do_symmetry_gz_one_gridfunction(ext,X,Y,Z,Symmetry,Azz,gridfunc_syms)
  end if
  !  write(*,*) finished with bssn symm. ghostzones
end subroutine fill_bssn_symmetry_gz_bssn_vars

!-----------------------------------------------------------------------------
! Fill in Us symmetry ghostzones
!-----------------------------------------------------------------------------
subroutine apply_U_syms(ext,cctk_nghostzones,Symmetry,X,Y,Z,U_syms,U)
  implicit none
  !INPUT PARAMETERS: 
  integer                                  :: Symmetry,cell_centering_enabled
  integer,dimension(3)                     :: ext,cctk_nghostzones,U_syms
  real*8,dimension(ext(1),ext(2),ext(3))   :: X,Y,Z,U
  !OTHER PARAMETERS: 
  integer                                  :: i,j,k,iend
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, AXISYM = 4)
  iend = cctk_nghostzones(1)
  if(Symmetry==AXISYM) iend=iend-1
  if((Symmetry==AXISYM .or. Symmetry==OCTANT) .and. X(1,1,1).lt.0.D0) then
     do i=1,iend
        U(i,:,:) = U_syms(1)*U(2*iend-i+1,:,:)
     end do
  end if
  if(Symmetry==OCTANT .and. Y(1,1,1).lt.0.D0) then
     do j=1,cctk_nghostzones(2)
        U(:,j,:) = U_syms(2)*U(:,2*cctk_nghostzones(2)-j+1,:)
     end do
  end if
  !Gotta be careful about Z(kmin) in axisymmetry, since we cant assume equatorial+axisymmetry in general:
  ! Still, theres a chance that the grid is split in such a way that (Z(1,1,1).lt.0.D0 .and. Z(1,1,cctk_nghostzones(3)+1).gt.0.D0) t
  if((Symmetry==EQUATORIAL .or. Symmetry==OCTANT .or. Symmetry==AXISYM) .and. &
       (Z(1,1,1).lt.0.D0 .and. Z(1,1,cctk_nghostzones(3)+1).gt.0.D0)) then
     do k=1,cctk_nghostzones(3)
        U(:,:,k) = U_syms(3)*U(:,:,2*cctk_nghostzones(3)-k+1)
     end do
  end if
end subroutine apply_U_syms

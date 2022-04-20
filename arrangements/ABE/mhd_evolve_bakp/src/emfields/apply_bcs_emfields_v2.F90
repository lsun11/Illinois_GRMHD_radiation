!---------------------------------------------------------------------!
! Update B^i outer boundaries					      !
!---------------------------------------------------------------------!
!
subroutine emfields_bc_newv2(ext,fake_ext, X,Y,Z, Bx,By,Bz, Bx_new, By_new, Bz_new, Symmetry, &
     have_bdry_min,have_bdry_max,bc)
  implicit none
  integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: Bx,By,Bz, Bx_new, By_new, Bz_new
  integer                                 :: Symmetry
  !
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  integer                            :: bc, FREEZE, EXTRAP, PERIODIC, COPY, QUAD
  integer                            :: PLANAR
  real*8, parameter		     :: TWO=2.d0, THREE = 3.d0
  real*8                             :: temp
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5)

  imin = ext(1)-fake_ext(1)+1
  jmin = ext(2)-fake_ext(2)+1
  kmin = ext(3)-fake_ext(3)+1

  imax = fake_ext(1)
  jmax = fake_ext(2)
  kmax = fake_ext(3)

  if(bc==FREEZE) then
     if(Symmetry.eq.AXISYM) then
        Bx_new(imax,2,:) = Bx(imax,2,:)
        By_new(imax,2,:) = By(imax,2,:)
        Bz_new(imax,2,:) = Bz(imax,2,:)
        Bx_new(:,2,kmax) = Bx(:,2,kmax)
        By_new(:,2,kmax) = By(:,2,kmax)
        Bz_new(:,2,kmax) = Bz(:,2,kmax)

        if (have_bdry_max(1)==1) then
           Bx_new(imax-1,2,:) = Bx(imax-1,2,:)
           By_new(imax-1,2,:) = By(imax-1,2,:)
           Bz_new(imax-1,2,:) = Bz(imax-1,2,:)
           Bx_new(imax-2,2,:) = Bx(imax-2,2,:)
           By_new(imax-2,2,:) = By(imax-2,2,:)
           Bz_new(imax-2,2,:) = Bz(imax-2,2,:)
        end if
        if (have_bdry_max(3)==1) then
           Bx_new(:,2,kmax-1) = Bx(:,2,kmax-1)
           By_new(:,2,kmax-1) = By(:,2,kmax-1)
           Bz_new(:,2,kmax-1) = Bz(:,2,kmax-1)
           Bx_new(:,2,kmax-2) = Bx(:,2,kmax-2)
           By_new(:,2,kmax-2) = By(:,2,kmax-2)
           Bz_new(:,2,kmax-2) = Bz(:,2,kmax-2)
        end if
        if (Z(1,1,kmin) .lt. 0.d0) then
           Bx_new(:,:,kmin) = Bx(:,:,kmin)
           By_new(:,:,kmin) = By(:,:,kmin)
           Bz_new(:,:,kmin) = Bz(:,:,kmin)
        end if
     else
        Bx_new(imax,:,:) = Bx(imax,:,:)
        By_new(imax,:,:) = By(imax,:,:)
        Bz_new(imax,:,:) = Bz(imax,:,:)
        Bx_new(:,:,kmax) = Bx(:,:,kmax)
        By_new(:,:,kmax) = By(:,:,kmax)
        Bz_new(:,:,kmax) = Bz(:,:,kmax)
        Bx_new(:,jmax,:) = Bx(:,jmax,:)
        By_new(:,jmax,:) = By(:,jmax,:)
        Bz_new(:,jmax,:) = Bz(:,jmax,:)

        if (have_bdry_max(1)==1) then
           Bx_new(imax-1,:,:) = Bx(imax-1,:,:)
           By_new(imax-1,:,:) = By(imax-1,:,:)
           Bz_new(imax-1,:,:) = Bz(imax-1,:,:)
           Bx_new(imax-2,:,:) = Bx(imax-2,:,:)
           By_new(imax-2,:,:) = By(imax-2,:,:)
           Bz_new(imax-2,:,:) = Bz(imax-2,:,:)
        end if
        if (have_bdry_max(3)==1) then
           Bx_new(:,:,kmax-1) = Bx(:,:,kmax-1)
           By_new(:,:,kmax-1) = By(:,:,kmax-1)
           Bz_new(:,:,kmax-1) = Bz(:,:,kmax-1)
           Bx_new(:,:,kmax-2) = Bx(:,:,kmax-2)
           By_new(:,:,kmax-2) = By(:,:,kmax-2)
           Bz_new(:,:,kmax-2) = Bz(:,:,kmax-2)
        end if
        if (have_bdry_max(2)==1) then
           Bx_new(:,jmax-1,:) = Bx(:,jmax-1,:)
           By_new(:,jmax-1,:) = By(:,jmax-1,:)
           Bz_new(:,jmax-1,:) = Bz(:,jmax-1,:)
           Bx_new(:,jmax-2,:) = Bx(:,jmax-2,:)
           By_new(:,jmax-2,:) = By(:,jmax-2,:)
           Bz_new(:,jmax-2,:) = Bz(:,jmax-2,:)
        end if
        if (have_bdry_min(1)==1) then
           Bx_new(imin,:,:) = Bx(imin,:,:)
           By_new(imin,:,:) = By(imin,:,:)
           Bz_new(imin,:,:) = Bz(imin,:,:)
        end if
        if (have_bdry_min(2)==1) then
           Bx_new(:,jmin,:) = Bx(:,jmin,:)
           By_new(:,jmin,:) = By(:,jmin,:)
           Bz_new(:,jmin,:) = Bz(:,jmin,:)
        end if
        if (have_bdry_min(3)==1) then
           Bx_new(:,:,kmin) = Bx(:,:,kmin)
           By_new(:,:,kmin) = By(:,:,kmin)
           Bz_new(:,:,kmin) = Bz(:,:,kmin)
        end if
     end if
  else if(bc==EXTRAP) then
     if (have_bdry_max(1)==1) then
        Bx_new(imax,:,:) = TWO*Bx_new(imax-1,:,:) & 
             - Bx_new(imax-2,:,:)
        By_new(imax,:,:) = TWO*By_new(imax-1,:,:) & 
             - By_new(imax-2,:,:)
        Bz_new(imax,:,:) = TWO*Bz_new(imax-1,:,:) & 
             - Bz_new(imax-2,:,:)
     end if

     if (Symmetry .ne. AXISYM .and. have_bdry_max(2)==1) then
        Bx_new(:,jmax,:) = TWO*Bx_new(:,jmax-1,:) &
             - Bx_new(:,jmax-2,:)
        By_new(:,jmax,:) = TWO*By_new(:,jmax-1,:) &
             - By_new(:,jmax-2,:)
        Bz_new(:,jmax,:) = TWO*Bz_new(:,jmax-1,:) &
             - Bz_new(:,jmax-2,:)
     end if

     if (have_bdry_max(3)==1) then
        Bx_new(:,:,kmax) = TWO*Bx_new(:,:,kmax-1) &
             - Bx_new(:,:,kmax-2)
        By_new(:,:,kmax) = TWO*By_new(:,:,kmax-1) &
             - By_new(:,:,kmax-2)
        Bz_new(:,:,kmax) = TWO*Bz_new(:,:,kmax-1) &
             - Bz_new(:,:,kmax-2)
     end if


     if (have_bdry_min(1)==1) then
        Bx_new(imin,:,:) = TWO*Bx_new(imin+1,:,:) &
             - Bx_new(imin+2,:,:)
        By_new(imin,:,:) = TWO*By_new(imin+1,:,:) &
             - By_new(imin+2,:,:)
        Bz_new(imin,:,:) = TWO*Bz_new(imin+1,:,:) &
             - Bz_new(imin+2,:,:)
     end if
     if (have_bdry_min(2)==1) then
        Bx_new(:,jmin,:) = TWO*Bx_new(:,jmin+1,:) &
             - Bx_new(:,jmin+2,:)
        By_new(:,jmin,:) = TWO*By_new(:,jmin+1,:) &
             - By_new(:,jmin+2,:)
        Bz_new(:,jmin,:) = TWO*Bz_new(:,jmin+1,:) &
             - Bz_new(:,jmin+2,:)
     end if
     if (have_bdry_min(3)==1) then
        Bx_new(:,:,kmin) = TWO*Bx_new(:,:,kmin+1) &
             - Bx_new(:,:,kmin+2)
        By_new(:,:,kmin) = TWO*By_new(:,:,kmin+1) &
             - By_new(:,:,kmin+2)
        Bz_new(:,:,kmin) = TWO*Bz_new(:,:,kmin+1) &
             - Bz_new(:,:,kmin+2)
     end if

     if (Symmetry==AXISYM .and. Z(1,1,kmin) .lt. 0.d0) then
        Bx_new(:,:,kmin) = TWO*Bx_new(:,:,kmin+1) &
             - Bx_new(:,:,kmin+2)
        By_new(:,:,kmin) = TWO*By_new(:,:,kmin+1) &
             - By_new(:,:,kmin+2)
        Bz_new(:,:,kmin) = TWO*Bz_new(:,:,kmin+1) &
             - Bz_new(:,:,kmin+2)
     end if

  else if(bc==PLANAR) then
! FIXME: DO BC==PLANAR
     write(*,*) "BC==PLANAR not yet supported!"
     stop

     do k=kmin+1,kmax-1
        do j=jmin,jmax
           do i=imin,imax
              if(i==imax) then
                 if(j==jmax) then
                    Bx_new(i,j,k) = Bx_new(i-1,j-1,k)
                    By_new(i,j,k) = By_new(i-1,j-1,k)
                    Bz_new(i,j,k) = Bz_new(i-1,j-1,k)
                 else if(j==jmax-1) then
                    Bx_new(i,j,k) = Bx_new(i-1,j,k)
                    By_new(i,j,k) = By_new(i-1,j,k)
                    Bz_new(i,j,k) = Bz_new(i-1,j,k)
                 else
                    Bx_new(i,j,k) = Bx_new(i-1,j+1,k)
                    By_new(i,j,k) = By_new(i-1,j+1,k)
                    Bz_new(i,j,k) = Bz_new(i-1,j+1,k)
                 end if
              else if(i==imin) then
                 if(j==jmin) then
                    Bx_new(i,j,k) = Bx_new(i+1,j+1,k)
                    By_new(i,j,k) = By_new(i+1,j+1,k)
                    Bz_new(i,j,k) = Bz_new(i+1,j+1,k)
                 else if(j==jmin+1) then
                    Bx_new(i,j,k) = Bx_new(i+1,j,k)
                    By_new(i,j,k) = By_new(i+1,j,k)
                    Bz_new(i,j,k) = Bz_new(i+1,j,k)
                 else
                    Bx_new(i,j,k) = Bx_new(i+1,j-1,k)
                    By_new(i,j,k) = By_new(i+1,j-1,k)
                    Bz_new(i,j,k) = Bz_new(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    Bx_new(i,j,k) = Bx_new(i,j+1,k)
                    By_new(i,j,k) = By_new(i,j+1,k)
                    Bz_new(i,j,k) = Bz_new(i,j+1,k)
                 else
                    Bx_new(i,j,k) = Bx_new(i-1,j+1,k)
                    By_new(i,j,k) = By_new(i-1,j+1,k)
                    Bz_new(i,j,k) = Bz_new(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    Bx_new(i,j,k) = Bx_new(i,j-1,k)
                    By_new(i,j,k) = By_new(i,j-1,k)
                    Bz_new(i,j,k) = Bz_new(i,j-1,k)
                 else
                    Bx_new(i,j,k) = Bx_new(i+1,j-1,k)
                    By_new(i,j,k) = By_new(i+1,j-1,k)
                    Bz_new(i,j,k) = Bz_new(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do

  else if (bc==COPY) then     
     if (have_bdry_max(1)==1) then
        Bx_new(imax,:,:) = Bx_new(imax-1,:,:)
        By_new(imax,:,:) = By_new(imax-1,:,:)
        Bz_new(imax,:,:) = Bz_new(imax-1,:,:)
     end if

     if (have_bdry_min(1)==1) then
        Bx_new(imin,:,:) = Bx_new(imin+1,:,:)
        By_new(imin,:,:) = By_new(imin+1,:,:)
        Bz_new(imin,:,:) = Bz_new(imin+1,:,:)
     end if

     if(have_bdry_max(2)==1) then
        Bx_new(:,jmax,:) = Bx_new(:,jmax-1,:)
        By_new(:,jmax,:) = By_new(:,jmax-1,:)
        Bz_new(:,jmax,:) = Bz_new(:,jmax-1,:)
     end if

     if (have_bdry_min(2)==1) then
        Bx_new(:,jmin,:) = Bx_new(:,jmin+1,:)
        By_new(:,jmin,:) = By_new(:,jmin+1,:)
        Bz_new(:,jmin,:) = Bz_new(:,jmin+1,:)
     end if

     if(have_bdry_max(3)==1) then
        Bx_new(:,:,kmax) = Bx_new(:,:,kmax-1)
        By_new(:,:,kmax) = By_new(:,:,kmax-1)
        Bz_new(:,:,kmax) = Bz_new(:,:,kmax-1)
     end if
     if(have_bdry_min(3)==1) then
        Bx_new(:,:,kmin) = Bx_new(:,:,kmin+1)
        By_new(:,:,kmin) = By_new(:,:,kmin+1)
        Bz_new(:,:,kmin) = Bz_new(:,:,kmin+1)
     end if
  else if(bc==QUAD) then

     if (have_bdry_max(1)==1) then
        Bx_new(imax,:,:) = THREE*Bx_new(imax-1,:,:) & 
             - THREE*Bx_new(imax-2,:,:) + Bx_new(imax-3,:,:)
        By_new(imax,:,:) = THREE*By_new(imax-1,:,:) & 
             - THREE*By_new(imax-2,:,:) + By_new(imax-3,:,:)
        Bz_new(imax,:,:) = THREE*Bz_new(imax-1,:,:) & 
             - THREE*Bz_new(imax-2,:,:) + Bz_new(imax-3,:,:)
     end if
     if (have_bdry_max(2)==1) then
        Bx_new(:,jmax,:) = THREE*Bx_new(:,jmax-1,:) &
             - THREE*Bx_new(:,jmax-2,:) + Bx_new(:,jmax-3,:)
        By_new(:,jmax,:) = THREE*By_new(:,jmax-1,:) &
             - THREE*By_new(:,jmax-2,:) + By_new(:,jmax-3,:)
        Bz_new(:,jmax,:) = THREE*Bz_new(:,jmax-1,:) &
             - THREE*Bz_new(:,jmax-2,:) + Bz_new(:,jmax-3,:)
     end if
     if (have_bdry_max(3)==1) then
        Bx_new(:,:,kmax) = THREE*Bx_new(:,:,kmax-1) &
             - THREE*Bx_new(:,:,kmax-2) + Bx_new(:,:,kmax-3)
        By_new(:,:,kmax) = THREE*By_new(:,:,kmax-1) &
             - THREE*By_new(:,:,kmax-2) + By_new(:,:,kmax-3)
        Bz_new(:,:,kmax) = THREE*Bz_new(:,:,kmax-1) &
             - THREE*Bz_new(:,:,kmax-2) + Bz_new(:,:,kmax-3)
     end if
     if (have_bdry_min(1)==1) then
        Bx_new(imin,:,:) = THREE*Bx_new(imin+1,:,:) &
             - THREE*Bx_new(imin+2,:,:) + Bx_new(imin+3,:,:)
        By_new(imin,:,:) = THREE*By_new(imin+1,:,:) &
             - THREE*By_new(imin+2,:,:) + By_new(imin+3,:,:)
        Bz_new(imin,:,:) = THREE*Bz_new(imin+1,:,:) &
             - THREE*Bz_new(imin+2,:,:) + Bz_new(imin+3,:,:)
     end if

     if (have_bdry_min(2)==1) then
        Bx_new(:,jmin,:) = THREE*Bx_new(:,jmin+1,:) &
             - THREE*Bx_new(:,jmin+2,:) + Bx_new(:,jmin+3,:)
        By_new(:,jmin,:) = THREE*By_new(:,jmin+1,:) &
             - THREE*By_new(:,jmin+2,:) + By_new(:,jmin+3,:)
        Bz_new(:,jmin,:) = THREE*Bz_new(:,jmin+1,:) &
             - THREE*Bz_new(:,jmin+2,:) + Bz_new(:,jmin+3,:)
     end if

     if (have_bdry_min(3)==1) then
        Bx_new(:,:,kmin) = THREE*Bx_new(:,:,kmin+1) &
             - THREE*Bx_new(:,:,kmin+2) + Bx_new(:,:,kmin+3)
        By_new(:,:,kmin) = THREE*By_new(:,:,kmin+1) &
             - THREE*By_new(:,:,kmin+2) + By_new(:,:,kmin+3)
        Bz_new(:,:,kmin) = THREE*Bz_new(:,:,kmin+1) &
             - THREE*Bz_new(:,:,kmin+2) + Bz_new(:,:,kmin+3)
     end if
  end if

end subroutine emfields_bc_newv2

subroutine emfields_bc_Ai(ext,fake_ext, X,Y,Z, Ax,Ay,Az, Ax_new, Ay_new, Az_new, Symmetry, &
     have_bdry_min,have_bdry_max,bc)
  implicit none
  integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: Ax,Ay,Az, Ax_new, Ay_new, Az_new
  integer                                 :: Symmetry
  !
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  integer                            :: bc, FREEZE, EXTRAP, PERIODIC, COPY, QUAD
  integer                            :: PLANAR
  real*8, parameter		     :: TWO=2.d0, THREE = 3.d0
  real*8                             :: temp
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5)

  imin = ext(1)-fake_ext(1)+1
  jmin = ext(2)-fake_ext(2)+1
  kmin = ext(3)-fake_ext(3)+1

  imax = fake_ext(1)
  jmax = fake_ext(2)
  kmax = fake_ext(3)

  if(bc==FREEZE) then
    write(*,*) 'Frozen BC not supported for A evolution.'
    stop
  else if(bc==COPY) then
     !COPY is the default boundary condition for A, as EM_BC is currently set in param.ccl
     ! Ax is defined on the semi-staggered grid (i,j+1/2,k+1/2)
     ! Ay is defined on the semi-staggered grid (i+1/2,j,k+1/2)
     ! Az is defined on the semi-staggered grid (i+1/2,j+1/2,k)
     ! We do linear extrapolation here, since we apply copy boundary conditions on Bi's,
     !  and we need curl Ai = Bi.  
     ! In other words, if we just used plain copy bc's on A, curl Ai=0 at the boundary, 
     !  which would be inconsistent with the value of Bi at the boundary.
     if (have_bdry_max(1)==1) then
        Ax_new(imax,:,:) = TWO*Ax_new(imax-1,:,:) & 
             - Ax_new(imax-2,:,:)
        Ay_new(imax,:,:) = TWO*Ay_new(imax-1,:,:) & 
             - Ay_new(imax-2,:,:)
        Az_new(imax,:,:) = TWO*Az_new(imax-1,:,:) & 
             - Az_new(imax-2,:,:)
     end if

     if (Symmetry .ne. AXISYM .and. have_bdry_max(2)==1) then
        Ax_new(:,jmax,:) = TWO*Ax_new(:,jmax-1,:) &
             - Ax_new(:,jmax-2,:)
        Ay_new(:,jmax,:) = TWO*Ay_new(:,jmax-1,:) &
             - Ay_new(:,jmax-2,:)
        Az_new(:,jmax,:) = TWO*Az_new(:,jmax-1,:) &
             - Az_new(:,jmax-2,:)
     end if

     if (have_bdry_max(3)==1) then
        Ax_new(:,:,kmax) = TWO*Ax_new(:,:,kmax-1) &
             - Ax_new(:,:,kmax-2)
        Ay_new(:,:,kmax) = TWO*Ay_new(:,:,kmax-1) &
             - Ay_new(:,:,kmax-2)
        Az_new(:,:,kmax) = TWO*Az_new(:,:,kmax-1) &
             - Az_new(:,:,kmax-2)
     end if


     if (have_bdry_min(1)==1) then
        Ax_new(imin,:,:) = TWO*Ax_new(imin+1,:,:) &
             - Ax_new(imin+2,:,:)
        Ay_new(imin,:,:) = TWO*Ay_new(imin+1,:,:) &
             - Ay_new(imin+2,:,:)
        Az_new(imin,:,:) = TWO*Az_new(imin+1,:,:) &
             - Az_new(imin+2,:,:)
     end if
     if (have_bdry_min(2)==1) then
        Ax_new(:,jmin,:) = TWO*Ax_new(:,jmin+1,:) &
             - Ax_new(:,jmin+2,:)
        Ay_new(:,jmin,:) = TWO*Ay_new(:,jmin+1,:) &
             - Ay_new(:,jmin+2,:)
        Az_new(:,jmin,:) = TWO*Az_new(:,jmin+1,:) &
             - Az_new(:,jmin+2,:)
     end if
     if (have_bdry_min(3)==1) then
        Ax_new(:,:,kmin) = TWO*Ax_new(:,:,kmin+1) &
             - Ax_new(:,:,kmin+2)
        Ay_new(:,:,kmin) = TWO*Ay_new(:,:,kmin+1) &
             - Ay_new(:,:,kmin+2)
        Az_new(:,:,kmin) = TWO*Az_new(:,:,kmin+1) &
             - Az_new(:,:,kmin+2)
     end if

     if (Symmetry==AXISYM .and. Z(1,1,kmin) .lt. 0.d0) then
        Ax_new(:,:,kmin) = TWO*Ax_new(:,:,kmin+1) &
             - Ax_new(:,:,kmin+2)
        Ay_new(:,:,kmin) = TWO*Ay_new(:,:,kmin+1) &
             - Ay_new(:,:,kmin+2)
        Az_new(:,:,kmin) = TWO*Az_new(:,:,kmin+1) &
             - Az_new(:,:,kmin+2)
     end if

  else if(bc==PLANAR) then
! FIXME: DO BC==PLANAR
     write(*,*) "BC==PLANAR not yet supported!"
     stop

     do k=kmin+1,kmax-1
        do j=jmin,jmax
           do i=imin,imax
              if(i==imax) then
                 if(j==jmax) then
                    Ax_new(i,j,k) = Ax_new(i-1,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j-1,k)
                    Az_new(i,j,k) = Az_new(i-1,j-1,k)
                 else if(j==jmax-1) then
                    Ax_new(i,j,k) = Ax_new(i-1,j,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j,k)
                    Az_new(i,j,k) = Az_new(i-1,j,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i-1,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j+1,k)
                    Az_new(i,j,k) = Az_new(i-1,j+1,k)
                 end if
              else if(i==imin) then
                 if(j==jmin) then
                    Ax_new(i,j,k) = Ax_new(i+1,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j+1,k)
                    Az_new(i,j,k) = Az_new(i+1,j+1,k)
                 else if(j==jmin+1) then
                    Ax_new(i,j,k) = Ax_new(i+1,j,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j,k)
                    Az_new(i,j,k) = Az_new(i+1,j,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i+1,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j-1,k)
                    Az_new(i,j,k) = Az_new(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    Ax_new(i,j,k) = Ax_new(i,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i,j+1,k)
                    Az_new(i,j,k) = Az_new(i,j+1,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i-1,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j+1,k)
                    Az_new(i,j,k) = Az_new(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    Ax_new(i,j,k) = Ax_new(i,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i,j-1,k)
                    Az_new(i,j,k) = Az_new(i,j-1,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i+1,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j-1,k)
                    Az_new(i,j,k) = Az_new(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do

  else if (bc==QUAD) then     
       write(*,*) 'QUAD BC not supported for A evolution'
       stop
  else if(bc==EXTRAP) then

     if (have_bdry_max(1)==1) then
        Ax_new(imax,:,:) = THREE*Ax_new(imax-1,:,:) & 
             - THREE*Ax_new(imax-2,:,:) + Ax_new(imax-3,:,:)
        Ay_new(imax,:,:) = THREE*Ay_new(imax-1,:,:) & 
             - THREE*Ay_new(imax-2,:,:) + Ay_new(imax-3,:,:)
        Az_new(imax,:,:) = THREE*Az_new(imax-1,:,:) & 
             - THREE*Az_new(imax-2,:,:) + Az_new(imax-3,:,:)
     end if
     if (have_bdry_max(2)==1) then
        Ax_new(:,jmax,:) = THREE*Ax_new(:,jmax-1,:) &
             - THREE*Ax_new(:,jmax-2,:) + Ax_new(:,jmax-3,:)
        Ay_new(:,jmax,:) = THREE*Ay_new(:,jmax-1,:) &
             - THREE*Ay_new(:,jmax-2,:) + Ay_new(:,jmax-3,:)
        Az_new(:,jmax,:) = THREE*Az_new(:,jmax-1,:) &
             - THREE*Az_new(:,jmax-2,:) + Az_new(:,jmax-3,:)
     end if
     if (have_bdry_max(3)==1) then
        Ax_new(:,:,kmax) = THREE*Ax_new(:,:,kmax-1) &
             - THREE*Ax_new(:,:,kmax-2) + Ax_new(:,:,kmax-3)
        Ay_new(:,:,kmax) = THREE*Ay_new(:,:,kmax-1) &
             - THREE*Ay_new(:,:,kmax-2) + Ay_new(:,:,kmax-3)
        Az_new(:,:,kmax) = THREE*Az_new(:,:,kmax-1) &
             - THREE*Az_new(:,:,kmax-2) + Az_new(:,:,kmax-3)
     end if
     if (have_bdry_min(1)==1) then
        Ax_new(imin,:,:) = THREE*Ax_new(imin+1,:,:) &
             - THREE*Ax_new(imin+2,:,:) + Ax_new(imin+3,:,:)
        Ay_new(imin,:,:) = THREE*Ay_new(imin+1,:,:) &
             - THREE*Ay_new(imin+2,:,:) + Ay_new(imin+3,:,:)
        Az_new(imin,:,:) = THREE*Az_new(imin+1,:,:) &
             - THREE*Az_new(imin+2,:,:) + Az_new(imin+3,:,:)
     end if

     if (have_bdry_min(2)==1) then
        Ax_new(:,jmin,:) = THREE*Ax_new(:,jmin+1,:) &
             - THREE*Ax_new(:,jmin+2,:) + Ax_new(:,jmin+3,:)
        Ay_new(:,jmin,:) = THREE*Ay_new(:,jmin+1,:) &
             - THREE*Ay_new(:,jmin+2,:) + Ay_new(:,jmin+3,:)
        Az_new(:,jmin,:) = THREE*Az_new(:,jmin+1,:) &
             - THREE*Az_new(:,jmin+2,:) + Az_new(:,jmin+3,:)
     end if

     if (have_bdry_min(3)==1) then
        Ax_new(:,:,kmin) = THREE*Ax_new(:,:,kmin+1) &
             - THREE*Ax_new(:,:,kmin+2) + Ax_new(:,:,kmin+3)
        Ay_new(:,:,kmin) = THREE*Ay_new(:,:,kmin+1) &
             - THREE*Ay_new(:,:,kmin+2) + Ay_new(:,:,kmin+3)
        Az_new(:,:,kmin) = THREE*Az_new(:,:,kmin+1) &
             - THREE*Az_new(:,:,kmin+2) + Az_new(:,:,kmin+3)
     end if
  end if

end subroutine emfields_bc_Ai

subroutine emfields_bc_Ai_psi6phi(ext,fake_ext, X,Y,Z, psi6phi,Ax,Ay,Az, psi6phi_new,Ax_new,Ay_new,Az_new,Symmetry, &
     have_bdry_min,have_bdry_max,bc)
  implicit none
  integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: Ax,Ay,Az, Ax_new, Ay_new, Az_new
  real*8, dimension(ext(1),ext(2),ext(3)) :: psi6phi,psi6phi_new
  integer                                 :: Symmetry
  !
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  integer                            :: bc, FREEZE, EXTRAP, PERIODIC, COPY, QUAD
  integer                            :: PLANAR
  real*8, parameter		     :: TWO=2.d0, THREE = 3.d0
  real*8                             :: temp
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5)

  imin = ext(1)-fake_ext(1)+1
  jmin = ext(2)-fake_ext(2)+1
  kmin = ext(3)-fake_ext(3)+1

  imax = fake_ext(1)
  jmax = fake_ext(2)
  kmax = fake_ext(3)

  if(bc==FREEZE) then
    write(*,*) 'Frozen BC not supported for A evolution.'
    stop
  else if(bc==COPY) then
     !COPY is the default boundary condition for A, as EM_BC is currently set in param.ccl
     ! Ax is defined on the semi-staggered grid (i,j+1/2,k+1/2)
     ! Ay is defined on the semi-staggered grid (i+1/2,j,k+1/2)
     ! Az is defined on the semi-staggered grid (i+1/2,j+1/2,k)
     ! We do linear extrapolation here, since we apply copy boundary conditions on Bi's,
     !  and we need curl Ai = Bi.  
     ! In other words, if we just used plain copy bc's on A, curl Ai=0 at the boundary, 
     !  which would be inconsistent with the value of Bi at the boundary.
     if (have_bdry_max(1)==1) then
        Ax_new(imax,:,:) = TWO*Ax_new(imax-1,:,:) & 
             - Ax_new(imax-2,:,:)
        Ay_new(imax,:,:) = TWO*Ay_new(imax-1,:,:) & 
             - Ay_new(imax-2,:,:)
        Az_new(imax,:,:) = TWO*Az_new(imax-1,:,:) & 
             - Az_new(imax-2,:,:)
	psi6phi_new(imax,:,:) = TWO*psi6phi_new(imax-1,:,:) &
             - psi6phi_new(imax-2,:,:)
     end if

     if (Symmetry .ne. AXISYM .and. have_bdry_max(2)==1) then
        Ax_new(:,jmax,:) = TWO*Ax_new(:,jmax-1,:) &
             - Ax_new(:,jmax-2,:)
        Ay_new(:,jmax,:) = TWO*Ay_new(:,jmax-1,:) &
             - Ay_new(:,jmax-2,:)
        Az_new(:,jmax,:) = TWO*Az_new(:,jmax-1,:) &
             - Az_new(:,jmax-2,:)
        psi6phi_new(:,jmax,:) = TWO*psi6phi_new(:,jmax-1,:) &
             - psi6phi_new(:,jmax-2,:)
     end if

     if (have_bdry_max(3)==1) then
        Ax_new(:,:,kmax) = TWO*Ax_new(:,:,kmax-1) &
             - Ax_new(:,:,kmax-2)
        Ay_new(:,:,kmax) = TWO*Ay_new(:,:,kmax-1) &
             - Ay_new(:,:,kmax-2)
        Az_new(:,:,kmax) = TWO*Az_new(:,:,kmax-1) &
             - Az_new(:,:,kmax-2)
	psi6phi_new(:,:,kmax) = TWO*psi6phi_new(:,:,kmax-1) &
             - psi6phi_new(:,:,kmax-2)
     end if


     if (have_bdry_min(1)==1) then
        Ax_new(imin,:,:) = TWO*Ax_new(imin+1,:,:) &
             - Ax_new(imin+2,:,:)
        Ay_new(imin,:,:) = TWO*Ay_new(imin+1,:,:) &
             - Ay_new(imin+2,:,:)
        Az_new(imin,:,:) = TWO*Az_new(imin+1,:,:) &
             - Az_new(imin+2,:,:)
	psi6phi_new(imin,:,:) = TWO*psi6phi_new(imin+1,:,:) &
             - psi6phi_new(imin+2,:,:)
     end if
     if (have_bdry_min(2)==1) then
        Ax_new(:,jmin,:) = TWO*Ax_new(:,jmin+1,:) &
             - Ax_new(:,jmin+2,:)
        Ay_new(:,jmin,:) = TWO*Ay_new(:,jmin+1,:) &
             - Ay_new(:,jmin+2,:)
        Az_new(:,jmin,:) = TWO*Az_new(:,jmin+1,:) &
             - Az_new(:,jmin+2,:)
        psi6phi_new(:,jmin,:) = TWO*psi6phi_new(:,jmin+1,:) &
             - psi6phi_new(:,jmin+2,:)
     end if
     if (have_bdry_min(3)==1) then
        Ax_new(:,:,kmin) = TWO*Ax_new(:,:,kmin+1) &
             - Ax_new(:,:,kmin+2)
        Ay_new(:,:,kmin) = TWO*Ay_new(:,:,kmin+1) &
             - Ay_new(:,:,kmin+2)
        Az_new(:,:,kmin) = TWO*Az_new(:,:,kmin+1) &
             - Az_new(:,:,kmin+2)
	psi6phi_new(:,:,kmin) = TWO*psi6phi_new(:,:,kmin+1) &
             - psi6phi_new(:,:,kmin+2)
     end if

     if (Symmetry==AXISYM .and. Z(1,1,kmin) .lt. 0.d0) then
        Ax_new(:,:,kmin) = TWO*Ax_new(:,:,kmin+1) &
             - Ax_new(:,:,kmin+2)
        Ay_new(:,:,kmin) = TWO*Ay_new(:,:,kmin+1) &
             - Ay_new(:,:,kmin+2)
        Az_new(:,:,kmin) = TWO*Az_new(:,:,kmin+1) &
             - Az_new(:,:,kmin+2)
        psi6phi_new(:,:,kmin) = TWO*psi6phi_new(:,:,kmin+1) &
             - psi6phi_new(:,:,kmin+2)
     end if

  else if(bc==PLANAR) then
! FIXME: DO BC==PLANAR
     write(*,*) "BC==PLANAR not yet supported!"
     stop

     do k=kmin+1,kmax-1
        do j=jmin,jmax
           do i=imin,imax
              if(i==imax) then
                 if(j==jmax) then
                    Ax_new(i,j,k) = Ax_new(i-1,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j-1,k)
                    Az_new(i,j,k) = Az_new(i-1,j-1,k)
                 else if(j==jmax-1) then
                    Ax_new(i,j,k) = Ax_new(i-1,j,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j,k)
                    Az_new(i,j,k) = Az_new(i-1,j,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i-1,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j+1,k)
                    Az_new(i,j,k) = Az_new(i-1,j+1,k)
                 end if
              else if(i==imin) then
                 if(j==jmin) then
                    Ax_new(i,j,k) = Ax_new(i+1,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j+1,k)
                    Az_new(i,j,k) = Az_new(i+1,j+1,k)
                 else if(j==jmin+1) then
                    Ax_new(i,j,k) = Ax_new(i+1,j,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j,k)
                    Az_new(i,j,k) = Az_new(i+1,j,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i+1,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j-1,k)
                    Az_new(i,j,k) = Az_new(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    Ax_new(i,j,k) = Ax_new(i,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i,j+1,k)
                    Az_new(i,j,k) = Az_new(i,j+1,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i-1,j+1,k)
                    Ay_new(i,j,k) = Ay_new(i-1,j+1,k)
                    Az_new(i,j,k) = Az_new(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    Ax_new(i,j,k) = Ax_new(i,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i,j-1,k)
                    Az_new(i,j,k) = Az_new(i,j-1,k)
                 else
                    Ax_new(i,j,k) = Ax_new(i+1,j-1,k)
                    Ay_new(i,j,k) = Ay_new(i+1,j-1,k)
                    Az_new(i,j,k) = Az_new(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do

  else if (bc==QUAD) then     
       write(*,*) 'QUAD BC not supported for A evolution'
       stop
  else if(bc==EXTRAP) then

     if (have_bdry_max(1)==1) then
        Ax_new(imax,:,:) = THREE*Ax_new(imax-1,:,:) & 
             - THREE*Ax_new(imax-2,:,:) + Ax_new(imax-3,:,:)
        Ay_new(imax,:,:) = THREE*Ay_new(imax-1,:,:) & 
             - THREE*Ay_new(imax-2,:,:) + Ay_new(imax-3,:,:)
        Az_new(imax,:,:) = THREE*Az_new(imax-1,:,:) & 
             - THREE*Az_new(imax-2,:,:) + Az_new(imax-3,:,:)
	psi6phi_new(imax,:,:) = THREE*psi6phi_new(imax-1,:,:) &
             - THREE*psi6phi_new(imax-2,:,:) + psi6phi_new(imax-3,:,:)
     end if
     if (have_bdry_max(2)==1) then
        Ax_new(:,jmax,:) = THREE*Ax_new(:,jmax-1,:) &
             - THREE*Ax_new(:,jmax-2,:) + Ax_new(:,jmax-3,:)
        Ay_new(:,jmax,:) = THREE*Ay_new(:,jmax-1,:) &
             - THREE*Ay_new(:,jmax-2,:) + Ay_new(:,jmax-3,:)
        Az_new(:,jmax,:) = THREE*Az_new(:,jmax-1,:) &
             - THREE*Az_new(:,jmax-2,:) + Az_new(:,jmax-3,:)
        psi6phi_new(:,jmax,:) = THREE*psi6phi_new(:,jmax-1,:) &
             - THREE*psi6phi_new(:,jmax-2,:) + psi6phi_new(:,jmax-3,:)
     end if
     if (have_bdry_max(3)==1) then
        Ax_new(:,:,kmax) = THREE*Ax_new(:,:,kmax-1) &
             - THREE*Ax_new(:,:,kmax-2) + Ax_new(:,:,kmax-3)
        Ay_new(:,:,kmax) = THREE*Ay_new(:,:,kmax-1) &
             - THREE*Ay_new(:,:,kmax-2) + Ay_new(:,:,kmax-3)
        Az_new(:,:,kmax) = THREE*Az_new(:,:,kmax-1) &
             - THREE*Az_new(:,:,kmax-2) + Az_new(:,:,kmax-3)
	psi6phi_new(:,:,kmax) = THREE*psi6phi_new(:,:,kmax-1) &
             - THREE*psi6phi_new(:,:,kmax-2) + psi6phi_new(:,:,kmax-3)
     end if
     if (have_bdry_min(1)==1) then
        Ax_new(imin,:,:) = THREE*Ax_new(imin+1,:,:) &
             - THREE*Ax_new(imin+2,:,:) + Ax_new(imin+3,:,:)
        Ay_new(imin,:,:) = THREE*Ay_new(imin+1,:,:) &
             - THREE*Ay_new(imin+2,:,:) + Ay_new(imin+3,:,:)
        Az_new(imin,:,:) = THREE*Az_new(imin+1,:,:) &
             - THREE*Az_new(imin+2,:,:) + Az_new(imin+3,:,:)
        psi6phi_new(imin,:,:) = THREE*psi6phi_new(imin+1,:,:) &
             - THREE*psi6phi_new(imin+2,:,:) + psi6phi_new(imin+3,:,:)
     end if

     if (have_bdry_min(2)==1) then
        Ax_new(:,jmin,:) = THREE*Ax_new(:,jmin+1,:) &
             - THREE*Ax_new(:,jmin+2,:) + Ax_new(:,jmin+3,:)
        Ay_new(:,jmin,:) = THREE*Ay_new(:,jmin+1,:) &
             - THREE*Ay_new(:,jmin+2,:) + Ay_new(:,jmin+3,:)
        Az_new(:,jmin,:) = THREE*Az_new(:,jmin+1,:) &
             - THREE*Az_new(:,jmin+2,:) + Az_new(:,jmin+3,:)
        psi6phi_new(:,jmin,:) = THREE*psi6phi_new(:,jmin+1,:) &
             - THREE*psi6phi_new(:,jmin+2,:) + psi6phi_new(:,jmin+3,:)
     end if

     if (have_bdry_min(3)==1) then
        Ax_new(:,:,kmin) = THREE*Ax_new(:,:,kmin+1) &
             - THREE*Ax_new(:,:,kmin+2) + Ax_new(:,:,kmin+3)
        Ay_new(:,:,kmin) = THREE*Ay_new(:,:,kmin+1) &
             - THREE*Ay_new(:,:,kmin+2) + Ay_new(:,:,kmin+3)
        Az_new(:,:,kmin) = THREE*Az_new(:,:,kmin+1) &
             - THREE*Az_new(:,:,kmin+2) + Az_new(:,:,kmin+3)
        psi6phi_new(:,:,kmin) = THREE*psi6phi_new(:,:,kmin+1) &
             - THREE*psi6phi_new(:,:,kmin+2) + psi6phi_new(:,:,kmin+3)
     end if
  end if

end subroutine emfields_bc_Ai_psi6phi

subroutine update_boundary_field_line_tracer(ext,fake_ext, Symmetry, &
     mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi, &
     mhd_psi_line_p, mhd_u_psi_p, mhd_chi_line_p, mhd_u_chi_p, &
     have_bdry_min,have_bdry_max,bc)
  implicit none
  integer, dimension(3)                   :: ext,fake_ext,have_bdry_min,have_bdry_max
  real*8, dimension(ext(1),ext(2),ext(3)) :: mhd_psi_line, mhd_u_psi, mhd_chi_line, mhd_u_chi
  real*8, dimension(ext(1),ext(2),ext(3)) :: mhd_psi_line_p, mhd_u_psi_p, mhd_chi_line_p, mhd_u_chi_p
  integer                                 :: Symmetry
  !
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  integer                            :: bc, FREEZE, EXTRAP, PERIODIC, COPY, QUAD
  integer                            :: PLANAR
  real*8, parameter                  :: TWO=2.d0, THREE = 3.d0
  real*8                             :: temp
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(COPY = 1, FREEZE = 2, EXTRAP = 3, QUAD = 4, PLANAR = 5)

  imin = ext(1)-fake_ext(1)+1
  jmin = ext(2)-fake_ext(2)+1
  kmin = ext(3)-fake_ext(3)+1

  imax = fake_ext(1)
  jmax = fake_ext(2)
  kmax = fake_ext(3)

  if (Symmetry==AXISYM) then
     write(*,*) 'Field line tracer does not support axisymmetry!'
     stop
  end if

  if(bc==FREEZE) then

     if (have_bdry_max(1)==1) then
        i = imax
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              mhd_psi_line(i,j,k) = mhd_psi_line_p(i,j,k)
              mhd_chi_line(i,j,k) = mhd_chi_line_p(i,j,k)
              mhd_u_psi(i,j,k) = mhd_u_psi_p(i,j,k)
              mhd_u_chi(i,j,k) = mhd_u_chi_p(i,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(2)==1) then
        j = jmax
        !$omp parallel do
        do k=1,ext(3)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = mhd_psi_line_p(i,j,k)
              mhd_chi_line(i,j,k) = mhd_chi_line_p(i,j,k)
              mhd_u_psi(i,j,k) = mhd_u_psi_p(i,j,k)
              mhd_u_chi(i,j,k) = mhd_u_chi_p(i,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(2)==1) then
        k = kmax
        !$omp parallel do
        do j=1,ext(2)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = mhd_psi_line_p(i,j,k)
              mhd_chi_line(i,j,k) = mhd_chi_line_p(i,j,k)
              mhd_u_psi(i,j,k) = mhd_u_psi_p(i,j,k)
              mhd_u_chi(i,j,k) = mhd_u_chi_p(i,j,k)
           end do
        end do
        !$omp end parallel do
     end if
     if (have_bdry_min(1)==1) then
        i = imin
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              mhd_psi_line(i,j,k) = mhd_psi_line_p(i,j,k)
              mhd_chi_line(i,j,k) = mhd_chi_line_p(i,j,k)
              mhd_u_psi(i,j,k) = mhd_u_psi_p(i,j,k)
              mhd_u_chi(i,j,k) = mhd_u_chi_p(i,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(2)==1) then
        j = jmin
        !$omp parallel do
        do k=1,ext(3)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = mhd_psi_line_p(i,j,k)
              mhd_chi_line(i,j,k) = mhd_chi_line_p(i,j,k)
              mhd_u_psi(i,j,k) = mhd_u_psi_p(i,j,k)
              mhd_u_chi(i,j,k) = mhd_u_chi_p(i,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(2)==1) then
        k = kmin
        !$omp parallel do
        do j=1,ext(2)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = mhd_psi_line_p(i,j,k)
              mhd_chi_line(i,j,k) = mhd_chi_line_p(i,j,k)
              mhd_u_psi(i,j,k) = mhd_u_psi_p(i,j,k)
              mhd_u_chi(i,j,k) = mhd_u_chi_p(i,j,k)
           end do
        end do
        !$omp end parallel do
     end if

  else if(bc==COPY) then
     if (have_bdry_max(1)==1) then
        i = imax
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              mhd_psi_line(i,j,k) = 2.d0*mhd_psi_line(i-1,j,k) - mhd_psi_line(i-2,j,k)
              mhd_chi_line(i,j,k) = 2.d0*mhd_chi_line(i-1,j,k) - mhd_chi_line(i-2,j,k)
              mhd_u_psi(i,j,k) = 2.d0*mhd_u_psi(i-1,j,k) - mhd_u_psi(i-2,j,k)
              mhd_u_chi(i,j,k) = 2.d0*mhd_u_chi(i-1,j,k) - mhd_u_chi(i-2,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(1)==1) then
        i = imin
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              mhd_psi_line(i,j,k) = 2.d0*mhd_psi_line(i+1,j,k) - mhd_psi_line(i+2,j,k)
              mhd_chi_line(i,j,k) = 2.d0*mhd_chi_line(i+1,j,k) - mhd_chi_line(i+2,j,k)
              mhd_u_psi(i,j,k) = 2.d0*mhd_u_psi(i+1,j,k) - mhd_u_psi(i+2,j,k)
              mhd_u_chi(i,j,k) = 2.d0*mhd_u_chi(i+1,j,k) - mhd_u_chi(i+2,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(2)==1) then
        j = jmax
        !$omp parallel do
        do k=1,ext(3)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 2.d0*mhd_psi_line(i,j-1,k) - mhd_psi_line(i,j-2,k)
              mhd_chi_line(i,j,k) = 2.d0*mhd_chi_line(i,j-1,k) - mhd_chi_line(i,j-2,k)
              mhd_u_psi(i,j,k) = 2.d0*mhd_u_psi(i,j-1,k) - mhd_u_psi(i,j-2,k)
              mhd_u_chi(i,j,k) = 2.d0*mhd_u_chi(i,j-1,k) - mhd_u_chi(i,j-2,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(2)==1) then
        j = jmin
        !$omp parallel do
        do k=1,ext(3)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 2.d0*mhd_psi_line(i,j+1,k) - mhd_psi_line(i,j+2,k)
              mhd_chi_line(i,j,k) = 2.d0*mhd_chi_line(i,j+1,k) - mhd_chi_line(i,j+2,k)
              mhd_u_psi(i,j,k) = 2.d0*mhd_u_psi(i,j+1,k) - mhd_u_psi(i,j+2,k)
              mhd_u_chi(i,j,k) = 2.d0*mhd_u_chi(i,j+1,k) - mhd_u_chi(i,j+2,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(3)==1) then
        k = kmax
        !$omp parallel do
        do j=1,ext(2)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 2.d0*mhd_psi_line(i,j,k-1) - mhd_psi_line(i,j,k-2)
              mhd_chi_line(i,j,k) = 2.d0*mhd_chi_line(i,j,k-1) - mhd_chi_line(i,j,k-2)
              mhd_u_psi(i,j,k) = 2.d0*mhd_u_psi(i,j,k-1) - mhd_u_psi(i,j,k-2)
              mhd_u_chi(i,j,k) = 2.d0*mhd_u_chi(i,j,k-1) - mhd_u_chi(i,j,k-2)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(3)==1) then
        k = kmin
        !$omp parallel do
        do j=1,ext(2)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 2.d0*mhd_psi_line(i,j,k+1) - mhd_psi_line(i,j,k+2)
              mhd_chi_line(i,j,k) = 2.d0*mhd_chi_line(i,j,k+1) - mhd_chi_line(i,j,k+2)
              mhd_u_psi(i,j,k) = 2.d0*mhd_u_psi(i,j,k+1) - mhd_u_psi(i,j,k+2)
              mhd_u_chi(i,j,k) = 2.d0*mhd_u_chi(i,j,k+1) - mhd_u_chi(i,j,k+2)
           end do
        end do
        !$omp end parallel do
     end if

  else if(bc==PLANAR) then
! FIXME: DO BC==PLANAR
     write(*,*) "BC==PLANAR not yet supported!"
     stop

  !!else if (bc==COPY) then

  !!   if (have_bdry_max(1)==1) then
  !!      i = imax
  !!      !$omp parallel do
  !!      do k=1,ext(3)
  !!         do j=1,ext(2)
  !!            mhd_psi_line(i,j,k) = mhd_psi_line(i-1,j,k)
  !!            mhd_chi_line(i,j,k) = mhd_chi_line(i-1,j,k)
  !!            mhd_u_psi(i,j,k) = mhd_u_psi(i-1,j,k)
  !!            mhd_u_chi(i,j,k) = mhd_u_chi(i-1,j,k)
  !!         end do
  !!      end do
  !!      !$omp end parallel do
  !!   end if

  !!   if (have_bdry_min(1)==1) then
  !!      i = imin
  !!      !$omp parallel do
  !!      do k=1,ext(3)
  !!         do j=1,ext(2)
  !!            mhd_psi_line(i,j,k) = mhd_psi_line(i+1,j,k)
  !!            mhd_chi_line(i,j,k) = mhd_chi_line(i+1,j,k)
  !!            mhd_u_psi(i,j,k) = mhd_u_psi(i+1,j,k)
  !!            mhd_u_chi(i,j,k) = mhd_u_chi(i+1,j,k)
  !!         end do
  !!      end do
  !!      !$omp end parallel do
  !!   end if

  !!   if (have_bdry_max(2)==1) then
  !!      j = jmax
  !!      !$omp parallel do
  !!      do k=1,ext(3)
  !!         do i=1,ext(1)
  !!            mhd_psi_line(i,j,k) = mhd_psi_line(i,j-1,k)
  !!            mhd_chi_line(i,j,k) = mhd_chi_line(i,j-1,k)
  !!            mhd_u_psi(i,j,k) = mhd_u_psi(i,j-1,k)
  !!            mhd_u_chi(i,j,k) = mhd_u_chi(i,j-1,k)
  !!         end do
  !!      end do
  !!      !$omp end parallel do
  !!   end if

  !!   if (have_bdry_min(2)==1) then
  !!      j = jmin
  !!      !$omp parallel do
  !!      do k=1,ext(3)
  !!         do i=1,ext(1)
  !!            mhd_psi_line(i,j,k) = mhd_psi_line(i,j+1,k)
  !!            mhd_chi_line(i,j,k) = mhd_chi_line(i,j+1,k)
  !!            mhd_u_psi(i,j,k) = mhd_u_psi(i,j+1,k)
  !!            mhd_u_chi(i,j,k) = mhd_u_chi(i,j+1,k)
  !!         end do
  !!      end do
  !!      !$omp end parallel do
  !!   end if

  !!   if (have_bdry_max(3)==1) then
  !!      k = kmax
  !!      !$omp parallel do
  !!      do j=1,ext(2)
  !!         do i=1,ext(1)
  !!            mhd_psi_line(i,j,k) = mhd_psi_line(i,j,k-1)
  !!            mhd_chi_line(i,j,k) = mhd_chi_line(i,j,k-1)
  !!            mhd_u_psi(i,j,k) = mhd_u_psi(i,j,k-1)
  !!            mhd_u_chi(i,j,k) = mhd_u_chi(i,j,k-1)
  !!         end do
  !!      end do
  !!      !$omp end parallel do
  !!   end if

  !!   if (have_bdry_min(3)==1) then
  !!      k = kmin
  !!      !$omp parallel do
  !!      do j=1,ext(2)
  !!         do i=1,ext(1)
  !!            mhd_psi_line(i,j,k) = mhd_psi_line(i,j,k+1)
  !!            mhd_chi_line(i,j,k) = mhd_chi_line(i,j,k+1)
  !!            mhd_u_psi(i,j,k) = mhd_u_psi(i,j,k+1)
  !!            mhd_u_chi(i,j,k) = mhd_u_chi(i,j,k+1)
  !!         end do
  !!      end do
  !!      !$omp end parallel do
  !!   end if

  else if(bc==QUAD .or. bc==EXTRAP) then
     if (have_bdry_max(1)==1) then
        i = imax
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              mhd_psi_line(i,j,k) = 3.d0*(mhd_psi_line(i-1,j,k) - mhd_psi_line(i-2,j,k)) + mhd_psi_line(i-3,j,k)
              mhd_chi_line(i,j,k) = 3.d0*(mhd_chi_line(i-1,j,k) - mhd_chi_line(i-2,j,k)) + mhd_chi_line(i-3,j,k)
              mhd_u_psi(i,j,k)    = 3.d0*(   mhd_u_psi(i-1,j,k) -    mhd_u_psi(i-2,j,k)) +    mhd_u_psi(i-3,j,k)
              mhd_u_chi(i,j,k)    = 3.d0*(   mhd_u_chi(i-1,j,k) -    mhd_u_chi(i-2,j,k)) +    mhd_u_chi(i-3,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(1)==1) then
        i = imin
        !$omp parallel do
        do k=1,ext(3)
           do j=1,ext(2)
              mhd_psi_line(i,j,k) = 3.d0*(mhd_psi_line(i+1,j,k) - mhd_psi_line(i+2,j,k)) + mhd_psi_line(i+3,j,k)
              mhd_chi_line(i,j,k) = 3.d0*(mhd_chi_line(i+1,j,k) - mhd_chi_line(i+2,j,k)) + mhd_chi_line(i+3,j,k)
              mhd_u_psi(i,j,k)    = 3.d0*(   mhd_u_psi(i+1,j,k) -    mhd_u_psi(i+2,j,k)) +    mhd_u_psi(i+3,j,k)
              mhd_u_chi(i,j,k)    = 3.d0*(   mhd_u_chi(i+1,j,k) -    mhd_u_chi(i+2,j,k)) +    mhd_u_chi(i+3,j,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(2)==1) then
        j = jmax
        !$omp parallel do
        do k=1,ext(3)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 3.d0*(mhd_psi_line(i,j-1,k) - mhd_psi_line(i,j-2,k)) + mhd_psi_line(i,j-3,k)
              mhd_chi_line(i,j,k) = 3.d0*(mhd_chi_line(i,j-1,k) - mhd_chi_line(i,j-2,k)) + mhd_chi_line(i,j-3,k)
              mhd_u_psi(i,j,k)    = 3.d0*(   mhd_u_psi(i,j-1,k) -    mhd_u_psi(i,j-2,k)) +    mhd_u_psi(i,j-3,k)
              mhd_u_chi(i,j,k)    = 3.d0*(   mhd_u_chi(i,j-1,k) -    mhd_u_chi(i,j-2,k)) +    mhd_u_chi(i,j-3,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(2)==1) then
        j = jmin
        !$omp parallel do
        do k=1,ext(3)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 3.d0*(mhd_psi_line(i,j+1,k) - mhd_psi_line(i,j+2,k)) + mhd_psi_line(i,j+3,k)
              mhd_chi_line(i,j,k) = 3.d0*(mhd_chi_line(i,j+1,k) - mhd_chi_line(i,j+2,k)) + mhd_chi_line(i,j+3,k)
              mhd_u_psi(i,j,k)    = 3.d0*(   mhd_u_psi(i,j+1,k) -    mhd_u_psi(i,j+2,k)) +    mhd_u_psi(i,j+3,k)
              mhd_u_chi(i,j,k)    = 3.d0*(   mhd_u_chi(i,j+1,k) -    mhd_u_chi(i,j+2,k)) +    mhd_u_chi(i,j+3,k)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_max(3)==1) then
        k = kmax
        !$omp parallel do
        do j=1,ext(2)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 3.d0*(mhd_psi_line(i,j,k-1) - mhd_psi_line(i,j,k-2)) + mhd_psi_line(i,j,k-3)
              mhd_chi_line(i,j,k) = 3.d0*(mhd_chi_line(i,j,k-1) - mhd_chi_line(i,j,k-2)) + mhd_chi_line(i,j,k-3)
              mhd_u_psi(i,j,k)    = 3.d0*(   mhd_u_psi(i,j,k-1) -    mhd_u_psi(i,j,k-2)) +    mhd_u_psi(i,j,k-3)
              mhd_u_chi(i,j,k)    = 3.d0*(   mhd_u_chi(i,j,k-1) -    mhd_u_chi(i,j,k-2)) +    mhd_u_chi(i,j,k-3)
           end do
        end do
        !$omp end parallel do
     end if

     if (have_bdry_min(3)==1) then
        k = kmin
        !$omp parallel do
        do j=1,ext(2)
           do i=1,ext(1)
              mhd_psi_line(i,j,k) = 3.d0*(mhd_psi_line(i,j,k+1) - mhd_psi_line(i,j,k+2)) + mhd_psi_line(i,j,k+3)
              mhd_chi_line(i,j,k) = 3.d0*(mhd_chi_line(i,j,k+1) - mhd_chi_line(i,j,k+2)) + mhd_chi_line(i,j,k+3)
              mhd_u_psi(i,j,k)    = 3.d0*(   mhd_u_psi(i,j,k+1) -    mhd_u_psi(i,j,k+2)) +    mhd_u_psi(i,j,k+3)
              mhd_u_chi(i,j,k)    = 3.d0*(   mhd_u_chi(i,j,k+1) -    mhd_u_chi(i,j,k+2)) +    mhd_u_chi(i,j,k+3)
           end do
        end do
        !$omp end parallel do
     end if

  end if

end subroutine update_boundary_field_line_tracer

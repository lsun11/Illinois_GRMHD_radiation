!-----------------------------------------------------------------------------
!
! Hydro evolution routines
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
! Compute coordinate velocities (see eq. (2.12) in Shibata, PRD 60, 104052)
!
!-----------------------------------------------------------------------------
subroutine find_v(ex,rho_star,shiftx,shifty,shiftz,                        &
        phi,lapse,w,h,st_x,st_y,st_z,                                 &
             gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,vx,vy,vz)
  implicit none
  integer, dimension(3)                             :: ex
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: rho_star
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: shiftx, shifty, shiftz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: phi, lapse, w, h
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: st_x, st_y, st_z
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(in)  :: gupyy,gupyz,gupzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(out) :: vx, vy, vz
!
! other variables
!
  integer              :: imin,jmin,kmin,imax,jmax,kmax
  integer              :: i, j, k
  real*8, parameter    :: ONE  = 1.D0, FOUR = 4.D0, Tiny = 1.D-10
  real*8               :: temp
  imin = lbound(phi,1)
  jmin = lbound(phi,2)
  kmin = lbound(phi,3)
  imax = ubound(phi,1)
  jmax = ubound(phi,2)
  kmax = ubound(phi,3)
  do k = kmin, kmax
     do j = jmin, jmax
        do i = imin, imax
           temp = h(i,j,k) + 1.d0
           temp = w(i,j,k) + 1.d0
!           if ((rho_star(i,j,k).gt.Tiny).and.(w(i,j,k).gt.Tiny)) then
           vx(i,j,k) = - shiftx(i,j,k) + exp(-FOUR*phi(i,j,k)) *         &
                (lapse(i,j,k) + ONE)/(h(i,j,k)*w(i,j,k)) *                 &
                ( gupxx(i,j,k) * st_x(i,j,k) + gupxy(i,j,k) * st_y(i,j,k)  &
                + gupxz(i,j,k) * st_z(i,j,k))
           vy(i,j,k) = - shifty(i,j,k) + exp(-FOUR*phi(i,j,k)) *        &
                (lapse(i,j,k) + ONE)/(h(i,j,k)*w(i,j,k)) *                 &
                ( gupxy(i,j,k) * st_x(i,j,k) + gupyy(i,j,k) * st_y(i,j,k)  &
                + gupyz(i,j,k) * st_z(i,j,k))
           vz(i,j,k) = - shiftz(i,j,k) + exp(-FOUR*phi(i,j,k)) *        &
                (lapse(i,j,k) + ONE)/(h(i,j,k)*w(i,j,k)) *                 &
                ( gupxz(i,j,k) * st_x(i,j,k) + gupyz(i,j,k) * st_y(i,j,k)  &
                + gupzz(i,j,k) * st_z(i,j,k))
!           else
!              vx(i,j,k) = 0.D0
!              vy(i,j,k) = 0.D0
!              vz(i,j,k) = 0.D0
!           endif
        end do
     end do
  end do
end subroutine find_v
!-----------------------------------------------------------------------------
!
! Update outer boundaries
!
!-----------------------------------------------------------------------------
subroutine hydro_matter_bc(ex, Z, &
     rho, tau, st_x, st_y, st_z, &
     rho_new, tau_new,st_x_new,st_y_new,st_z_new,Symmetry,bc)
 !    proc_imax,proc_jmax,proc_kmax,glob_imax,glob_jmax,glob_kmax,bc)
  implicit none
  integer, dimension(3)                :: ex
  real*8, dimension(ex(1),ex(2),ex(3)) :: Z
  real*8, dimension(ex(1),ex(2),ex(3)) :: rho,tau,st_x
  real*8, dimension(ex(1),ex(2),ex(3)) :: st_y,st_z
  real*8, dimension(ex(1),ex(2),ex(3)) :: rho_new,tau_new,st_x_new
  real*8, dimension(ex(1),ex(2),ex(3)) :: st_y_new,st_z_new
  integer                              :: Symmetry
!
! Other variables:
! 
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  integer                            :: proc_imax,proc_jmax,proc_kmax
  integer                            :: glob_imax,glob_jmax,glob_kmax
  real*8, parameter                  :: ZERO = 0.D0
  real*8, parameter                  :: TWO = 2.D0, THREE = 3.D0
  real*8                             :: nx,ny,nz
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT
  integer                            :: PI_SYMM, AXISYM
  integer                            :: bc, FREEZE, EXTRAP, PERIODIC, OUTF, COPY, QUAD
  integer                            :: PLANAR
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(OUTF = 1, FREEZE = 2, COPY = 3, EXTRAP = 4, QUAD = 5, PLANAR = 6)
!
! Input translation
!
  imin = lbound(tau_new,1)
  jmin = lbound(tau_new,2)
  kmin = lbound(tau_new,3)
  imax = ubound(tau_new,1)
  jmax = ubound(tau_new,2)
  kmax = ubound(tau_new,3)
  if(bc==FREEZE) then
     if(Symmetry.eq.AXISYM) then
        rho_new(imax,2,:) = rho(imax,2,:)
        tau_new(imax,2,:) = tau(imax,2,:)
        st_x_new(imax,2,:) = st_x(imax,2,:)
        st_y_new(imax,2,:) = st_y(imax,2,:)
        st_z_new(imax,2,:) = st_z(imax,2,:)
        rho_new(:,2,kmax) = rho(:,2,kmax)
        tau_new(:,2,kmax) = tau(:,2,kmax)
        st_x_new(:,2,kmax) = st_x(:,2,kmax)
        st_y_new(:,2,kmax) = st_y(:,2,kmax)
        st_z_new(:,2,kmax) = st_z(:,2,kmax)
        if(proc_imax==glob_imax) then
           rho_new(imax-1,2,:) = rho(imax-1,2,:)
           tau_new(imax-1,2,:) = tau(imax-1,2,:)
           st_x_new(imax-1,2,:) = st_x(imax-1,2,:)
           st_y_new(imax-1,2,:) = st_y(imax-1,2,:)
           st_z_new(imax-1,2,:) = st_z(imax-1,2,:)
           rho_new(imax-2,2,:) = rho(imax-2,2,:)
           tau_new(imax-2,2,:) = tau(imax-2,2,:)
           st_x_new(imax-2,2,:) = st_x(imax-2,2,:)
           st_y_new(imax-2,2,:) = st_y(imax-2,2,:)
           st_z_new(imax-2,2,:) = st_z(imax-2,2,:)
           rho_new(imax-3,2,:) = rho(imax-3,2,:)
           tau_new(imax-3,2,:) = tau(imax-3,2,:)
           st_x_new(imax-3,2,:) = st_x(imax-3,2,:)
           st_y_new(imax-3,2,:) = st_y(imax-3,2,:)
           st_z_new(imax-3,2,:) = st_z(imax-3,2,:)
        end if
        if(proc_kmax==glob_kmax) then
           rho_new(:,2,kmax-1) = rho(:,2,kmax-1)
           tau_new(:,2,kmax-1) = tau(:,2,kmax-1)
           st_x_new(:,2,kmax-1) = st_x(:,2,kmax-1)
           st_y_new(:,2,kmax-1) = st_y(:,2,kmax-1)
           st_z_new(:,2,kmax-1) = st_z(:,2,kmax-1)
           rho_new(:,2,kmax-2) = rho(:,2,kmax-2)
           tau_new(:,2,kmax-2) = tau(:,2,kmax-2)
           st_x_new(:,2,kmax-2) = st_x(:,2,kmax-2)
           st_y_new(:,2,kmax-2) = st_y(:,2,kmax-2)
           st_z_new(:,2,kmax-2) = st_z(:,2,kmax-2)
           rho_new(:,2,kmax-3) = rho(:,2,kmax-3)
           tau_new(:,2,kmax-3) = tau(:,2,kmax-3)
           st_x_new(:,2,kmax-3) = st_x(:,2,kmax-3)
           st_y_new(:,2,kmax-3) = st_y(:,2,kmax-3)
           st_z_new(:,2,kmax-3) = st_z(:,2,kmax-3)
        end if
     else
        rho_new(imax,:,:) = rho(imax,:,:)
        tau_new(imax,:,:) = tau(imax,:,:)
        st_x_new(imax,:,:) = st_x(imax,:,:)
        st_y_new(imax,:,:) = st_y(imax,:,:)
        st_z_new(imax,:,:) = st_z(imax,:,:)
        rho_new(:,:,kmax) = rho(:,:,kmax)
        tau_new(:,:,kmax) = tau(:,:,kmax)
        st_x_new(:,:,kmax) = st_x(:,:,kmax)
        st_y_new(:,:,kmax) = st_y(:,:,kmax)
        st_z_new(:,:,kmax) = st_z(:,:,kmax)
        rho_new(:,jmax,:) = rho(:,jmax,:)
        tau_new(:,jmax,:) = tau(:,jmax,:)
        st_x_new(:,jmax,:) = st_x(:,jmax,:)
        st_y_new(:,jmax,:) = st_y(:,jmax,:)
        st_z_new(:,jmax,:) = st_z(:,jmax,:)
        if(proc_imax==glob_imax) then
           rho_new(imax-1,:,:) = rho(imax-1,:,:)
           tau_new(imax-1,:,:) = tau(imax-1,:,:)
           st_x_new(imax-1,:,:) = st_x(imax-1,:,:)
           st_y_new(imax-1,:,:) = st_y(imax-1,:,:)
           st_z_new(imax-1,:,:) = st_z(imax-1,:,:)
           rho_new(imax-2,:,:) = rho(imax-2,:,:)
           tau_new(imax-2,:,:) = tau(imax-2,:,:)
           st_x_new(imax-2,:,:) = st_x(imax-2,:,:)
           st_y_new(imax-2,:,:) = st_y(imax-2,:,:)
           st_z_new(imax-2,:,:) = st_z(imax-2,:,:)
        end if
        if(proc_kmax==glob_kmax) then
           rho_new(:,:,kmax-1) = rho(:,:,kmax-1)
           tau_new(:,:,kmax-1) = tau(:,:,kmax-1)
           st_x_new(:,:,kmax-1) = st_x(:,:,kmax-1)
           st_y_new(:,:,kmax-1) = st_y(:,:,kmax-1)
           st_z_new(:,:,kmax-1) = st_z(:,:,kmax-1)
           rho_new(:,:,kmax-2) = rho(:,:,kmax-2)
           tau_new(:,:,kmax-2) = tau(:,:,kmax-2)
           st_x_new(:,:,kmax-2) = st_x(:,:,kmax-2)
           st_y_new(:,:,kmax-2) = st_y(:,:,kmax-2)
           st_z_new(:,:,kmax-2) = st_z(:,:,kmax-2)
        end if
        if(proc_jmax==glob_jmax) then
           rho_new(:,jmax-1,:) = rho(:,jmax-1,:)
           tau_new(:,jmax-1,:) = tau(:,jmax-1,:)
           st_x_new(:,jmax-1,:) = st_x(:,jmax-1,:)
           st_y_new(:,jmax-1,:) = st_y(:,jmax-1,:)
           st_z_new(:,jmax-1,:) = st_z(:,jmax-1,:)
           rho_new(:,jmax-2,:) = rho(:,jmax-2,:)
           tau_new(:,jmax-2,:) = tau(:,jmax-2,:)
           st_x_new(:,jmax-2,:) = st_x(:,jmax-2,:)
           st_y_new(:,jmax-2,:) = st_y(:,jmax-2,:)
           st_z_new(:,jmax-2,:) = st_z(:,jmax-2,:)
        end if
        if (Symmetry .ne. OCTANT) then
           rho_new(imin,:,:) = rho(imin,:,:)
           tau_new(imin,:,:) = tau(imin,:,:)
           st_x_new(imin,:,:) = st_x(imin,:,:)
           st_y_new(imin,:,:) = st_y(imin,:,:)
           st_z_new(imin,:,:) = st_z(imin,:,:)
           if (Symmetry .ne. PI_SYMM) then
              rho_new(:,jmin,:) = rho(:,jmin,:)
              tau_new(:,jmin,:) = tau(:,jmin,:)
              st_x_new(:,jmin,:) = st_x(:,jmin,:)
              st_y_new(:,jmin,:) = st_y(:,jmin,:)
              st_z_new(:,jmin,:) = st_z(:,jmin,:)
              if (Symmetry == NO_SYMM .or. Z(1,1,kmin+1) .lt. 0.d0) then
                 rho_new(:,:,kmin) = rho(:,:,kmin)
                 tau_new(:,:,kmin) = tau(:,:,kmin)
                 st_x_new(:,:,kmin) = st_x(:,:,kmin)
                 st_y_new(:,:,kmin) = st_y(:,:,kmin)
                 st_z_new(:,:,kmin) = st_z(:,:,kmin)
              end if
           endif
        endif
     end if
  else if(bc==OUTF) then
     if(Symmetry.eq.AXISYM) then
        rho_new(imax,:,:) = rho_new(imax-1,:,:)
        tau_new(imax,:,:) = tau_new(imax-1,:,:)
        st_x_new(imax,:,:) = st_x_new(imax-1,:,:)
        st_y_new(imax,:,:) = st_y_new(imax-1,:,:)
        st_z_new(imax,:,:) = st_z_new(imax-1,:,:)
        do k=kmin,kmax
           do j=jmin,jmax
              if(st_x_new(imax,j,k).lt.ZERO) st_x_new(imax,j,k) = ZERO
           end do
        end do
        rho_new(:,:,kmax) = rho_new(:,:,kmax-1)
        tau_new(:,:,kmax) = tau_new(:,:,kmax-1)
        st_x_new(:,:,kmax) = st_x_new(:,:,kmax-1)
        st_y_new(:,:,kmax) = st_y_new(:,:,kmax-1)
        st_z_new(:,:,kmax) = st_z_new(:,:,kmax-1)
        do j=jmin,jmax
           do i=imin,imax
              if(st_z_new(i,j,kmax).lt.ZERO) st_z_new(i,j,kmax) = ZERO
           end do
        end do
        if (Z(1,1,kmin+1) .lt. 0.d0) then
           rho_new(:,:,kmin) = rho_new(:,:,kmin+1)
           tau_new(:,:,kmin) = tau_new(:,:,kmin+1)
           st_x_new(:,:,kmin) = st_x_new(:,:,kmin+1)
           st_y_new(:,:,kmin) = st_y_new(:,:,kmin+1)
           st_z_new(:,:,kmin) = st_z_new(:,:,kmin+1)
           do i=imin,imax
              do j=jmin,jmax
                 if(st_z_new(i,j,kmin).gt.ZERO) st_z_new(i,j,kmin) = ZERO
              end do
           end do
        end if
     else
        rho_new(imax,:,:) = rho_new(imax-1,:,:)
        tau_new(imax,:,:) = tau_new(imax-1,:,:)
        st_x_new(imax,:,:) = st_x_new(imax-1,:,:)
        st_y_new(imax,:,:) = st_y_new(imax-1,:,:)
        st_z_new(imax,:,:) = st_z_new(imax-1,:,:)
        do j=jmin,jmax
           do k=kmin,kmax
              if(st_x_new(imax,j,k).lt.ZERO) st_x_new(imax,j,k) = ZERO
           end do
        end do
        rho_new(:,jmax,:) = rho_new(:,jmax-1,:)
        tau_new(:,jmax,:) = tau_new(:,jmax-1,:)
        st_x_new(:,jmax,:) = st_x_new(:,jmax-1,:)
        st_y_new(:,jmax,:) = st_y_new(:,jmax-1,:)
        st_z_new(:,jmax,:) = st_z_new(:,jmax-1,:)
        do i=imin,imax
           do k=kmin,kmax
              if(st_y_new(i,jmax,k).lt.ZERO) st_y_new(i,jmax,k) = ZERO
           end do
        end do
        rho_new(:,:,kmax) = rho_new(:,:,kmax-1)
        tau_new(:,:,kmax) = tau_new(:,:,kmax-1)
        st_x_new(:,:,kmax) = st_x_new(:,:,kmax-1)
        st_y_new(:,:,kmax) = st_y_new(:,:,kmax-1)
        st_z_new(:,:,kmax) = st_z_new(:,:,kmax-1)
        do i=imin,imax
           do j=jmin,jmax
              if(st_z_new(i,j,kmax).lt.ZERO) st_z_new(i,j,kmax) = ZERO
           end do
        end do
        if (Symmetry .ne. OCTANT) then
           rho_new(imin,:,:) = rho_new(imin+1,:,:)
           tau_new(imin,:,:) = tau_new(imin+1,:,:)
           st_x_new(imin,:,:) = st_x_new(imin+1,:,:)
           st_y_new(imin,:,:) = st_y_new(imin+1,:,:)
           st_z_new(imin,:,:) = st_z_new(imin+1,:,:)
           do k=kmin,kmax
              do j=jmin,jmax
                 if(st_x_new(imin,j,k).gt.ZERO) st_x_new(imin,j,k) = ZERO
              end do
           end do
           if (Symmetry .ne. PI_SYMM) then
              rho_new(:,jmin,:) = rho_new(:,jmin+1,:)
              tau_new(:,jmin,:) = tau_new(:,jmin+1,:)
              st_x_new(:,jmin,:) = st_x_new(:,jmin+1,:)
              st_y_new(:,jmin,:) = st_y_new(:,jmin+1,:)
              st_z_new(:,jmin,:) = st_z_new(:,jmin+1,:)
              do k=kmin,kmax
                 do i=imin,imax
                    if(st_y_new(i,jmin,k).gt.ZERO) st_y_new(i,jmin,k) = ZERO
                 end do
              end do
              if (Symmetry == NO_SYMM .or. Z(1,1,kmin+1) .lt. 0.d0) then
                 rho_new(:,:,kmin) = rho_new(:,:,kmin+1)
                 tau_new(:,:,kmin) = tau_new(:,:,kmin+1)
                 st_x_new(:,:,kmin) = st_x_new(:,:,kmin+1)
                 st_y_new(:,:,kmin) = st_y_new(:,:,kmin+1)
                 st_z_new(:,:,kmin) = st_z_new(:,:,kmin+1)
                 do i=imin,imax
                    do j=jmin,jmax
                       if(st_z_new(i,j,kmin).gt.ZERO) st_z_new(i,j,kmin) = ZERO
                    end do
                 end do
              end if
           endif
        endif
     end if
  else if(bc==EXTRAP) then
     rho_new(imax,:,:) = TWO*rho_new(imax-1,:,:) &
          - rho_new(imax-2,:,:)
     tau_new(imax,:,:) = TWO*tau_new(imax-1,:,:) &
          - tau_new(imax-2,:,:)
     st_x_new(imax,:,:) = TWO*st_x_new(imax-1,:,:) &
          - st_x_new(imax-2,:,:)
     st_y_new(imax,:,:) = TWO*st_y_new(imax-1,:,:) &
          - st_y_new(imax-2,:,:)
     st_z_new(imax,:,:) = TWO*st_z_new(imax-1,:,:) &
          - st_z_new(imax-2,:,:)
     rho_new(:,jmax,:) = TWO*rho_new(:,jmax-1,:) &
          - rho_new(:,jmax-2,:)
     tau_new(:,jmax,:) = TWO*tau_new(:,jmax-1,:) &
          - tau_new(:,jmax-2,:)
     st_x_new(:,jmax,:) = TWO*st_x_new(:,jmax-1,:) &
          - st_x_new(:,jmax-2,:)
     st_y_new(:,jmax,:) = TWO*st_y_new(:,jmax-1,:) &
          - st_y_new(:,jmax-2,:)
     st_z_new(:,jmax,:) = TWO*st_z_new(:,jmax-1,:) &
          - st_z_new(:,jmax-2,:)
     rho_new(:,:,kmax) = TWO*rho_new(:,:,kmax-1) &
          - rho_new(:,:,kmax-2)
     tau_new(:,:,kmax) = TWO*tau_new(:,:,kmax-1) &
          - tau_new(:,:,kmax-2)
     st_x_new(:,:,kmax) = TWO*st_x_new(:,:,kmax-1) &
          - st_x_new(:,:,kmax-2)
     st_y_new(:,:,kmax) = TWO*st_y_new(:,:,kmax-1) &
          - st_y_new(:,:,kmax-2)
     st_z_new(:,:,kmax) = TWO*st_z_new(:,:,kmax-1) &
          - st_z_new(:,:,kmax-2)
     if (Symmetry .ne. OCTANT .and. Symmetry .ne. AXISYM) then
        rho_new(imin,:,:) = TWO*rho_new(imin+1,:,:) &
             - rho_new(imin+2,:,:)
        tau_new(imin,:,:) = TWO*tau_new(imin+1,:,:) &
             - tau_new(imin+2,:,:)
        st_x_new(imin,:,:) = TWO*st_x_new(imin+1,:,:) &
             - st_x_new(imin+2,:,:)
        st_y_new(imin,:,:) = TWO*st_y_new(imin+1,:,:) &
             - st_y_new(imin+2,:,:)
        st_z_new(imin,:,:) = TWO*st_z_new(imin+1,:,:) &
             - st_z_new(imin+2,:,:)
        if (Symmetry .ne. PI_SYMM) then
          rho_new(:,jmin,:) = TWO*rho_new(:,jmin+1,:) &
             - rho_new(:,jmin+2,:)
          tau_new(:,jmin,:) = TWO*tau_new(:,jmin+1,:) &
             - tau_new(:,jmin+2,:)
          st_x_new(:,jmin,:) = TWO*st_x_new(:,jmin+1,:) &
             - st_x_new(:,jmin+2,:)
          st_y_new(:,jmin,:) = TWO*st_y_new(:,jmin+1,:) &
             - st_y_new(:,jmin+2,:)
          st_z_new(:,jmin,:) = TWO*st_z_new(:,jmin+1,:) &
             - st_z_new(:,jmin+2,:)
          if (Symmetry == NO_SYMM .or. Z(1,1,kmin+1) .lt. 0.d0) then
            rho_new(:,:,kmin) = TWO*rho_new(:,:,kmin+1) &
             - rho_new(:,:,kmin+2)
            tau_new(:,:,kmin) = TWO*tau_new(:,:,kmin+1) &
             - tau_new(:,:,kmin+2)
            st_x_new(:,:,kmin) = TWO*st_x_new(:,:,kmin+1) &
             - st_x_new(:,:,kmin+2)
            st_y_new(:,:,kmin) = TWO*st_y_new(:,:,kmin+1) &
             - st_y_new(:,:,kmin+2)
            st_z_new(:,:,kmin) = TWO*st_z_new(:,:,kmin+1) &
             - st_z_new(:,:,kmin+2)
          end if
        endif
     endif
  elseif(bc==PLANAR) then
     do k=kmin+1,kmax-1
        do j=jmin,jmax
           do i=imin,imax
              if(i==imax) then
                 if(j==jmax) then
                    rho_new(i,j,k) = rho_new(i-1,j-1,k)
                    tau_new(i,j,k) = tau_new(i-1,j-1,k)
                    st_x_new(i,j,k) = st_x_new(i-1,j-1,k)
                    st_y_new(i,j,k) = st_y_new(i-1,j-1,k)
                    st_z_new(i,j,k) = st_z_new(i-1,j-1,k)
                 else if(j==jmax-1) then
                    rho_new(i,j,k) = rho_new(i-1,j,k)
                    tau_new(i,j,k) = tau_new(i-1,j,k)
                    st_x_new(i,j,k) = st_x_new(i-1,j,k)
                    st_y_new(i,j,k) = st_y_new(i-1,j,k)
                    st_z_new(i,j,k) = st_z_new(i-1,j,k)
                 else
                    rho_new(i,j,k) = rho_new(i-1,j+1,k)
                    tau_new(i,j,k) = tau_new(i-1,j+1,k)
                    st_x_new(i,j,k) = st_x_new(i-1,j+1,k)
                    st_y_new(i,j,k) = st_y_new(i-1,j+1,k)
                    st_z_new(i,j,k) = st_z_new(i-1,j+1,k)
                 end if
              else if(i==imin) then
                 if(j==jmin) then
                    rho_new(i,j,k) = rho_new(i+1,j+1,k)
                    tau_new(i,j,k) = tau_new(i+1,j+1,k)
                    st_x_new(i,j,k) = st_x_new(i+1,j+1,k)
                    st_y_new(i,j,k) = st_y_new(i+1,j+1,k)
                    st_z_new(i,j,k) = st_z_new(i+1,j+1,k)
                 else if(j==jmin+1) then
                    rho_new(i,j,k) = rho_new(i+1,j,k)
                    tau_new(i,j,k) = tau_new(i+1,j,k)
                    st_x_new(i,j,k) = st_x_new(i+1,j,k)
                    st_y_new(i,j,k) = st_y_new(i+1,j,k)
                    st_z_new(i,j,k) = st_z_new(i+1,j,k)
                 else
                    rho_new(i,j,k) = rho_new(i+1,j-1,k)
                    tau_new(i,j,k) = tau_new(i+1,j-1,k)
                    st_x_new(i,j,k) = st_x_new(i+1,j-1,k)
                    st_y_new(i,j,k) = st_y_new(i+1,j-1,k)
                    st_z_new(i,j,k) = st_z_new(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    rho_new(i,j,k) = rho_new(i,j+1,k)
                    tau_new(i,j,k) = tau_new(i,j+1,k)
                    st_x_new(i,j,k) = st_x_new(i,j+1,k)
                    st_y_new(i,j,k) = st_y_new(i,j+1,k)
                    st_z_new(i,j,k) = st_z_new(i,j+1,k)
                 else
                    rho_new(i,j,k) = rho_new(i-1,j+1,k)
                    tau_new(i,j,k) = tau_new(i-1,j+1,k)
                    st_x_new(i,j,k) = st_x_new(i-1,j+1,k)
                    st_y_new(i,j,k) = st_y_new(i-1,j+1,k)
                    st_z_new(i,j,k) = st_z_new(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    rho_new(i,j,k) = rho_new(i,j-1,k)
                    tau_new(i,j,k) = tau_new(i,j-1,k)
                    st_x_new(i,j,k) = st_x_new(i,j-1,k)
                    st_y_new(i,j,k) = st_y_new(i,j-1,k)
                    st_z_new(i,j,k) = st_z_new(i,j-1,k)
                 else
                    rho_new(i,j,k) = rho_new(i+1,j-1,k)
                    tau_new(i,j,k) = tau_new(i+1,j-1,k)
                    st_x_new(i,j,k) = st_x_new(i+1,j-1,k)
                    st_y_new(i,j,k) = st_y_new(i+1,j-1,k)
                    st_z_new(i,j,k) = st_z_new(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do
  else if (bc==COPY) then
     rho_new(imax,:,:) = rho_new(imax-1,:,:)
     tau_new(imax,:,:) = tau_new(imax-1,:,:)
     st_x_new(imax,:,:) = st_x_new(imax-1,:,:)
     st_y_new(imax,:,:) = st_y_new(imax-1,:,:)
     st_z_new(imax,:,:) = st_z_new(imax-1,:,:)
     if(Symmetry.ne.OCTANT .and. Symmetry.ne.AXISYM) then
        rho_new(imin,:,:) = rho_new(imin+1,:,:)
        tau_new(imin,:,:) = tau_new(imin+1,:,:)
        st_x_new(imin,:,:) = st_x_new(imin+1,:,:)
        st_y_new(imin,:,:) = st_y_new(imin+1,:,:)
        st_z_new(imin,:,:) = st_z_new(imin+1,:,:)
     end if
     if(Symmetry .ne. AXISYM) then
        rho_new(:,jmax,:) = rho_new(:,jmax-1,:)
        tau_new(:,jmax,:) = tau_new(:,jmax-1,:)
        st_x_new(:,jmax,:) = st_x_new(:,jmax-1,:)
        st_y_new(:,jmax,:) = st_y_new(:,jmax-1,:)
        st_z_new(:,jmax,:) = st_z_new(:,jmax-1,:)
        if(Symmetry.ne.OCTANT .and. Symmetry.ne.PI_SYMM) then
           rho_new(:,jmin,:) = rho_new(:,jmin+1,:)
           tau_new(:,jmin,:) = tau_new(:,jmin+1,:)
           st_x_new(:,jmin,:) = st_x_new(:,jmin+1,:)
           st_y_new(:,jmin,:) = st_y_new(:,jmin+1,:)
           st_z_new(:,jmin,:) = st_z_new(:,jmin+1,:)
        end if
     end if
     rho_new(:,:,kmax) = rho_new(:,:,kmax-1)
     tau_new(:,:,kmax) = tau_new(:,:,kmax-1)
     st_x_new(:,:,kmax) = st_x_new(:,:,kmax-1)
     st_y_new(:,:,kmax) = st_y_new(:,:,kmax-1)
     st_z_new(:,:,kmax) = st_z_new(:,:,kmax-1)
     if(Symmetry==NO_SYMM .or. Z(1,1,kmin+1) .lt. 0.d0) then
        rho_new(:,:,kmin) = rho_new(:,:,kmin+1)
        tau_new(:,:,kmin) = tau_new(:,:,kmin+1)
        st_x_new(:,:,kmin) = st_x_new(:,:,kmin+1)
        st_y_new(:,:,kmin) = st_y_new(:,:,kmin+1)
        st_z_new(:,:,kmin) = st_z_new(:,:,kmin+1)
     end if
  else if(bc==QUAD) then
     rho_new(imax,:,:) = THREE*rho_new(imax-1,:,:) &
          - THREE*rho_new(imax-2,:,:) + rho_new(imax-3,:,:)
     tau_new(imax,:,:) = THREE*tau_new(imax-1,:,:) &
          - THREE*tau_new(imax-2,:,:) + tau_new(imax-3,:,:)
     st_x_new(imax,:,:) = THREE*st_x_new(imax-1,:,:) &
          - THREE*st_x_new(imax-2,:,:) + st_x_new(imax-3,:,:)
     st_y_new(imax,:,:) = THREE*st_y_new(imax-1,:,:) &
          - THREE*st_y_new(imax-2,:,:) + st_y_new(imax-3,:,:)
     st_z_new(imax,:,:) = THREE*st_z_new(imax-1,:,:) &
          - THREE*st_z_new(imax-2,:,:) + st_z_new(imax-3,:,:)
     rho_new(:,jmax,:) = THREE*rho_new(:,jmax-1,:) &
          - THREE*rho_new(:,jmax-2,:) + rho_new(:,jmax-3,:)
     tau_new(:,jmax,:) = THREE*tau_new(:,jmax-1,:) &
          - THREE*tau_new(:,jmax-2,:) + tau_new(:,jmax-3,:)
     st_x_new(:,jmax,:) = THREE*st_x_new(:,jmax-1,:) &
          - THREE*st_x_new(:,jmax-2,:) + st_x_new(:,jmax-3,:)
     st_y_new(:,jmax,:) = THREE*st_y_new(:,jmax-1,:) &
          - THREE*st_y_new(:,jmax-2,:) + st_y_new(:,jmax-3,:)
     st_z_new(:,jmax,:) = THREE*st_z_new(:,jmax-1,:) &
          - THREE*st_z_new(:,jmax-2,:) + st_z_new(:,jmax-3,:)
     rho_new(:,:,kmax) = THREE*rho_new(:,:,kmax-1) &
          - THREE*rho_new(:,:,kmax-2) + rho_new(:,:,kmax-3)
     tau_new(:,:,kmax) = THREE*tau_new(:,:,kmax-1) &
          - THREE*tau_new(:,:,kmax-2) + tau_new(:,:,kmax-3)
     st_x_new(:,:,kmax) = THREE*st_x_new(:,:,kmax-1) &
          - THREE*st_x_new(:,:,kmax-2) + st_x_new(:,:,kmax-3)
     st_y_new(:,:,kmax) = THREE*st_y_new(:,:,kmax-1) &
          - THREE*st_y_new(:,:,kmax-2) + st_y_new(:,:,kmax-3)
     st_z_new(:,:,kmax) = THREE*st_z_new(:,:,kmax-1) &
          - THREE*st_z_new(:,:,kmax-2) + st_z_new(:,:,kmax-3)
     if (Symmetry .ne. OCTANT) then
      if (Symmetry .ne. AXISYM) then
        rho_new(imin,:,:) = THREE*rho_new(imin+1,:,:) &
             - THREE*rho_new(imin+2,:,:) + rho_new(imin+3,:,:)
        tau_new(imin,:,:) = THREE*tau_new(imin+1,:,:) &
             - THREE*tau_new(imin+2,:,:) + tau_new(imin+3,:,:)
        st_x_new(imin,:,:) = THREE*st_x_new(imin+1,:,:) &
             - THREE*st_x_new(imin+2,:,:) + st_x_new(imin+3,:,:)
        st_y_new(imin,:,:) = THREE*st_y_new(imin+1,:,:) &
             - THREE*st_y_new(imin+2,:,:) + st_y_new(imin+3,:,:)
        st_z_new(imin,:,:) = THREE*st_z_new(imin+1,:,:) &
             - THREE*st_z_new(imin+2,:,:) + st_z_new(imin+3,:,:)
       end if
        if (Symmetry .ne. PI_SYMM) then
           rho_new(:,jmin,:) = THREE*rho_new(:,jmin+1,:) &
                - THREE*rho_new(:,jmin+2,:) + rho_new(:,jmin+3,:)
           tau_new(:,jmin,:) = THREE*tau_new(:,jmin+1,:) &
                - THREE*tau_new(:,jmin+2,:) + tau_new(:,jmin+3,:)
           st_x_new(:,jmin,:) = THREE*st_x_new(:,jmin+1,:) &
                - THREE*st_x_new(:,jmin+2,:) + st_x_new(:,jmin+3,:)
           st_y_new(:,jmin,:) = THREE*st_y_new(:,jmin+1,:) &
                - THREE*st_y_new(:,jmin+2,:) + st_y_new(:,jmin+3,:)
           st_z_new(:,jmin,:) = THREE*st_z_new(:,jmin+1,:) &
                - THREE*st_z_new(:,jmin+2,:) + st_z_new(:,jmin+3,:)
           if (Symmetry == NO_SYMM .or. Z(1,1,kmin+1) .lt. 0.d0) then
              rho_new(:,:,kmin) = THREE*rho_new(:,:,kmin+1) &
                   - THREE*rho_new(:,:,kmin+2) + rho_new(:,:,kmin+3)
              tau_new(:,:,kmin) = THREE*tau_new(:,:,kmin+1) &
                   - THREE*tau_new(:,:,kmin+2) + tau_new(:,:,kmin+3)
              st_x_new(:,:,kmin) = THREE*st_x_new(:,:,kmin+1) &
                   - THREE*st_x_new(:,:,kmin+2) + st_x_new(:,:,kmin+3)
              st_y_new(:,:,kmin) = THREE*st_y_new(:,:,kmin+1) &
                   - THREE*st_y_new(:,:,kmin+2) + st_y_new(:,:,kmin+3)
              st_z_new(:,:,kmin) = THREE*st_z_new(:,:,kmin+1) &
                   - THREE*st_z_new(:,:,kmin+2) + st_z_new(:,:,kmin+3)
           end if
        endif
     endif
  end if
end subroutine hydro_matter_bc
!-----------------------------------------------------------------------------
!
! Fix fluid nablas on outer boundaries
!
!-----------------------------------------------------------------------------
subroutine nablas_bc(ex,nrho,nP,nvx,nvy,nvz,nBx,nBy,nBz,Symmetry)
  implicit none
  integer, dimension(3)                :: ex
  real*8, dimension(ex(1),ex(2),ex(3)) :: nrho,nP,nvx,nvy,nvz
  real*8, dimension(ex(1),ex(2),ex(3)) :: nBx,nBy,nBz
  integer                              :: Symmetry
  integer                            :: bc, EXTRAP, PLANAR, COPY, QUAD
  parameter(COPY = 1, EXTRAP = 2, QUAD = 3, PLANAR = 4)
  interface
     subroutine one_nabla_bc(ex,f_new,Symmetry,bc)
       implicit none
       integer, dimension(3)                 :: ex
       real*8, dimension(ex(1),ex(2),ex(3))  :: f_new
       integer                               :: Symmetry,bc
     end subroutine one_nabla_bc
  end interface
!  
! Choose boundary condition
!
  bc = COPY
!
! Apply to each nabla
!  
  call one_nabla_bc(ex,nrho,Symmetry,bc)
  call one_nabla_bc(ex,nP,Symmetry,bc)
  call one_nabla_bc(ex,nvx,Symmetry,bc)
  call one_nabla_bc(ex,nvy,Symmetry,bc)
  call one_nabla_bc(ex,nvz,Symmetry,bc)
  call one_nabla_bc(ex,nBx,Symmetry,bc)
  call one_nabla_bc(ex,nBy,Symmetry,bc)
  call one_nabla_bc(ex,nBz,Symmetry,bc)
end subroutine nablas_bc
subroutine one_nabla_bc(ex,f_new,Symmetry,bc)
  implicit none
  integer, dimension(3)                 :: ex
  real*8, dimension(ex(1),ex(2),ex(3))  :: f_new
  integer                               :: Symmetry,bc
  integer                            :: i,j,k
  integer                            :: imin, jmin, kmin, imax, jmax, kmax
  real*8, parameter                  :: ZERO = 0.D0
  real*8, parameter                  :: TWO = 2.D0, THREE = 3.d0
  real*8                             :: nx,ny,nz
  integer                            :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM
  integer                            :: AXISYM
  integer                            :: EXTRAP, PLANAR, COPY, QUAD
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  parameter(COPY = 1, EXTRAP = 2, QUAD = 3, PLANAR = 4)
  imin = lbound(f_new,1)
  jmin = lbound(f_new,2)
  kmin = lbound(f_new,3)
  imax = ubound(f_new,1)
  jmax = ubound(f_new,2)
  kmax = ubound(f_new,3)
  if(bc==EXTRAP) then
     f_new(imax,:,:) = TWO*f_new(imax-1,:,:) &
          - f_new(imax-2,:,:)
     f_new(:,jmax,:) = TWO*f_new(:,jmax-1,:) &
          - f_new(:,jmax-2,:)
     f_new(:,:,kmax) = TWO*f_new(:,:,kmax-1) &
          - f_new(:,:,kmax-2)
     if (Symmetry .ne. OCTANT .and. Symmetry .ne. AXISYM) then
        f_new(imin,:,:) = TWO*f_new(imin+1,:,:) &
             - f_new(imin+2,:,:)
        if (Symmetry .ne. PI_SYMM) then
          f_new(:,jmin,:) = TWO*f_new(:,jmin+1,:) &
             - f_new(:,jmin+2,:)
          if (Symmetry == NO_SYMM) then
            f_new(:,:,kmin) = TWO*f_new(:,:,kmin+1) &
             - f_new(:,:,kmin+2)
          end if
        endif
     endif
  elseif(bc==PLANAR) then
     do k=kmin+1,kmax-1
        do j=jmin,jmax
           do i=imin,imax
              if(i==imax) then
                 if(j==jmax) then
                    f_new(i,j,k) = f_new(i-1,j-1,k)
                 else if(j==jmax-1) then
                    f_new(i,j,k) = f_new(i-1,j,k)
                 else
                    f_new(i,j,k) = f_new(i-1,j+1,k)
                 end if
              else if(i==imin) then
                 if(j==jmin) then
                    f_new(i,j,k) = f_new(i+1,j+1,k)
                 else if(j==jmin+1) then
                    f_new(i,j,k) = f_new(i+1,j,k)
                 else
                    f_new(i,j,k) = f_new(i+1,j-1,k)
                 end if
              else if(j==jmin) then
                 if(i==imin+1) then
                    f_new(i,j,k) = f_new(i,j+1,k)
                 else
                    f_new(i,j,k) = f_new(i-1,j+1,k)
                 end if
              else if(j==jmax) then
                 if(i==imax-1) then
                    f_new(i,j,k) = f_new(i,j-1,k)
                 else
                    f_new(i,j,k) = f_new(i+1,j-1,k)
                 end if
              end if
           end do
        end do
     end do
  else if (bc==COPY) then
     f_new(imin,:,:) = f_new(imin+1,:,:)
     f_new(imax,:,:) = f_new(imax-1,:,:)
     f_new(:,jmin,:) = f_new(:,jmin+1,:)
     f_new(:,jmax,:) = f_new(:,jmax-1,:)
     f_new(:,:,kmin) = f_new(:,:,kmin+1)
     f_new(:,:,kmax) = f_new(:,:,kmax-1)
  else if(bc==QUAD) then
     f_new(imax,:,:) = THREE*f_new(imax-1,:,:) &
          - THREE*f_new(imax-2,:,:) + f_new(imax-3,:,:)
     f_new(:,jmax,:) = THREE*f_new(:,jmax-1,:) &
          - THREE*f_new(:,jmax-2,:) + f_new(:,jmax-3,:)
     f_new(:,:,kmax) = THREE*f_new(:,:,kmax-1) &
          - THREE*f_new(:,:,kmax-2) + f_new(:,:,kmax-3)
     if (Symmetry .ne. OCTANT .and. Symmetry .ne. AXISYM) then
        f_new(imin,:,:) = THREE*f_new(imin+1,:,:) &
             - THREE*f_new(imin+2,:,:) + f_new(imin+3,:,:)
        if (Symmetry .ne. PI_SYMM) then
          f_new(:,jmin,:) = THREE*f_new(:,jmin+1,:) &
             - THREE*f_new(:,jmin+2,:) + f_new(:,jmin+3,:)
          if (Symmetry == NO_SYMM) then
            f_new(:,:,kmin) = THREE*f_new(:,:,kmin+1) &
             - THREE*f_new(:,:,kmin+2) + f_new(:,:,kmin+3)
          end if
        endif
     endif
  end if
end subroutine one_nabla_bc
! ----------------------------------------------------------------------
!                                                                      !
!    Compute K = P/rho_b^Gamma                                         !
!                                                                      !
! ---------------------------------------------------------------------!
!
subroutine polyk(ex,k,P,rho_b,rho_b_atm,K_init,n)
  implicit none
  integer, dimension(3)                       :: ex
  real*8, dimension(ex(1),ex(2),ex(3))        :: k,P,rho_b
  real*8                                      :: n, gamma, K_init,rho_b_atm
!
  gamma=1.d0+1.d0/n
  where (rho_b > 1.0001d0*rho_b_atm)
    k=P/rho_b**gamma
  elsewhere
    k=K_init
  end where
end subroutine polyk

subroutine compute_optical_depth(ext, X,Y,Z, dx, dy, dz, optd, optd_x, optd_y, optd_z, ka_gf, ks_gf)

  implicit none

  integer, dimension(3)                   :: ext
  real*8, dimension(ext(1),ext(2),ext(3)) :: X,Y,Z
  real*8, dimension(ext(1),ext(2),ext(3)) :: optd, optd_x, optd_y, optd_z, ka_gf, ks_gf
  real*8                                  :: dx,dy,dz
  integer                                 :: i,j,k, ii,jj,kk

print *, "^^^^^^^^^^^^^^^^^^^^^^ Start compute_optical_depth ^^^^^^^^^^^^^^^^^^^^^^^^"


  do k=1, ext(3)
     do j=1, ext(2)
        do i=1, ext(1)

!!!Initialize optd_x,y,z                                                                                                                                                                                                                                                   
           optd_x(i,j,k)=0.d0
           optd_y(i,j,k)=0.d0
           optd_z(i,j,k)=0.d0


           if (X(i,j,k) .ge. 0.0) then
              if (Y(i,j,k) .ge. 0.0) then
                 if (Z(i,j,k) .ge. 0.0) then

!!! 1st octant X>=0, Y>=0, Z>=0                                                                                                                                                                                                                                            
                    do ii = i, ext(1)
                       optd_x(i,j,k) = optd_x(i,j,k) + (ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                    end do
                    do jj = j, ext(2)
                       optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                    end do
                    do kk = k, ext(3)
                       optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                    end do
!!!-------------------------------------------------                                                                                                                                                                                                                       

                  else if (Z(i,j,k) .lt. 0.0) then
!!! 2nd octant X>=0, Y>=0, Z<0                                                                                                                                                                                                                                             
                     do ii = i,ext(1)
                        optd_x(i,j,k) = optd_x(i,j,k) + (ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = j,ext(2)
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = 1,k
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  end if

               else if (Y(i,j,k) .lt. 0.0) then

                  if (Z(i,j,k) .ge. 0.0) then
!!! 3rd octant X>=0, Y<0, Z>=0                                                                                                                                                                                                                                             
                     do ii = i,ext(1)
                        optd_x(i,j,k) = optd_x(i,j,k) +(ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = 1,j
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = k,ext(3)
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  else if (Z(i,j,k) .lt. 0.0) then
!!! 4nd octant X>=0, Y<0, Z<0                                                                                                                                                                                                                                              
                     do ii = i,ext(1)
                        optd_x(i,j,k) = optd_x(i,j,k) +(ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = 1,j
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = 1,k
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  end if

               end if

            else if (X(i,j,k) .lt. 0.0) then
               if (Y(i,j,k) .ge. 0.0) then
                  if (Z(i,j,k) .ge. 0.0) then
!!! 5th octant X<0, Y>=0, Z>=0                                                                                                                                                                                                                                             
                     do ii = 1,i
                        optd_x(i,j,k) = optd_x(i,j,k) + (ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = j, ext(2)
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = k, ext(3)
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  else if (Z(i,j,k) .lt. 0.0) then
!!! 6th octant X<0, Y>=0, Z<0                                                                                                                                                                                                                                              
                     do ii = 1,i
                        optd_x(i,j,k) = optd_x(i,j,k) + (ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = j, ext(2)
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = 1,k
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  end if

               else if (Y(i,j,k) .lt. 0.0) then

                  if (Z(i,j,k) .ge. 0.0) then
!!! 7th octant X<0, Y<0, Z>=0                                                                                                                                                                                                                                              
                     do ii = 1,i
                        optd_x(i,j,k) = optd_x(i,j,k) + (ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = 1,j
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = k,ext(3)
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  else if (Z(i,j,k) .lt. 0.0) then
!!! 8th octant X<0, Y<0, Z<0                                                                                                                                                                                                                                               
                     do ii = 1,i
                        optd_x(i,j,k) = optd_x(i,j,k) + (ka_gf(ii,j,k)+ks_gf(ii,j,k))*dx
                     end do
                     do jj = 1,j
                        optd_y(i,j,k) = optd_y(i,j,k) + (ka_gf(i,jj,k)+ks_gf(i,jj,k))*dy
                     end do
                     do kk = 1,k
                        optd_z(i,j,k) = optd_z(i,j,k) + (ka_gf(i,j,kk)+ks_gf(i,j,kk))*dz
                     end do
!!!-------------------------------------------------                                                                                                                                                                                                                       
                  end if

               end if
            end if



            optd(i,j,k) = min(optd_x(i,j,k), optd_y(i,j,k), optd_z(i,j,k))

        end do
     end do
  end do




 print *, "^^^^^^^^^^^^^^^^^^^^^^ Done with compute_optical_depth ^^^^^^^^^^^^^^^^^^^^^^^^"

end subroutine compute_optical_depth

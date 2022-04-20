c     $Header: /numrelcvs/AEIThorns/Dissipation/src/apply_dissipation.F77,v 1.8 2007/12/05 23:14:13 schnetter Exp $

#include "cctk.h"

      subroutine apply_dissipation (var, rhs, ni, nj, nk, dx, order, epsdis)
      implicit none
      
      integer   ni, nj, nk
      CCTK_REAL var(ni,nj,nk), rhs(ni,nj,nk)
      CCTK_REAL dx(3)
      CCTK_INT  order
      CCTK_REAL epsdis(ni,nj,nk)
      
      integer   i, j, k
      
      if (order .eq. 1) then
         
c$omp parallel do
         do k = 2, nk-1
            do j = 2, nj-1
               do i = 2, ni-1
                  rhs(i,j,k) = rhs(i,j,k) + epsdis(i,j,k)
     $                 * (+ (var(i-1,j,k) - 2*var(i,j,k) + var(i+1,j,k)) / dx(1)
     $                    + (var(i,j-1,k) - 2*var(i,j,k) + var(i,j+1,k)) / dx(2)
     $                    + (var(i,j,k-1) - 2*var(i,j,k) + var(i,j,k+1)) / dx(3))
                  
               end do
            end do
         end do
c$omp end parallel do
         
      else if (order .eq. 3) then
         
c$omp parallel do
         do k = 3, nk-2
            do j = 3, nj-2
               do i = 3, ni-2
                  
                  rhs(i,j,k) = rhs(i,j,k) - epsdis(i,j,k) / 16
     $                 * (+ (var(i-2,j,k) - 4*var(i-1,j,k) + 6*var(i,j,k) - 4*var(i+1,j,k) + var(i+2,j,k)) / dx(1)
     $                    + (var(i,j-2,k) - 4*var(i,j-1,k) + 6*var(i,j,k) - 4*var(i,j+1,k) + var(i,j+2,k)) / dx(2)
     $                    + (var(i,j,k-2) - 4*var(i,j,k-1) + 6*var(i,j,k) - 4*var(i,j,k+1) + var(i,j,k+2)) / dx(3))
                  
               end do
            end do
         end do
c$omp end parallel do
         
      else if (order .eq. 5) then
         
c$omp parallel do
        do k = 4, nk-3
          do j = 4, nj-3
            do i = 4, ni-3
               
              rhs(i,j,k) = rhs(i,j,k) + epsdis(i,j,k) / 64
     $           * (+ (var(i-3,j,k) - 6*var(i-2,j,k) + 15*var(i-1,j,k) - 20*var(i,j,k) + 15*var(i+1,j,k) - 6*var(i+2,j,k) + var(i+3,j,k)) / dx(1)
     $              + (var(i,j-3,k) - 6*var(i,j-2,k) + 15*var(i,j-1,k) - 20*var(i,j,k) + 15*var(i,j+1,k) - 6*var(i,j+2,k) + var(i,j+3,k) ) / dx(2)
     $              + (var(i,j,k-3) - 6*var(i,j,k-2) + 15*var(i,j,k-1) - 20*var(i,j,k) + 15*var(i,j,k+1) - 6*var(i,j,k+2) + var(i,j,k+3) ) / dx(3))
              
            end do
          end do
        end do
c$omp end parallel do
        
      else if (order .eq. 7) then
         
c$omp parallel do
        do k = 5, nk-4
          do j = 5, nj-4
            do i = 5, ni-4
               
              rhs(i,j,k) = rhs(i,j,k) - epsdis(i,j,k) / 256
     $           * (+ (var(i-4,j,k) - 8*var(i-3,j,k) + 28*var(i-2,j,k) - 56*var(i-1,j,k) + 70*var(i,j,k) - 56*var(i+1,j,k) + 28*var(i+2,j,k) - 8*var(i+3,j,k) + var(i+4,j,k)) / dx(1)
     $              + (var(i,j-4,k) - 8*var(i,j-3,k) + 28*var(i,j-2,k) - 56*var(i,j-1,k) + 70*var(i,j,k) - 56*var(i,j+1,k) + 28*var(i,j+2,k) - 8*var(i,j+3,k) + var(i,j+4,k)) / dx(2)
     $              + (var(i,j,k-4) - 8*var(i,j,k-3) + 28*var(i,j,k-2) - 56*var(i,j,k-1) + 70*var(i,j,k) - 56*var(i,j,k+1) + 28*var(i,j,k+2) - 8*var(i,j,k+3) + var(i,j,k+4)) / dx(3))
              
            end do
          end do
        end do
c$omp end parallel do
        
      else
         
         call CCTK_WARN (0, "internal error")
         
      end if
      
      end

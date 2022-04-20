#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
!
subroutine scalarwaveMoL_postregrid(CCTK_ARGUMENTS)

  implicit none

  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer :: levelnumber,i,j,k

  levelnumber = int( log(cctk_levfac(1)*1.d0)/log(2.D0)+1.D0+1.d-10 )

  do k=1,cctk_lsh(3)
     do j=1,cctk_lsh(2)
        do i=1,cctk_lsh(1)
	   if ( abs(X(i,j,k))+abs(Y(i,j,k)+1.8d0)+abs(Z(i,j,k)) .lt. 1.d-10) then
	      write(*,*) 'Scalar wave post regrid: '
	      write(*,*) levelnumber,i,j,k,phi_stagger(i,j,k)
	   end if
	end do
     end do 
  end do

end subroutine scalarwaveMoL_postregrid

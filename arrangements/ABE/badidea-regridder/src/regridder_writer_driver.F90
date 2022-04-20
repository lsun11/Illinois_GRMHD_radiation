
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine bhns_regridder_writer_driver(CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  integer               :: current_line,which_processor,checksum
  real*8 		:: xmin_interp,ymin_interp,zmin_interp,dx_interp,dy_interp,dz_interp
  integer               :: Nx_interp,Ny_interp,Nz_interp

  integer               :: i,j,k
  integer               :: gridfunc_varindex,numpoints,which_point
  real*8, allocatable, dimension(:,:)          :: pointcoords
  real*8, allocatable, dimension(:)            :: output

  do current_line=1,10000
     
     ! Read the grid parameters
     call bhns_regrid_read_one_line_of_inputfile(current_line, &
          xmin_interp,ymin_interp,zmin_interp, &
          dx_interp,dy_interp,dz_interp, &
          Nx_interp,Ny_interp,Nz_interp, &
          which_processor,checksum)
     
     if(which_processor .ne. -1) then
        numpoints = Nx_interp*Ny_interp*Nz_interp
        allocate(output(numpoints))
        allocate(pointcoords(numpoints,3))

        which_point=0
        do k=1,Nz_interp
           do j=1,Ny_interp
              do i=1,Nx_interp
                 pointcoords(which_point,1) = xmin_interp + i*dx_interp
                 pointcoords(which_point,2) = ymin_interp + j*dy_interp
                 pointcoords(which_point,3) = zmin_interp + k*dz_interp
                 which_point = which_point + 1
              end do
           end do
        end do
        
        ! Magnetic fields
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::Ax");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::Ay");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::Az");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::psi6phi");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)

        ! Hydro
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::rho_star");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::tau");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::mhd_st_x");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::mhd_st_y");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"mhd_evolve::mhd_st_z");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)

        ! Metric
        call CCTK_VarIndex(gridfunc_varindex,"bssn::phi");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::chi");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::trK");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)

        call CCTK_VarIndex(gridfunc_varindex,"bssn::gxx");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::gxy");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::gxz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::gyy");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::gyz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::gzz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)

        call CCTK_VarIndex(gridfunc_varindex,"bssn::Axx");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Axy");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Axz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Ayy");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Ayz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Azz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)

        call CCTK_VarIndex(gridfunc_varindex,"bssn::Gammax");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Gammay");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        call CCTK_VarIndex(gridfunc_varindex,"bssn::Gammaz");
        call interp_driver_carp(cctkGH,numpoints,pointcoords,gridfunc_varindex,output); call bhns_regridder_out_1gf(numpoints,output,checksum)
        
        deallocate(pointcoords)
        deallocate(output)
     end if

  end do

end subroutine bhns_regridder_writer_driver

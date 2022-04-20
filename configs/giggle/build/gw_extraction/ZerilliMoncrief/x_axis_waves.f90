!----------------------------------------------------------------------------- 
!
! Contains wave extraction algorithms (uses MPI extensively)
!
!----------------------------------------------------------------------------- 
!===================================================================
! Computes h+ and hx on the upper x-axis boundary assuming TT gauge
!===================================================================
subroutine x_axis_waves(cctkGH,cctk_nghostzones,xout,ext,global_ext, &
     dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglobal, &
     phi,gxx,gxy,gxz,gyy,gyz,gzz,PhysR,dRdr,amplitudes)
  implicit none
  REAL*8  excision_radius
  REAL*8  run_time
  INTEGER*4 Symmetry
  INTEGER*4 bssn_enable
  INTEGER*4 cowling_enable
  INTEGER*4 excision_enable
  INTEGER*4 fisheye_enable
  INTEGER*4 iter_count
  INTEGER*4 number_of_mol_ministeps
  INTEGER*4 rot_metric
  INTEGER*4 trA_detg_enforce
  COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count&
  &,number_of_mol_ministeps,rot_metric,trA_detg_enforce
  REAL*8  dR_GW
  REAL*8  ddR_GW
  REAL*8  phi_GW
  REAL*8  radius_GW
  REAL*8  radius_GW_phys
  REAL*8  theta_GW
  COMMON /gw_extractionrest/dR_GW,ddR_GW,phi_GW,radius_GW,radius_GW_phys,theta_GW
  REAL*8  compute_Psi4_max_radius
  REAL*8  compute_Psi4_min_radius
  REAL*8  radius_GW_Psi4(101)
  integer*8  psif_vec
  INTEGER*4 enable_interp_onepoint
  INTEGER*4 nmodes_Psi4
  INTEGER*4 nmodes_ZM
  INTEGER*4 num_extraction_radii
  INTEGER*4 radius_power
  INTEGER*4 scale_with_radius
  INTEGER*4 use_Rij_from_compute_ricci
  COMMON /gw_extractionpriv/compute_Psi4_max_radius,compute_Psi4_min_radius,radius_GW_Psi4,psif_vec,enable_interp_onepoint,nmodes_P&
  &si4,nmodes_ZM,num_extraction_radii,radius_power,scale_with_radius,use_Rij_from_compute_ricci
  REAL*8  CCTKH0
  REAL*8  CCTKH6
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  CCTKH9
  REAL*8  CCTKH13
  REAL*8  rbr
  REAL*8  CCTKH17
  REAL*8  surf_radius
  INTEGER*4 CCTKH1
  INTEGER*4 CCTKH2
  INTEGER*4 CCTKH3
  INTEGER*4 CCTKH4
  INTEGER*4 CCTKH5
  INTEGER*4 CCTKH10
  INTEGER*4 CCTKH11
  INTEGER*4 CCTKH12
  INTEGER*4 CCTKH14
  INTEGER*4 CCTKH15
  INTEGER*4 CCTKH16
  INTEGER*4 numphi
  INTEGER*4 numtheta
  INTEGER*4 CCTKH18
  INTEGER*4 CCTKH19
  COMMON /DIAGNOSTICS_VACUUMrest/CCTKH0,CCTKH6,CCTKH7,CCTKH8,CCTKH9,CCTKH13,rbr,CCTKH17,surf_radius,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCT&
  &KH5,CCTKH10,CCTKH11,CCTKH12,CCTKH14,CCTKH15,CCTKH16,numphi,numtheta,CCTKH18,CCTKH19
  REAL*8  CCTKH25
  REAL*8  CCTKH39
  REAL*8  CCTKH47
  REAL*8  CCTKH49
  REAL*8  CCTKH51
  REAL*8  CCTKH53
  REAL*8  CCTKH55
  REAL*8  CCTKH57
  REAL*8  CCTKH59
  REAL*8  CCTKH61
  REAL*8  CCTKH63
  integer*8  CCTKH22
  integer*8  CCTKH23
  integer*8  CCTKH26
  integer*8  CCTKH29
  integer*8  CCTKH30
  integer*8  CCTKH31
  integer*8  CCTKH34
  integer*8  CCTKH35
  integer*8  CCTKH40
  integer*8  CCTKH41
  integer*8  CCTKH43
  integer*8  CCTKH65
  integer*8  CCTKH67
  integer*8  CCTKH69
  integer*8  CCTKH71
  integer*8  CCTKH72
  integer*8  CCTKH75
  INTEGER*4 CCTKH20
  INTEGER*4 CCTKH21
  INTEGER*4 CCTKH24
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH32
  INTEGER*4 CCTKH33
  INTEGER*4 CCTKH36
  INTEGER*4 CCTKH37
  INTEGER*4 CCTKH38
  INTEGER*4 out_every
  INTEGER*4 CCTKH42
  INTEGER*4 CCTKH44
  INTEGER*4 CCTKH45
  INTEGER*4 CCTKH46
  INTEGER*4 CCTKH48
  INTEGER*4 CCTKH50
  INTEGER*4 CCTKH52
  INTEGER*4 CCTKH54
  INTEGER*4 CCTKH56
  INTEGER*4 CCTKH58
  INTEGER*4 CCTKH60
  INTEGER*4 CCTKH62
  INTEGER*4 CCTKH64
  INTEGER*4 CCTKH66
  INTEGER*4 CCTKH68
  INTEGER*4 CCTKH70
  INTEGER*4 CCTKH73
  INTEGER*4 CCTKH74
  COMMON /IOrest/CCTKH25,CCTKH39,CCTKH47,CCTKH49,CCTKH51,CCTKH53,CCTKH55,CCTKH57,CCTKH59,CCTKH61,CCTKH63,CCTKH22,CCTKH23,CCTKH26,CC&
  &TKH29,CCTKH30,CCTKH31,CCTKH34,CCTKH35,CCTKH40,CCTKH41,CCTKH43,CCTKH65,CCTKH67,CCTKH69,CCTKH71,CCTKH72,CCTKH75,CCTKH20,CCTKH21,CC&
  &TKH24,CCTKH27,CCTKH28,CCTKH32,CCTKH33,CCTKH36,CCTKH37,CCTKH38,out_every,CCTKH42,CCTKH44,CCTKH45,CCTKH46,CCTKH48,CCTKH50,CCTKH52,&
  &CCTKH54,CCTKH56,CCTKH58,CCTKH60,CCTKH62,CCTKH64,CCTKH66,CCTKH68,CCTKH70,CCTKH73,CCTKH74
  
  external     CCTK_PointerTo
  integer*8 CCTK_PointerTo
  interface
  integer function CCTK_Equals (arg1, arg2)
  implicit none
  integer*8 arg1
  character(*) arg2
  end function CCTK_Equals
  integer function CCTK_MyProc (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_MyProc
  integer function CCTK_nProcs (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_nProcs
  integer function CCTK_IsThornActive (name)
  implicit none
  character(*) name
  end function CCTK_IsThornActive
  integer*8 function CCTK_NullPointer ()
  implicit none
  end function CCTK_NullPointer
  end interface
  
  integer*8 :: cctkgh
  integer, dimension(3)                    :: ext,cctk_nghostzones,global_ext
  integer                                  :: interpolation_order
  real*8,dimension(ext(1))                 :: Xlocal1d
  real*8,dimension(ext(2))                 :: Ylocal1d
  real*8,dimension(ext(3))                 :: Zlocal1d
  real*8,dimension(ext(1),ext(2),ext(3))   :: phi,gxx,gxy,gxz,gyy,gyz,gzz
  real*8,dimension(ext(1),ext(2),ext(3))   :: PhysR,dRdr
  real*8,dimension(global_ext(1))          :: Xglobal
  real*8,dimension(global_ext(2))          :: Yglobal
  real*8,dimension(global_ext(3))          :: Zglobal
  real*8                                   :: xout,dX,dY,dZ
  real*8,dimension(16,3)                   :: pointcoords
  real*8,dimension(2)                      :: amplitudes
  real*8,dimension(16)                     :: gzz_f,gyy_f,gyz_f
  real*8,dimension(16)                     :: gxx_f,gxy_f,gxz_f
  real*8,dimension(16)                     :: fishr_f,PhysR_f,dRdr_f
  real*8,dimension(16)                     :: x_f,y_f,z_f
  real*8,dimension(16)                     :: fac,Jxx,Jxy,Jxz,Jyy,Jyz,Jzz
  real*8,dimension(16)                     :: gxxi,gxyi,gxzi,gyyi,gyzi,gzzi
  real*8,dimension(4,4)                    :: hzz,hyy,hyz,haveit,Rp
  real*8                                   :: symy,symz,r
  real*8                                   :: hzz_ll, hzz_l, hzz_u, hzz_uu
  real*8                                   :: hyz_ll, hyz_l, hyz_u, hyz_uu
  real*8                                   :: hyy_ll, hyy_l, hyy_u, hyy_uu
  real*8                                   :: hzz_cent, hyz_cent, hyy_cent
  real*8                                   :: r_ll, r_l, r_u, r_uu
  real*8                                   :: h_plus, h_cross
  integer                                  :: i,j
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  haveit = 0.D0
  hzz = 0.D0
  hyy = 0.D0
  hyz = 0.D0
!  pointcoords(:,1) = xout
  pointcoords(:,1) = Xglobal((xout-Xglobal(1))/dX+1)
  if(Symmetry==OCTANT) then
     do i=1,4
        do j=1,4
           pointcoords((i-1)*4 + j,2) = -1.5D0*dY + (i-1)*dY
           pointcoords((i-1)*4 + j,3) = -1.5D0*dZ + (j-1)*dZ
           if(abs(pointcoords((i-1)*4 + j,2)+0.5D0*dY).lt.1.D-5*dY) then
              pointcoords((i-1)*4 + j,2) = 0.5D0*dY
           else if(abs(pointcoords((i-1)*4 + j,2)+1.5D0*dY).lt.1.D-5*dY) then
              pointcoords((i-1)*4 + j,2) = 1.5D0*dY
           end if
           if(abs(pointcoords((i-1)*4 + j,3)+0.5D0*dZ).lt.1.D-5*dZ) then
              pointcoords((i-1)*4 + j,3) = 0.5D0*dZ
           else if(abs(pointcoords((i-1)*4 + j,3)+1.5D0*dZ).lt.1.D-5*dZ) then
              pointcoords((i-1)*4 + j,3) = 1.5D0*dZ
           end if
        end do
     end do
  else if(Symmetry == AXISYM) then
     do i=1,4
        do j=1,4
           pointcoords((i-1)*4 + j,2) = 0.D0
           pointcoords((i-1)*4 + j,3) = -1.5D0*dZ + (j-1)*dZ
           if(abs(pointcoords((i-1)*4 + j,3)+0.5D0*dZ).lt.1.D-5*dZ) then
              pointcoords((i-1)*4 + j,3) = 0.5D0*dZ
           else if(abs(pointcoords((i-1)*4 + j,3)+1.5D0*dZ).lt.1.D-5*dZ) then
              pointcoords((i-1)*4 + j,3) = 1.5D0*dZ
           end if
        end do
     end do
  else if (Symmetry == EQUATORIAL) then
     do i=1,4
        do j=1,4
           pointcoords((i-1)*4 + j,2) = -1.5D0*dY + (i-1)*dY
           pointcoords((i-1)*4 + j,3) = -1.5D0*dZ + (j-1)*dZ
           if(abs(pointcoords((i-1)*4 + j,3)+0.5D0*dZ).lt.1.D-5*dZ) then
              pointcoords((i-1)*4 + j,3) = 0.5D0*dZ
           else if(abs(pointcoords((i-1)*4 + j,3)+1.5D0*dZ).lt.1.D-5*dZ) then
              pointcoords((i-1)*4 + j,3) = 1.5D0*dZ
           end if
        end do
     end do
  else
     do i=1,4
        do j=1,4
           pointcoords((i-1)*4 + j,2) = -1.5D0*dY + (i-1)*dY
           pointcoords((i-1)*4 + j,3) = -1.5D0*dZ + (j-1)*dZ
        end do
     end do
  end if
  call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglo&
  &bal,pointcoords,gyy,gyy_f)
  call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglo&
  &bal,pointcoords,gyz,gyz_f)
  call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Zglo&
  &bal,pointcoords,gzz,gzz_f)
  ! Fisheye stuff
  if (fisheye_enable==1) then
     x_f(:) = pointcoords(:,1)
     y_f(:) = pointcoords(:,2)
     z_f(:) = pointcoords(:,3)
     fishr_f = sqrt(x_f*x_f + y_f*y_f + z_f*z_f)
     call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Z&
  &global,pointcoords,gxx,gxx_f)
     call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Z&
  &global,pointcoords,gxy,gxy_f)
     call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Z&
  &global,pointcoords,gxz,gxz_f)
     call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Z&
  &global,pointcoords,PhysR,PhysR_f)
     call return_funcvals_at_points(cctkGH,cctk_nghostzones,16,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Yglobal,Z&
  &global,pointcoords,dRdr,dRdr_f)
     ! Transform gyy, gyz and gzz
     fac = dRdr_f - PhysR_f/fishr_f
     Jxx = PhysR_f/fishr_f + x_f/fishr_f * x_f/fishr_f * fac
     Jxy = x_f/fishr_f * y_f/fishr_f * fac
     Jxz = x_f/fishr_f * z_f/fishr_f * fac
     Jyy = PhysR_f/fishr_f + y_f/fishr_f * y_f/fishr_f * fac
     Jyz = y_f/fishr_f * z_f/fishr_f * fac
     Jzz = PhysR_f/fishr_f + z_f/fishr_f * z_f/fishr_f * fac
     fac = (PhysR_f/fishr_f)**(-4.d0/3.d0) * dRdr_f**(-2.d0/3.d0)
     gxxi = gxx_f
     gxyi = gxy_f
     gxzi = gxz_f
     gyyi = gyy_f
     gyzi = gyz_f
     gzzi = gzz_f
     gyy_f = (Jxy*Jxy*gxxi + 2.d0*Jxy*Jyy*gxyi + &
                2.d0*Jxy*Jyz*gxzi + Jyy*Jyy*gyyi + &
                2.d0*Jyy*Jyz*gyzi + Jyz*Jyz*gzzi)*fac
     gyz_f = (Jxy*Jxz*gxxi + (Jxy*Jyz+Jyy*Jxz)*gxyi + &
                (Jxy*Jzz+Jyz*Jxz)*gxzi + Jyy*Jyz*gyyi + &
                (Jyy*Jzz+Jyz*Jyz)*gyzi + Jyz*Jzz*gzzi)*fac
     gzz_f = (Jxz*Jxz*gxxi + 2.d0*Jxz*Jyz*gxyi + &
                2.d0*Jxz*Jzz*gxzi + Jyz*Jyz*gyyi + &
                2.d0*Jyz*Jzz*gyzi + Jzz*Jzz*gzzi)*fac
   end if
  do i=1,4
     do j=1,4
        if (Symmetry == OCTANT .or. Symmetry == PI_SYMM) then
           if(i.le.2) then
              symy = -1.D0
           else
              symy = 1.D0
           end if
           if(j.le.2) then
              symz = -1.D0
           else
              symz = 1.D0
           end if
        else if (Symmetry == EQUATORIAL) then
           symy = 1.D0
           if(j.le.2) then
              symz = -1.D0
           else
              symz = 1.D0
           end if
        else
           symy = 1.D0
           symz = 1.D0
        end if
        !backwards definition
        hyy(j,i) = gyy_f((i-1)*4+j)
        hyz(j,i) = symy*symz*gyz_f((i-1)*4+j)
        hzz(j,i) = gzz_f((i-1)*4+j)
        if (fisheye_enable==1) then
           Rp(j,i) = PhysR_f((i-1)*4+j)
           end if
     end do
  end do
!  r = Xglobal((pointcoords(1,1)-Xglobal(1))/dX+1)
  r = pointcoords(1,1)
  !==================================================================
  ! Finally, compute h_+ and h_x
  !================================================================== 
  hzz_ll = -0.0625D0*(hzz(1,1)+hzz(1,4)) + 0.5625D0*(hzz(1,2)+hzz(1,3))
  hzz_l = -0.0625D0*(hzz(2,1)+hzz(2,4)) + 0.5625D0*(hzz(2,2)+hzz(2,3))
  hzz_u = -0.0625D0*(hzz(3,1)+hzz(3,4)) + 0.5625D0*(hzz(3,2)+hzz(3,3))
  hzz_uu = -0.0625D0*(hzz(4,1)+hzz(4,4)) + 0.5625D0*(hzz(4,2)+hzz(4,3))
  hyz_ll = -0.0625D0*(hyz(1,1)+hyz(1,4)) + 0.5625D0*(hyz(1,2)+hyz(1,3))
  hyz_l = -0.0625D0*(hyz(2,1)+hyz(2,4)) + 0.5625D0*(hyz(2,2)+hyz(2,3))
  hyz_u = -0.0625D0*(hyz(3,1)+hyz(3,4)) + 0.5625D0*(hyz(3,2)+hyz(3,3))
  hyz_uu = -0.0625D0*(hyz(4,1)+hyz(4,4)) + 0.5625D0*(hyz(4,2)+hyz(4,3))
  hyy_ll = -0.0625D0*(hyy(1,1)+hyy(1,4)) + 0.5625D0*(hyy(1,2)+hyy(1,3))
  hyy_l = -0.0625D0*(hyy(2,1)+hyy(2,4)) + 0.5625D0*(hyy(2,2)+hyy(2,3))
  hyy_u = -0.0625D0*(hyy(3,1)+hyy(3,4)) + 0.5625D0*(hyy(3,2)+hyy(3,3))
  hyy_uu = -0.0625D0*(hyy(4,1)+hyy(4,4)) + 0.5625D0*(hyy(4,2)+hyy(4,3))
  hzz_cent = -0.0625D0*(hzz_ll+hzz_uu) + 0.5625D0*(hzz_l+hzz_u)
  hyz_cent = -0.0625D0*(hyz_ll+hyz_uu) + 0.5625D0*(hyz_l+hyz_u)
  hyy_cent = -0.0625D0*(hyy_ll+hyy_uu) + 0.5625D0*(hyy_l+hyy_u)
  if (fisheye_enable==1) then
     r_ll = -0.0625D0*(Rp(1,1)+Rp(1,4)) + 0.5625D0*(Rp(1,2)+Rp(1,3))
     r_l = -0.0625D0*(Rp(2,1)+Rp(2,4)) + 0.5625D0*(Rp(2,2)+Rp(2,3))
     r_u = -0.0625D0*(Rp(3,1)+Rp(3,4)) + 0.5625D0*(Rp(3,2)+Rp(3,3))
     r_uu = -0.0625D0*(Rp(4,1)+Rp(4,4)) + 0.5625D0*(Rp(4,2)+Rp(4,3))
     r = -0.0625D0*(r_ll+r_uu) + 0.5625D0*(r_l+r_u)
  end if
  amplitudes(1) = r * (hzz_cent - hyy_cent)/2.0
  amplitudes(2) = r * hyz_cent
end subroutine x_axis_waves
!-------------------------------------------------------------------------------
!
! Simplified version of interpolate_pointset: simply returns set of function 
!   values on stencil.
!
!------------------------------------------------------------------------------
subroutine return_funcvals_at_points(cctkGH,nghostzones,numpoints,ext,global_ext,dX,dY,dZ,Xlocal1d,Ylocal1d,Zlocal1d,Xglobal,Ygloba&
  &l,Zglobal,pointcoords,func,outputvals)
  implicit none
  REAL*8  excision_radius
  REAL*8  run_time
  INTEGER*4 Symmetry
  INTEGER*4 bssn_enable
  INTEGER*4 cowling_enable
  INTEGER*4 excision_enable
  INTEGER*4 fisheye_enable
  INTEGER*4 iter_count
  INTEGER*4 number_of_mol_ministeps
  INTEGER*4 rot_metric
  INTEGER*4 trA_detg_enforce
  COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count&
  &,number_of_mol_ministeps,rot_metric,trA_detg_enforce
  REAL*8  dR_GW
  REAL*8  ddR_GW
  REAL*8  phi_GW
  REAL*8  radius_GW
  REAL*8  radius_GW_phys
  REAL*8  theta_GW
  COMMON /gw_extractionrest/dR_GW,ddR_GW,phi_GW,radius_GW,radius_GW_phys,theta_GW
  REAL*8  compute_Psi4_max_radius
  REAL*8  compute_Psi4_min_radius
  REAL*8  radius_GW_Psi4(101)
  integer*8  psif_vec
  INTEGER*4 enable_interp_onepoint
  INTEGER*4 nmodes_Psi4
  INTEGER*4 nmodes_ZM
  INTEGER*4 num_extraction_radii
  INTEGER*4 radius_power
  INTEGER*4 scale_with_radius
  INTEGER*4 use_Rij_from_compute_ricci
  COMMON /gw_extractionpriv/compute_Psi4_max_radius,compute_Psi4_min_radius,radius_GW_Psi4,psif_vec,enable_interp_onepoint,nmodes_P&
  &si4,nmodes_ZM,num_extraction_radii,radius_power,scale_with_radius,use_Rij_from_compute_ricci
  REAL*8  CCTKH0
  REAL*8  CCTKH6
  REAL*8  CCTKH7
  REAL*8  CCTKH8
  REAL*8  CCTKH9
  REAL*8  CCTKH13
  REAL*8  rbr
  REAL*8  CCTKH17
  REAL*8  surf_radius
  INTEGER*4 CCTKH1
  INTEGER*4 CCTKH2
  INTEGER*4 CCTKH3
  INTEGER*4 CCTKH4
  INTEGER*4 CCTKH5
  INTEGER*4 CCTKH10
  INTEGER*4 CCTKH11
  INTEGER*4 CCTKH12
  INTEGER*4 CCTKH14
  INTEGER*4 CCTKH15
  INTEGER*4 CCTKH16
  INTEGER*4 numphi
  INTEGER*4 numtheta
  INTEGER*4 CCTKH18
  INTEGER*4 CCTKH19
  COMMON /DIAGNOSTICS_VACUUMrest/CCTKH0,CCTKH6,CCTKH7,CCTKH8,CCTKH9,CCTKH13,rbr,CCTKH17,surf_radius,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCT&
  &KH5,CCTKH10,CCTKH11,CCTKH12,CCTKH14,CCTKH15,CCTKH16,numphi,numtheta,CCTKH18,CCTKH19
  REAL*8  CCTKH25
  REAL*8  CCTKH39
  REAL*8  CCTKH47
  REAL*8  CCTKH49
  REAL*8  CCTKH51
  REAL*8  CCTKH53
  REAL*8  CCTKH55
  REAL*8  CCTKH57
  REAL*8  CCTKH59
  REAL*8  CCTKH61
  REAL*8  CCTKH63
  integer*8  CCTKH22
  integer*8  CCTKH23
  integer*8  CCTKH26
  integer*8  CCTKH29
  integer*8  CCTKH30
  integer*8  CCTKH31
  integer*8  CCTKH34
  integer*8  CCTKH35
  integer*8  CCTKH40
  integer*8  CCTKH41
  integer*8  CCTKH43
  integer*8  CCTKH65
  integer*8  CCTKH67
  integer*8  CCTKH69
  integer*8  CCTKH71
  integer*8  CCTKH72
  integer*8  CCTKH75
  INTEGER*4 CCTKH20
  INTEGER*4 CCTKH21
  INTEGER*4 CCTKH24
  INTEGER*4 CCTKH27
  INTEGER*4 CCTKH28
  INTEGER*4 CCTKH32
  INTEGER*4 CCTKH33
  INTEGER*4 CCTKH36
  INTEGER*4 CCTKH37
  INTEGER*4 CCTKH38
  INTEGER*4 out_every
  INTEGER*4 CCTKH42
  INTEGER*4 CCTKH44
  INTEGER*4 CCTKH45
  INTEGER*4 CCTKH46
  INTEGER*4 CCTKH48
  INTEGER*4 CCTKH50
  INTEGER*4 CCTKH52
  INTEGER*4 CCTKH54
  INTEGER*4 CCTKH56
  INTEGER*4 CCTKH58
  INTEGER*4 CCTKH60
  INTEGER*4 CCTKH62
  INTEGER*4 CCTKH64
  INTEGER*4 CCTKH66
  INTEGER*4 CCTKH68
  INTEGER*4 CCTKH70
  INTEGER*4 CCTKH73
  INTEGER*4 CCTKH74
  COMMON /IOrest/CCTKH25,CCTKH39,CCTKH47,CCTKH49,CCTKH51,CCTKH53,CCTKH55,CCTKH57,CCTKH59,CCTKH61,CCTKH63,CCTKH22,CCTKH23,CCTKH26,CC&
  &TKH29,CCTKH30,CCTKH31,CCTKH34,CCTKH35,CCTKH40,CCTKH41,CCTKH43,CCTKH65,CCTKH67,CCTKH69,CCTKH71,CCTKH72,CCTKH75,CCTKH20,CCTKH21,CC&
  &TKH24,CCTKH27,CCTKH28,CCTKH32,CCTKH33,CCTKH36,CCTKH37,CCTKH38,out_every,CCTKH42,CCTKH44,CCTKH45,CCTKH46,CCTKH48,CCTKH50,CCTKH52,&
  &CCTKH54,CCTKH56,CCTKH58,CCTKH60,CCTKH62,CCTKH64,CCTKH66,CCTKH68,CCTKH70,CCTKH73,CCTKH74
  
  external     CCTK_PointerTo
  integer*8 CCTK_PointerTo
  interface
  integer function CCTK_Equals (arg1, arg2)
  implicit none
  integer*8 arg1
  character(*) arg2
  end function CCTK_Equals
  integer function CCTK_MyProc (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_MyProc
  integer function CCTK_nProcs (cctkGH)
  implicit none
  integer*8 cctkGH
  end function CCTK_nProcs
  integer function CCTK_IsThornActive (name)
  implicit none
  character(*) name
  end function CCTK_IsThornActive
  integer*8 function CCTK_NullPointer ()
  implicit none
  end function CCTK_NullPointer
  end interface
  
  integer*8 :: cctkgh
  integer, dimension(3)                    :: ext,nghostzones,global_ext
  integer                                  :: numpoints
  real*8, dimension(3)                     :: offset
  real*8, dimension(numpoints,3)           :: pointcoords
  real*8, dimension(3)                     :: localpointcoord
  real*8,dimension(ext(1),ext(2),ext(3))   :: func
  real*8,dimension(ext(1))                 :: Xlocal1d
  real*8,dimension(ext(2))                 :: Ylocal1d
  real*8,dimension(ext(3))                 :: Zlocal1d
  real*8,dimension(global_ext(1))          :: Xglobal
  real*8,dimension(global_ext(2))          :: Yglobal
  real*8,dimension(global_ext(3))          :: Zglobal
  real*8                                   :: dX,dY,dZ,teeny,oneodX,oneodY,oneodZ
  real*8                                   :: interpvaluej,interpvaluek
  real*8,dimension(numpoints)              :: outputvals,funcvalues
  integer                                  :: n,ii,jj,kk
  integer, dimension(1)                    :: ivalue,jvalue,kvalue
  integer, dimension(1)                    :: ivalueglobal,jvalueglobal,kvalueglobal
  integer                                  :: imin,jmin,kmin,imax,jmax,kmax
  integer                                  :: handle,index,ierr
  integer                                  :: nsymghostzonesx,nsymghostzonesy,nsymghostzonesz
  integer                                  :: NO_SYMM, EQUATORIAL, OCTANT, PI_SYMM, AXISYM
  parameter(NO_SYMM = 0, EQUATORIAL = 1, OCTANT = 2, PI_SYMM = 3, AXISYM = 4)
  imin = 0
  jmin = 0
  kmin = 0
  imax = ext(1)
  jmax = ext(2)
  kmax = ext(3)
  ivalue = 1
  jvalue = 1
  kvalue = 1
  ivalueglobal = 1
  jvalueglobal = 1
  kvalueglobal = 1
  teeny = dX*1.D-5
  oneodx = 1.D0/dX
  oneody = 1.D0/dY
  oneodz = 1.D0/dZ
  nsymghostzonesx = 0
  nsymghostzonesy = 0
  nsymghostzonesz = 0
  if(Symmetry.eq.EQUATORIAL) then
!     nsymghostzonesz=1
     nsymghostzonesz=nghostzones(3)
  else if(Symmetry.eq.AXISYM) then
     nsymghostzonesx=nghostzones(1)-1
     nsymghostzonesz=nghostzones(3)
  else if(Symmetry.eq.OCTANT) then
     nsymghostzonesx=nghostzones(1)
     nsymghostzonesy=nghostzones(2)
     nsymghostzonesz=nghostzones(3)
  end if
  !By default, assume upper bound is processor boundary and not outer boundary.
  !We subtract # of ghostzones*2 from ext(i) to avoid overlap in our stencil points
  imax = ext(1)-nghostzones(1)*2
  jmax = ext(2)-nghostzones(2)*2
  kmax = ext(3)-nghostzones(3)*2
  !check if upper proc bound is global boundary and fix ijkmax
  if(abs(Xlocal1d(ext(1))-Xglobal(global_ext(1))).lt.teeny) imax = ext(1)
  if(abs(Ylocal1d(ext(2))-Yglobal(global_ext(2))).lt.teeny) jmax = ext(2)
  if(abs(Zlocal1d(ext(3))-Zglobal(global_ext(3))).lt.teeny) kmax = ext(3)
  do n=1,numpoints
     !We define ijkvalueglobal(1) = corner of stencil closest to (-inf,-inf,-inf)
     ivalueglobal(1) = (pointcoords(n,1)-Xglobal(1))*oneodx+1
     jvalueglobal(1) = (pointcoords(n,2)-Yglobal(1))*oneody+1
     kvalueglobal(1) = (pointcoords(n,3)-Zglobal(1))*oneodz+1
     !If point touches outside of global domain, move it back into global domain
     if(ivalueglobal(1).gt.global_ext(1)) then
        ivalueglobal(1) = global_ext(1)
        write(*,*) "OUTSIDE X BOUND"
     end if
     if(jvalueglobal(1).gt.global_ext(2)) then
        jvalueglobal(1) = global_ext(2)
        write(*,*) "OUTSIDE Y BOUND"
     end if
     if(kvalueglobal(1).gt.global_ext(3)) then
        kvalueglobal(1) = global_ext(3)
        write(*,*) "OUTSIDE Z BOUND"
     end if
     !Now that we know global point, hunt for that point on local processor!
     localpointcoord(1) = Xglobal(ivalueglobal(1))+teeny
     localpointcoord(2) = Yglobal(jvalueglobal(1))+teeny
     localpointcoord(3) = Zglobal(kvalueglobal(1))+teeny
     !Next search for coordinate on local domain!
     ivalue(1) = (localpointcoord(1)-Xlocal1d(1))*oneodx+1
     jvalue(1) = (localpointcoord(2)-Ylocal1d(1))*oneody+1
     kvalue(1) = (localpointcoord(3)-Zlocal1d(1))*oneodz+1
     !If the glove doesnt fit, you must acquit!
     !I.e., when search returns <0 or >ijkmax, that means the point is off the local processor grid!
     if(ivalue(1).lt.1 .or. ivalue(1).gt.imax) ivalue(1) = -1
     if(jvalue(1).lt.1 .or. jvalue(1).gt.jmax) jvalue(1) = -1
     if(kvalue(1).lt.1 .or. kvalue(1).gt.kmax) kvalue(1) = -1
!     write(*,*) hix.,ivalue(1),n,pointcoords(n,1),Xglobal(ivalueglobal(1))
!     write(*,*) hiy.,jvalue(1),n,pointcoords(n,2),Yglobal(jvalueglobal(1))
!     write(*,*) hiz.,kvalue(1),n,pointcoords(n,3),Zglobal(kvalueglobal(1))
     funcvalues(n) = 1.D307
!     write(*,*) hif,ivalue(1),jvalue(1),kvalue(1)
     if(ivalue(1).gt.0 .and. jvalue(1).gt.0 .and. kvalue(1).gt.0) then
!        write(*,*) hif2,ivalue(1),jvalue(1),kvalue(1)
        ii = ivalue(1)
        jj = jvalue(1)
        kk = kvalue(1)
        funcvalues(n) = func(ii,jj,kk)
     end if
  end do
  outputvals = 0.D0
  call CCTK_ReductionHandle(handle,"minimum")
  call CCTK_ReduceLocArrayToArray1D(ierr,cctkGH,-1,handle,&
       funcvalues,outputvals,numpoints,107)
end subroutine return_funcvals_at_points
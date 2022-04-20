/*@@
  @file      AHFinder_output.F
  @date      Fri 2 Aug 2001
  @author    Thomas Radke
  @desc 
             Does the output of AHFinder data into files.
  @enddesc 
@@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

      subroutine AHFinder_InitOutput(CCTK_ARGUMENTS)
      use AHFinder_dat
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      
      character(len=3) ::   ext
      character(len=200) :: almf,areaf,massf,radf,circeqf,merip1f,merip2f
      character(len=200) :: asymxf,asymyf,asymzf
      character(len=200) :: disf
      character(len=200) :: gaussf, gauss0f, gauss1f, gauss2f

      myproc = CCTK_MyProc(cctkGH)
      if (myproc.ne.0) then
        return
      end if

!     ****************************
!     ***   BUILD FILE NAMES   ***
!     ****************************

      call CCTK_FortranString(nfile,out_dir,filestr)

c     some compilers cannot trim an empty string, so we add at least one char
      if (find3) then
         ext = "_" // char(ichar('0') + mfind - 1) // "."
      else
         ext = "."
      end if

      almf = filestr(1:nfile) // "/ahf_coeff" // trim(ext) // "alm"
      areaf = filestr(1:nfile) // "/ahf_area" // trim(ext) // "tl"
      massf = filestr(1:nfile) // "/ahf_mass" // trim(ext) // "tl"
      radf  = filestr(1:nfile) // "/ahf_rad" // trim(ext) // "tl"

      circeqf = filestr(1:nfile) // "/ahf_circ_eq" // trim(ext) // "tl"
      merip1f = filestr(1:nfile) // "/ahf_meri_p1" // trim(ext) // "tl"
      merip2f = filestr(1:nfile) // "/ahf_meri_p2" // trim(ext) // "tl"

      asymxf = filestr(1:nfile) // "/ahf_asymx" // trim(ext) // "tl"
      asymyf = filestr(1:nfile) // "/ahf_asymy" // trim(ext) // "tl"
      asymzf = filestr(1:nfile) // "/ahf_asymz" // trim(ext) // "tl"

      disf = filestr(1:nfile) // "/ahf_d12" // trim(ext) // "tl"

      gaussf = filestr(1:nfile) // "/ahf.gauss"
      gauss0f = filestr(1:nfile) // "/ahf_0.gauss"
      gauss1f = filestr(1:nfile) // "/ahf_1.gauss"
      gauss2f = filestr(1:nfile) // "/ahf_2.gauss"


      open(11,file=massf,form='formatted',status='replace')
      write(11,"(A14)") '" Horizon mass'
      close(11)

      open(11,file=radf,form='formatted',status='replace')
      write(11,"(A26)") '" Horizon radius'
      close(11)

      open(11,file=circeqf,form='formatted',status='replace')
      write(11,"(A26)") '" Equatorial circumference'
      close(11)

      open(11,file=merip1f,form='formatted',status='replace')
      write(11,"(A27)") '" Length of meridian, phi=0'
      close(11)

      open(11,file=merip2f,form='formatted',status='replace')
      write(11,"(A30)") '" Length of meridian, phi=pi/2'
      close(11)

      open(11,file=areaf,form='formatted',status='replace')
      write(11,"(A14)") '" Horizon area'
      close(11)

      open(11,file=asymxf,form='formatted',status='replace')
      write(11,"(A29)") '" Asymmetries on x reflection'
      close(11)

      open(11,file=asymyf,form='formatted',status='replace')
      write(11,"(A29)") '" Asymmetries on y reflection'
      close(11)

      open(11,file=asymzf,form='formatted',status='replace')
      write(11,"(A29)") '" Asymmetries on z reflection'
      close(11)

      open(11,file=almf,form='formatted',status='replace')
      write(11,"(A21)")    '# Radial coefficients'
      write(11,"(A1)")     '#'
      close(11)

      open(11,file=disf,form='formatted',status='replace')
      write(11,"(A18)") '" Horizon distance'
      close(11)

      if (find3) then
        call AHFinder_WriteGaussHeader(CCTK_PASS_FTOF, gauss0f)
        call AHFinder_WriteGaussHeader(CCTK_PASS_FTOF, gauss1f)
        call AHFinder_WriteGaussHeader(CCTK_PASS_FTOF, gauss2f)
      else
        call AHFinder_WriteGaussHeader(CCTK_PASS_FTOF, gaussf)
      end if

      end subroutine AHFinder_InitOutput



      
      subroutine AHFinder_WriteGaussHeader(CCTK_ARGUMENTS, filename)
      use AHFinder_dat
      implicit none
      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      character(len=200) :: filename

      call AHFinder_SetReflections(CCTK_ARGUMENTS)

      open(11,file=filename,form='formatted', status='replace')
      
      write(11,"(A20)") '# GAUSSIAN CURVATURE'
      write(11,"(A1)")  '#'
      write(11,"(A35)") '# The data is written in a loop as:'
      write(11,"(A1)")  '#'
      write(11,"(A17)") '# do i=1,ntheta'
      write(11,"(A18)") '#    do j=1,nphi'
      write(11,"(A27)") '#       write gaussian(i,j)'
      write(11,"(A11)") '#    end do'
      write(11,"(A8)")  '# end do'
      write(11,"(A1)")  '#'
      write(11,"(A40)") '# theta and phi are subdivided uniformly'
      write(11,"(A26)") '# according to symmetries:'
      write(11,"(A1)")  '#'
      write(11,"(A36)") '# phi=[0,2 pi]   (refx=refy=.false.)'
      write(11,"(A44)") '# phi=[0,pi]     (refx=.false., refy=.true.)'
      write(11,"(A35)") '# phi=[0,pi/2]   (refx=refy=.true.)'
      write(11,"(A1)")  '#'
      write(11,"(A31)") '# theta=[0,pi]   (refz=.false.)'
      write(11,"(A30)") '# theta=[0,pi/2] (refz=.true.)'
      write(11,"(A1)")  '#'

      write(11,"(A9,L1)")  '# refx = ',refx
      write(11,"(A9,L1)")  '# refy = ',refy
      write(11,"(A9,L1)")  '# refz = ',refz
      write(11,"(A1)")     '#'
      write(11,"(A11,I7)") '# ntheta = ',(ahf_ntheta+1)
      write(11,"(A11,I7)") '# nphi   = ',(ahf_nphi+1)

      close(11)

      end subroutine AHFinder_WriteGaussHeader




      subroutine AHFinder_Output(CCTK_ARGUMENTS, report, status, horizon, mtype, intarea_h)

      use AHFinder_dat

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer mtype
      logical report, status, horizon
      logical massf_exists, radf_exists, circeqf_exists
      logical merip1f_exists, merip2f_exists, areaf_exists
      logical asymxf_exists, asymyf_exists, asymzf_exists, almf_exists
      CCTK_REAL intarea_h

      integer l, m, ierror
      character(len=3) ::   ext
      character(len=200) :: almf,areaf,massf,radf,circeqf,merip1f,merip2f
      character(len=200) :: asymxf,asymyf,asymzf
      character(len=200) :: options


!     ****************************
!     ***   BUILD FILE NAMES   ***
!     ****************************

c     some compilers cannot trim an empty string, so we add at least one char
      if (find3) then
         ext = "_" // char(ichar('0') + mfind - 1) // "."
      else
         ext = "."
      end if

      almf = filestr(1:nfile) // "/ahf_coeff" // trim(ext) // "alm"
      areaf = filestr(1:nfile) // "/ahf_area" // trim(ext) // "tl"
      massf = filestr(1:nfile) // "/ahf_mass" // trim(ext) // "tl"
      radf  = filestr(1:nfile) // "/ahf_rad" // trim(ext) // "tl"

      circeqf = filestr(1:nfile) // "/ahf_circ_eq" // trim(ext) // "tl"
      merip1f = filestr(1:nfile) // "/ahf_meri_p1" // trim(ext) // "tl"
      merip2f = filestr(1:nfile) // "/ahf_meri_p2" // trim(ext) // "tl"

      asymxf = filestr(1:nfile) // "/ahf_asymx" // trim(ext) // "tl"
      asymyf = filestr(1:nfile) // "/ahf_asymy" // trim(ext) // "tl"
      asymzf = filestr(1:nfile) // "/ahf_asymz" // trim(ext) // "tl"


!     *********************************
!     ***   OUTPUT GRID FUNCTIONS   ***
!     *********************************

      if (ahf_2Doutput.ne.0) then
         if (.not.find3) then
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahfgrid","IOFlexIO_2D")
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahf_exp","IOFlexIO_2D")
         else if (mfind.eq.3) then
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahfgrid3","IOFlexIO_2D")
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahf_exp3","IOFlexIO_2D")
         end if
         if (ierror < 0) then
           call CCTK_WARN(2, "Failed to output 2D data! For this to work, you must have IOFlexIO enabled")
         end if
      end if

      if (ahf_3Doutput.ne.0) then
         if (.not.find3) then
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahfgrid","IOFlexIO_3D")
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahf_exp","IOFlexIO_3D")
         else if (mfind.eq.3) then
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahfgrid","IOFlexIO_3D")
            call CCTK_OutputVarByMethod(ierror,cctkGH,
     .           "ahfinder::ahf_exp","IOFlexIO_3D")
         end if
         if (ierror < 0) then
           call CCTK_WARN(2, "Failed to output 3D data! For this to work, you must have IOFlexIO enabled")
         end if
      end if


!     **************************************************
!     ***   OUTPUT ASCII DATA ONLY ON PROCESSOR 0    ***
!     **************************************************

      if (myproc.eq.0) then

!        Mass.
        
        inquire(file=massf, exist=massf_exists)
        if (massf_exists) then
          open(11,file=massf,form='formatted',status='old',position='append')
        else
          open(11,file=massf,form='formatted',status='replace')
          write(11,"(A14)") '" Horizon mass'
        end if

        if (status.and.report) then
           out_mass = 0.141047396D0 * sqrt(intarea_h)
        else
           out_mass = 0.0
        end if
        write(11,"(2ES14.6)") cctk_time, out_mass

        close(11)

!        Radius.

        inquire(file=radf, exist=radf_exists)
        if (radf_exists) then
          open(11,file=radf,form='formatted',status='old',position='append')
        else
          open(11,file=radf,form='formatted',status='replace')
          write(11,"(A26)") '" Horizon radius'
        end if

        if (status.and.report) then
           out_radius = c0(0)
        else
           out_radius = 0.0
        end if
        write(11,"(2ES14.6)") cctk_time, out_radius

        close(11)

!        Equatorial circumference.

        inquire(file=circeqf, exist=circeqf_exists)
        if (circeqf_exists) then
          open(11,file=circeqf,form='formatted',status='old',position='append')
        else
          open(11,file=circeqf,form='formatted',status='replace')
          write(11,"(A26)") '" Equatorial circumference'
        end if

        if (status.and.report) then
           out_perimeter = circ_eq
        else
           out_perimeter = 0.0
        end if
        write(11,"(2ES14.6)") cctk_time, out_perimeter

        close(11)

!        Meridians.

        inquire(file=merip1f, exist=merip1f_exists)
        if (merip1f_exists) then
          open(11,file=merip1f,form='formatted',status='old',position='append')
        else
          open(11,file=merip1f,form='formatted',status='replace')
          write(11,"(A27)") '" Length of meridian, phi=0'
        end if

        if (status.and.report) then
           out_meridian1 = meri_p1
        else
           out_meridian1 = 0.0
        end if
        write(11,"(2ES14.6)") cctk_time, out_meridian1

        close(11)

        inquire(file=merip2f, exist=merip2f_exists)
        if (merip2f_exists) then
          open(11,file=merip2f,form='formatted',status='old', position='append')
        else
          open(11,file=merip2f,form='formatted',status='replace')
          write(11,"(A30)") '" Length of meridian, phi=pi/2'
        end if

        if (status.and.report) then
           out_meridian2 = meri_p2
        else
           out_meridian2 = 0.0
        end if
        write(11,"(2ES14.6)") cctk_time, out_meridian2

        close(11)

!        Area.
        
        inquire(file=areaf, exist=areaf_exists)
        if (areaf_exists) then
          open(11,file=areaf,form='formatted',status='old',position='append')
        else
          open(11,file=areaf,form='formatted',status='replace')
          write(11,"(A14)") '" Horizon area'
        end if

        if (status.and.report) then
           out_area = intarea_h
        else
           out_area = 0.0
        end if
        write(11,"(2ES14.6)") cctk_time, out_area

        close(11)

!        Asymmetries in x.  Coefficients cc with odd m, and
!        coefficients cs with even m indicate asymmetries
!        on refelction x<->-x.

        inquire(file=asymxf, exist=asymxf_exists)
        if (asymxf_exists) then
          open(11,file=asymxf,form='formatted',status='old', position='append')
        else
          open(11,file=asymxf,form='formatted',status='replace')
          write(11,"(A29)") '" Asymmetries on x reflection'
        end if

        out_asymx = 0.0D0
        if (status.and.report.and..not.refx) then
           do l=1,lmax
              do m=1,l
                 if (mod(m,2).ne.0) then
                    out_asymx = out_asymx + abs(cc(l,m))
                 else
                    out_asymx = out_asymx + abs(cs(l,m))
                 end if
              end do
           end do
           out_asymx = out_asymx / c0(0)
        end if
        write(11,"(2ES14.6)") cctk_time, out_asymx

        close(11)

!        Asymmetries in y.  Any cs coefficient different from
!        zero indicates asymmetries on reflection y<->-y.

        inquire(file=asymyf, exist=asymyf_exists)
        if (asymyf_exists) then
          open(11,file=asymyf,form='formatted',status='old', position='append')
        else
          open(11,file=asymyf,form='formatted',status='replace')
          write(11,"(A29)") '" Asymmetries on y reflection'
        end if

        out_asymy = 0.0D0
        if (status.and.report.and..not.refy) then
           do l=1,lmax
              do m=1,l
                 out_asymy = out_asymy + abs(cs(l,m))
              end do
           end do
           out_asymy = out_asymy / c0(0)
        end if
        write(11,"(2ES14.6)") cctk_time, out_asymy

        close(11)

!        Asymmetries in z. Any coefficient c0 with odd l,
!        or any coefficient cc and cs with odd (l-m) indicate
!        asymmetries on reflection z<->-z.

        inquire(file=asymzf, exist=asymzf_exists)
        if (asymzf_exists) then
          open(11,file=asymzf,form='formatted',status='old',position='append')
        else
          open(11,file=asymzf,form='formatted',status='replace')
          write(11,"(A29)") '" Asymmetries on z reflection'
        end if

        out_asymz = 0.0D0
        if (status.and.report.and..not.refz) then
           do l=1,lmax
              if (mod(l,2).ne.0) out_asymz = out_asymz + abs(c0(l))
              do m=1,l
                 if (mod(l-m,2).ne.0) then
                    out_asymz = out_asymz + abs(cc(l,m)) + abs(cs(l,m))
                 end if
              end do
           end do
           out_asymz = out_asymz / c0(0)
        end if
        write(11,"(2ES14.6)") cctk_time, out_asymz

        close(11)

!       Expansion coefficients.

        inquire(file=almf, exist=almf_exists)
        if (almf_exists) then
          open(11,file=almf,form='formatted',status='old',position='append')
        else
          open(11,file=almf,form='formatted',status='replace')
          write(11,"(A21)")    '# Radial coefficients'
          write(11,"(A1)")     '#'
        end if

        write(11,"(A12,I4)") '# Time step ',cctk_iteration
        write(11,"(A7,ES14.6)") '# Time ',cctk_time
        write(11,"(A7,I4)")  '# Call ',ahf_ncall

        if (status.and.report) then
           if (mtype.eq.1) then
              if (horizon) then
                 write(11,"(A30,I4)") '# Surface found: Outer horizon'
              else
                 write(11,"(A31,I4)") '# Surface found: Outer horizon?'
              end if
           else if (mtype.eq.2) then
              if (horizon) then
                 write(11,"(A30,I4)") '# Surface found: Inner horizon'
              else
                 write(11,"(A31,I4)") '# Surface found: Inner horizon?'
              end if
           else if (mtype.eq.3) then
              if (horizon) then
                 write(11,"(A34,I4)") '# Surface found: Marginal horizon?'
              else
                 write(11,"(A32,I4)") '# Surface found: Trapped surface'
              end if
           else if (mtype.eq.4) then
              if (horizon) then
                 write(11,"(A34,I4)") '# Surface found: Marginal horizon?'
              else
                 write(11,"(A30,I4)") '# Surface found: Not a horizon'
              end if
           end if

           out_centerx = xc
           out_centery = yc
           out_centerz = zc

           write(11,"(A15,3ES14.6)") '# centered on: ', out_centerx, out_centery, out_centerz
  
           write(11,"(A1)")  '#'
           write(11,"(A22)") '#     a_lm       l   m'
           write(11,"(A1)")  '#'

           write(11,"(ES14.6,2I4)") c0(0),0,0

           if (lmax .gt. 0) then
              out_c0 = c0

              out_cs = cs
              out_cc = cc

              do l=1,lmax
                 do m=l,1,-1
                    write(11,"(ES14.6,2I4)") out_cs(l,m),l,-m
                 end do
                 write(11,"(ES14.6,2I4)") out_c0(l),l,0
                 do m=1,l
                    write(11,"(ES14.6,2I4)") out_cc(l,m),l,m
                 end do
              end do
           end if
          write(11,*)
        else
           write(11,"(A18,I4)") '# No surface found'
        end if

        close(11)

      end if

!     ***********************************************
!     ***   OUTPUT AHFINDER DATA IN HDF5 FORMAT   ***
!     ***********************************************

      options = "mfind=" // char(48 + mfind) // " mtype=" // char(48 + mtype)
      if (status) then
         options = options(1:15) // " status=1"
      else
         options = options(1:15) // " status=0"
      end if
      options = options // " "
      if (horizon) then
         options = options(1:24) // " horizon=1"
      else
         options = options(1:24) // " horizon=0"
      end if

      if (ahf_HDF5output.ne.0) then
        call CCTK_OutputVarByMethod (ierror, cctkGH, options, "IOAHFinderHDF5")
        if (ierror .ne. 0) then
           call CCTK_WARN (1, 'Error calling I/O method "IOAHFinderHDF5" (I/O method not activated ?)')
        end if
      end if


      end subroutine AHFinder_Output

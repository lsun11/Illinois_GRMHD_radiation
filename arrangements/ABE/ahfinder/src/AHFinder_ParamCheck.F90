!  @@
!  @file      ParamCheck.F90
!  @date      Mon Apr 29 17:51:25 2002
!  @author    Peter Diener
!  @desc 
!  Parameter checking stuff for AHFinder
!  @enddesc
!  @version $Header: /cactus/CactusEinstein/AHFinder/src/AHFinder_ParamCheck.F90,v 1.2 2003/10/27 15:31:27 schnetter Exp $
! @@

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! CCTK_FILEVERSION(CactusEinstein_AHFinder_ParamCheck_f)

!  @@
!  @routine    AHFinder_ParamCheck
!  @date       Thu Apr 29 17:51:25 2002
!  @author     Peter Diener
!  @desc 
!  Scheduled routine to detect invalid parameter settings.
!  @enddesc 
!  @calls     
!  @calledby   
!  @history 
!
!  @endhistory 
!
! @@

subroutine AHFinder_ParamCheck(CCTK_ARGUMENTS)
 
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  if ( ( .not. CCTK_EQUALS(metric_type,'physical') ) .and. &
       ( .not. CCTK_EQUALS(metric_type,'static conformal') ) ) then
    call CCTK_PARAMWARN('Unknown ADMBase::metric_type - known types are "physical" and "static conformal"')
  end if
end subroutine AHFinder_ParamCheck

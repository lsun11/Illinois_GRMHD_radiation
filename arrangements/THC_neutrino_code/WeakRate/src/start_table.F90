#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

integer function StartTableReader()

  implicit none

  DECLARE_CCTK_PARAMETERS

  INTEGER :: fortran_length, weak_table_reader, err

  character*200 fortran_filename

  call CCTK_FortranString(fortran_length,table_filename,fortran_filename)
  write(*,*) 'start_table01'
  err = weak_table_reader(fortran_filename)

  StartTableReader = err

end function StartTableReader

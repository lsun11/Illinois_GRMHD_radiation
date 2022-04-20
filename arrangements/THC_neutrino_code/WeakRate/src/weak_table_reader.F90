#include <weaktableindex.h>
! #########################################################
! TABLE UNITS:
!        density              g/cm^3
!        temperature          MeV
!        ye                   number fraction per baryon
! #########################################################

INTEGER function weak_table_reader(myFilename)
  ! This routine reads the table and initializes
  ! all variables in the module.

  use hdf5
  use table3d_mod
  use lk_hdfOpenField_mod

  implicit none
  character*200,  INTENT(IN) :: myFilename
  character*200 ::  myFieldName
  integer :: error

  write(6,*) "Reading table with weak interaction rates: ", myFilename

  !DENSITY
  myFieldName = "density"
  error = lk_hdfOpenField(myFilename,myFieldName,logrho)
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, rho"
    weak_table_reader = -1
    return
  endif
  logrho = log10(logrho)

  !TEMPERATURE
  myFieldName = "temperature"
  error = lk_hdfOpenField(myFilename,myFieldName,logtemp)
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, temp"
    weak_table_reader = -1
    return
  endif
  logtemp = log10(logtemp)

  !Ye
  myFieldName = "ye"
  error = lk_hdfOpenField(myFilename,myFieldName,yeTable)
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, ye"
    weak_table_reader = -1
    return
  endif

  !Mass factor in MeV
  myFieldName = "mass_factor"
  error = lk_hdfOpenField(myFilename,myFieldName,mass_fact)
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, mass_factor"
    weak_table_reader = -1
    return
  endif
  !mass_fact = 9.223158894119980d+02

  nrho = size(logrho,1)
  ntemp = size(logtemp,1)
  nye = size(yeTable,1)
  write(6,*) 'Number of points in the table:: ', nrho, ntemp, nye

  !Allocate the number of variables needed to compute the weak rates
  nvars = 8

  allocate(allvariables(nrho,ntemp,nye,nvars))

  myFieldName = "mu_n"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,MU_N))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, MU_N"
    weak_table_reader = -1
    return
  endif

  myFieldName = "mu_p"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,MU_P))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, MU_P"
    weak_table_reader = -1
    return
  endif
  write(6,*) "mu_p",  allvariables(214,136,20,MU_P)

  myFieldName = "mu_e"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,MU_E))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, MU_E"
    weak_table_reader = -1
    return
  endif
  write(6,*) "mu_e", allvariables(214,136,20,MU_E)

  myFieldName = "Abar"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,ABAR))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, ABAR"
    weak_table_reader = -1
    return
  endif
  write(6,*) "Abar", allvariables(214,136,20,ABAR)

  myFieldName = "Zbar"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,ZBAR))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, ZBAR"
    weak_table_reader = -1
    return
  endif
  write(6,*) "Zbar", allvariables(214,136,20,ZBAR)

  myFieldName = "Xp"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,XP))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, XP"
    weak_table_reader = -1
    return
  endif
  write(6,*) "Xp", allvariables(214,136,20,XP)

  myFieldName = "Xn"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,XN))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, XN"
    weak_table_reader = -1
    return
  endif
  write(6,*) "Xn", allvariables(214,136,20,XN)

  myFieldName = "Xh"
  error = lk_hdfOpenField(myFilename,myFieldName,allvariables(:,:,:,XH))
  if (error.ne.0) then
    write(6,*) "WeakRates :: problem reading the table, XH"
    weak_table_reader = -1
    return
  endif
  write(6,*) "Xh", allvariables(214,136,20,XH)

  ! set min-max values:
  eos_lrhomin = logrho(1)
  eos_lrhomax = logrho(nrho)
  eos_rhomin = 10.**logrho(1)
  eos_rhomax = 10.**logrho(nrho)

  eos_ltempmin = logtemp(1)
  eos_ltempmax = logtemp(ntemp)
  eos_tempmin = 10.**eos_ltempmin
  eos_tempmax = 10.**eos_ltempmax

  eos_yemin = yeTable(1)
  eos_yemax = yeTable(nye)

  ! Set the spacings
  dlrho = logrho(nrho)-logrho(nrho-1)
  dye = yeTable(nye) - yeTable(nye-1)
  dltemp = logtemp(ntemp) - logtemp(ntemp-1)

  write(6,*) "WeakRates: Done reading table"
  write(6,*)

  weak_table_reader = 0
end function weak_table_reader

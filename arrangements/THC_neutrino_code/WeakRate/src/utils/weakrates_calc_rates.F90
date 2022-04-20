#include <table3d_mod.F90>
#include <weak_constants.F90>

program weakrates_calc_rates
    ! -----------------------------------------------------------------
    ! Variable declarations
    ! -----------------------------------------------------------------
    use table3d_mod
    use weak_constants

    implicit none
    character*200 table_file_name
    integer ierr
    integer weak_table_reader
    integer Emissions_cgs
    integer Opacities_cgs
    integer NeutrinoDens_cgs

    real*8 rho, temp, ye
    real*8 R_nue, R_nua, R_nux
    real*8 Q_nue, Q_nua, Q_nux
    real*8 kappa_0_nue, kappa_0_nua, kappa_0_nux
    real*8 kappa_1_nue, kappa_1_nua, kappa_1_nux
    real*8 n_nue, n_nua, n_nux
    real*8 e_nue, e_nua, e_nux

    real*8 r_cgs2cactus, q_cgs2cactus, edens_cgs2cactus
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Read input file
    ! -----------------------------------------------------------------
    open(111, file='weakrates_calc_rates.inp', status='unknown', &
            form='formatted', action='read')
    read(111, 200) table_file_name
    read(111, 201) rho
    read(111, 201) temp
    read(111, 201) ye
    close(111)

    write(*,*) 'Table file name : ', trim(table_file_name)
    write(*,*) 'rho  [g/cm^3]   = ', rho
    write(*,*) 'temp [MeV]      = ', temp
    write(*,*) 'Y_e             = ', ye

200 format(a200)
201 format(e16.9)
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Read table and compute values
    ! -----------------------------------------------------------------
    ierr = weak_table_reader(table_file_name)
    if(ierr.ne.0) then
        write(6,*) 'Failed reading ', trim(table_file_name)
        stop
    end if

    ierr = Emissions_cgs(rho, temp, ye, &
        R_nue, R_nua, R_nux, Q_nue, Q_nua, Q_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the emission rates'
        stop
    end if

    ierr = Opacities_cgs(rho, temp, ye, &
            kappa_0_nue, kappa_0_nua, kappa_0_nux, &
            kappa_1_nue, kappa_1_nua, kappa_1_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the opacities'
        stop
    end if

    ierr = NeutrinoDens_cgs(rho, temp, ye, &
           n_nue, n_nua, n_nux, e_nue, e_nua, e_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the neutrino densities'
        stop
    end if
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Output results CGS
    ! -----------------------------------------------------------------
    write(*,*) '------------ CGS -------------------------------------'

    write(*,*) 'Number emission rates [1/(cm^3 s)]'
    write(*,*) '        nue     = ', R_nue
    write(*,*) '        nua     = ', R_nua
    write(*,*) '        nux     = ', R_nux

    write(*,*) 'Energy emission rates [erg/(cm^3 s)]'
    write(*,*) '        nue     = ', Q_nue * mev_to_erg
    write(*,*) '        nua     = ', Q_nua * mev_to_erg
    write(*,*) '        nux     = ', Q_nux * mev_to_erg

    write(*,*) 'Number opacities [1/cm]'
    write(*,*) '        nue     = ', kappa_0_nue
    write(*,*) '        nua     = ', kappa_0_nua
    write(*,*) '        nux     = ', kappa_0_nux

    write(*,*) 'Energy opacities [1/cm]'
    write(*,*) '        nue     = ', kappa_1_nue
    write(*,*) '        nua     = ', kappa_1_nua
    write(*,*) '        nux     = ', kappa_1_nux

    write(*,*) 'Neutrino number density [1/cm^3]'
    write(*,*) '        nue     = ', n_nue
    write(*,*) '        nua     = ', n_nua
    write(*,*) '        nux     = ', n_nux

    write(*,*) 'Neutrino energy density [erg/cm^3]'
    write(*,*) '        nue     = ', e_nue * mev_to_erg
    write(*,*) '        nua     = ', e_nua * mev_to_erg
    write(*,*) '        nux     = ', e_nux * mev_to_erg

    write(*,*) '------------ Cactus ----------------------------------'

    write(*,*) 'Number emission rates'
    r_cgs2cactus = 1. / (cgs2cactusTime*cgs2cactusLength**3)
    write(*,*) '        nue     = ', R_nue * r_cgs2cactus
    write(*,*) '        nua     = ', R_nua * r_cgs2cactus
    write(*,*) '        nux     = ', R_nux * r_cgs2cactus

    write(*,*) 'Energy emission rates'
    q_cgs2cactus = mev_to_erg*cgs2cactusenergy/ &
                   (cgs2cactusTime * cgs2cactusLength**3)
    write(*,*) '        nue     = ', Q_nue * q_cgs2cactus
    write(*,*) '        nua     = ', Q_nue * q_cgs2cactus
    write(*,*) '        nux     = ', Q_nue * q_cgs2cactus

    write(*,*) 'Number opacities'
    write(*,*) '        nue     = ', kappa_0_nue / cgs2cactusLength
    write(*,*) '        nua     = ', kappa_0_nua / cgs2cactusLength
    write(*,*) '        nux     = ', kappa_0_nux / cgs2cactusLength

    write(*,*) 'Energy opacities'
    write(*,*) '        nue     = ', kappa_1_nue / cgs2cactusLength
    write(*,*) '        nua     = ', kappa_1_nua / cgs2cactusLength
    write(*,*) '        nux     = ', kappa_1_nux / cgs2cactusLength

    write(*,*) 'Neutrino number density'
    write(*,*) '        nue     = ', n_nue / (cgs2cactusLength**3)
    write(*,*) '        nua     = ', n_nue / (cgs2cactusLength**3)
    write(*,*) '        nux     = ', n_nue / (cgs2cactusLength**3)

    write(*,*) 'Neutrino energy density'
    edens_cgs2cactus = mev_to_erg * cgs2cactusEnergy / cgs2cactusLength**3
    write(*,*) '        nue     = ', e_nue * edens_cgs2cactus
    write(*,*) '        nua     = ', e_nua * edens_cgs2cactus
    write(*,*) '        nux     = ', e_nux * edens_cgs2cactus
    ! -----------------------------------------------------------------
end program

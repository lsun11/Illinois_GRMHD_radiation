program weakrates_test_rates
    ! -----------------------------------------------------------------
    ! Variable declarations
    ! -----------------------------------------------------------------
    implicit none
    character*200 table_file_name
    integer ierr
    integer weak_table_reader
    integer Emissions_cgs
    integer Opacities_cgs
    integer NeutrinoDens_cgs

    ! Error and tolerance
    real*8 err, eps

    ! Reference point
    real*8 ref_rho, ref_temp, ref_ye

    ! Reference value
    real*8 ref_Rbeta_nue, ref_Rbeta_nua, ref_Rbeta_nux
    real*8 ref_Qbeta_nue, ref_Qbeta_nua, ref_Qbeta_nux
    real*8 ref_Rpair_nue, ref_Rpair_nua, ref_Rpair_nux
    real*8 ref_Qpair_nue, ref_Qpair_nua, ref_Qpair_nux
    real*8 ref_Rplsm_nue, ref_Rplsm_nua, ref_Rplsm_nux
    real*8 ref_Qplsm_nue, ref_Qplsm_nua, ref_Qplsm_nux
    real*8 ref_R_nue, ref_R_nua, ref_R_nux
    real*8 ref_Q_nue, ref_Q_nua, ref_Q_nux

    real*8 ref_kappa_0_nue, ref_kappa_0_nua, ref_kappa_0_nux
    real*8 ref_kappa_1_nue, ref_kappa_1_nua, ref_kappa_1_nux

    real*8 ref_n_nue, ref_n_nua, ref_n_nux
    real*8 ref_e_nue, ref_e_nua, ref_e_nux

    ! Values computed by WeakRates
    real*8 my_R_nue, my_R_nua, my_R_nux
    real*8 my_Q_nue, my_Q_nua, my_Q_nux

    real*8 my_kappa_0_nue, my_kappa_0_nua, my_kappa_0_nux
    real*8 my_kappa_1_nue, my_kappa_1_nua, my_kappa_1_nux

    real*8 my_n_nue, my_n_nua, my_n_nux
    real*8 my_e_nue, my_e_nua, my_e_nux
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Reference values
    ! -----------------------------------------------------------------
    ! Reference point
    ref_rho  = 7.582629253507815d+14 ! g/cm^3
    ref_temp = 16.852590447507527  ! MeV
    ref_ye   = 0.234693877551020

    ! Number emissivities (1/(sec cm^3))
    ref_Rbeta_nue = 1.379324165823913d+44
    ref_Rbeta_nua = 5.671639161779945d+36
    ref_Rbeta_nux = 0.0e0

    ref_Rpair_nue = 6.680945831638925d+35
    ref_Rpair_nua = 6.680945831638925d+35
    ref_Rpair_nux = 6.680945831638925d+35

    ref_Rplsm_nue = 5.980430160212835d+36
    ref_Rplsm_nua = 5.980430160212835d+36
    ref_Rplsm_nux = 1.125673694269156d+37

    ref_R_nue = ref_Rbeta_nue + ref_Rpair_nue + ref_Rplsm_nue
    ref_R_nua = ref_Rbeta_nua + ref_Rpair_nua + ref_Rplsm_nua
    ref_R_nux = ref_Rbeta_nux + ref_Rpair_nux + ref_Rplsm_nux

    ! Energy emissivities (MeV/(sec cm^3))
    ref_Qbeta_nue = 3.490444173353645d+46
    ref_Qbeta_nua = 4.779090609608149d+38
    ref_Qbeta_nux = 0.0e0

    ref_Qpair_nue = 1.028101444972817d+38
    ref_Qpair_nua = 1.028101444972817d+38
    ref_Qpair_nux = 1.030850507254617d+38

    ref_Qplsm_nue = 1.243103898095211d+38
    ref_Qplsm_nua = 1.243103898095211d+38
    ref_Qplsm_nux = 2.339847335127855d+38

    ref_Q_nue = ref_Qbeta_nue + ref_Qpair_nue + ref_Qplsm_nue
    ref_Q_nua = ref_Qbeta_nua + ref_Qpair_nua + ref_Qplsm_nua
    ref_Q_nux = ref_Qbeta_nux + ref_Qpair_nux + ref_Qplsm_nux

    ! Opacities (1/cm)
    ref_kappa_0_nue = 0.132606720395860
    ref_kappa_0_nua = 0.034458169701843
    ref_kappa_0_nux = 0.037144174652465
    ref_kappa_1_nue = 0.161306100056073
    ref_kappa_1_nua = 0.057426366767866
    ref_kappa_1_nux = 0.059745262109698

    ! Neutrino densities (1/cm^3)
    ref_n_nue = 4.930208545514009d+36
    ref_n_nua = 4.077159500359407d+31
    ref_n_nux = 2.276099545195539d+35

    ! Neutrino energy density (erg/cm^3)
    ref_e_nue = 8.563762504824248e+32
    ref_e_nua = 3.302885005308371d+27
    ref_e_nux = 1.936764547826890d+31
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Compute values
    ! -----------------------------------------------------------------
    ! Initialize Table
    table_file_name = "WeakRates.hdf5"
    ierr = weak_table_reader(table_file_name)
    if(ierr.ne.0) then
        write(6,*) 'Failed reading WeakRates.hdf5'
        stop
    end if

    ierr = Emissions_cgs(ref_rho, ref_temp, ref_ye, &
        my_R_nue, my_R_nua, my_R_nux, my_Q_nue, my_Q_nua, my_Q_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the emission rates'
        stop
    end if

    ierr = Opacities_cgs(ref_rho, ref_temp, ref_ye, &
            my_kappa_0_nue, my_kappa_0_nua, my_kappa_0_nux, &
            my_kappa_1_nue, my_kappa_1_nua, my_kappa_1_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the opacities'
        stop
    end if

    ierr = NeutrinoDens_cgs(ref_rho, ref_temp, ref_ye, &
           my_n_nue, my_n_nua, my_n_nux,  &
           my_e_nue, my_e_nua, my_e_nux)
    if(ierr.ne.0) then
        write(6,*) 'Failed computing the neutrino densities'
        stop
    end if
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! Check values
    ! -----------------------------------------------------------------
    eps = 1.0d-7

    err = abs(my_R_nue - ref_R_nue)/ref_R_nue
    if(err.gt.eps) then
        write(6,*) 'R_nue   =', my_R_nue
        write(6,*) 'ref     =', ref_R_nue
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_R_nua - ref_R_nua)/ref_R_nua
    if(err.gt.eps) then
        write(6,*) 'R_nua   =', my_R_nue
        write(6,*) 'ref     =', ref_R_nua
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_R_nux - ref_R_nux)/ref_R_nux
    if(err.gt.eps) then
        write(6,*) 'R_nux   =', my_R_nux
        write(6,*) 'ref     =', ref_R_nux
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_Q_nue - ref_Q_nue)/ref_Q_nue
    if(err.gt.eps) then
        write(6,*) 'Q_nue   =', my_Q_nue
        write(6,*) 'ref     =', ref_Q_nue
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_Q_nua - ref_Q_nua)/ref_Q_nua
    if(err.gt.eps) then
        write(6,*) 'Q_nua   =', my_Q_nue
        write(6,*) 'ref     =', ref_Q_nua
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_Q_nux - ref_Q_nux)/ref_Q_nux
    if(err.gt.eps) then
        write(6,*) 'Q_nux   =', my_Q_nux
        write(6,*) 'ref     =', ref_Q_nux
        write(6,*) 'rel.err =', err
        write(6,*)
    end if

    err = abs(my_kappa_0_nue - ref_kappa_0_nue)/ref_kappa_0_nue
    if(err.gt.eps) then
        write(6,*) 'kappa_0_nue =', my_kappa_0_nue
        write(6,*) 'ref         =', ref_kappa_0_nue
        write(6,*) 'rel.err     =', err
        write(6,*)
    end if
    err = abs(my_kappa_0_nua - ref_kappa_0_nua)/ref_kappa_0_nua
    if(err.gt.eps) then
        write(6,*) 'kappa_0_nua =', my_kappa_0_nua
        write(6,*) 'ref         =', ref_kappa_0_nua
        write(6,*) 'rel.err     =', err
        write(6,*)
    end if
    err = abs(my_kappa_0_nux - ref_kappa_0_nux)/ref_kappa_0_nux
    if(err.gt.eps) then
        write(6,*) 'kappa_0_nux =', my_kappa_0_nux
        write(6,*) 'ref         =', ref_kappa_0_nux
        write(6,*) 'rel.err     =', err
        write(6,*)
    end if
    err = abs(my_kappa_1_nue - ref_kappa_1_nue)/ref_kappa_1_nue
    if(err.gt.eps) then
        write(6,*) 'kappa_1_nue =', my_kappa_1_nue
        write(6,*) 'ref         =', ref_kappa_1_nue
        write(6,*) 'rel.err     =', err
        write(6,*)
    end if
    err = abs(my_kappa_1_nua - ref_kappa_1_nua)/ref_kappa_1_nua
    if(err.gt.eps) then
        write(6,*) 'kappa_1_nua =', my_kappa_1_nua
        write(6,*) 'ref         =', ref_kappa_1_nua
        write(6,*) 'rel.err     =', err
        write(6,*)
    end if
    err = abs(my_kappa_1_nux - ref_kappa_1_nux)/ref_kappa_1_nux
    if(err.gt.eps) then
        write(6,*) 'kappa_1_nux =', my_kappa_1_nux
        write(6,*) 'ref         =', ref_kappa_1_nux
        write(6,*) 'rel.err     =', err
        write(6,*)
    end if

    err = abs(my_n_nue - ref_n_nue)/ref_n_nue
    if(err.gt.eps) then
        write(6,*) 'n_nue   =', my_n_nue
        write(6,*) 'ref     =', ref_n_nue
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_n_nua - ref_n_nua)/ref_n_nua
    if(err.gt.eps) then
        write(6,*) 'n_nua   =', my_n_nua
        write(6,*) 'ref     =', ref_n_nua
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_n_nux - ref_n_nux)/ref_n_nux
    if(err.gt.eps) then
        write(6,*) 'n_nux   =', my_n_nux
        write(6,*) 'ref     =', ref_n_nux
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_e_nue - ref_e_nue)/ref_e_nue
    if(err.gt.eps) then
        write(6,*) 'e_nue   =', my_e_nue
        write(6,*) 'ref     =', ref_e_nue
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_e_nua - ref_e_nua)/ref_e_nua
    if(err.gt.eps) then
        write(6,*) 'e_nua   =', my_e_nua
        write(6,*) 'ref     =', ref_e_nua
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    err = abs(my_e_nux - ref_e_nux)/ref_e_nux
    if(err.gt.eps) then
        write(6,*) 'e_nux   =', my_e_nux
        write(6,*) 'ref     =', ref_e_nux
        write(6,*) 'rel.err =', err
        write(6,*)
    end if
    ! -----------------------------------------------------------------
end program

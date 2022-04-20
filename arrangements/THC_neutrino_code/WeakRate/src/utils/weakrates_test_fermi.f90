program weakrates_test_fermi
    implicit none
    integer i
    real*8 err, eps
    real*8 fermi0, fermi1, fermi2, fermi3, fermi4, fermi5
    real*8 fermi4o2, fermi4o3, fermi5o3, fermi5o4
    real*8, dimension(4) :: rfermi0, rfermi1, rfermi2, rfermi3, rfermi4, rfermi5
    real*8, dimension(4) :: rfermi4o2, rfermi4o3, rfermi5o3, rfermi5o4
    real*8, dimension(4) :: eta

    eps = 1.0d-7

    eta       = (/ 1.0d-4, 1.0d0, 100.0d0, 1000.0d0 /)
    ! Reference values (from an independent implementation)
    rfermi0   = (/ 0.6931971818099453d0, 1.313261687518223d0,  &
                   100d0,                1000d0               /)
    rfermi1   = (/ 0.8225053673325695d0, 1.809505094932662d0,  &
                   5001.644899999999d0,  500001.6449d0        /)
    rfermi2   = (/ 1.803265838339036d0,  4.319966849146549d0,  &
                   333662.3233333333d0,  333336623.2333333d0  /)
    rfermi3   = (/ 5.682897261903042d0,  14.40330317089867d0,  &
                   25049359.3644d0,      250004934811.3644d0  /)
    rfermi4   = (/ 23.33268992689878d0,  60.91826854366128d0,  &
                   2006584245.76d0,      200006579745457.6d0  /)
    rfermi5   = (/ 118.2732202851887d0,  314.7867623129019d0,  &
                   167490273342.1989d0,  1.666748914803108d17 /)
    rfermi4o2 = (/ 12.93912934567108d0,  14.10155926444117d0,  &
                   6013.817280039122d0,  600013.8172800005d0  /)
    rfermi4o3 = (/ 4.105773666421222d0,  4.229465131772296d0,  &
                   80.10521213615331d0,  800.0105273776618d0  /)
    rfermi5o3 = (/ 20.81213416932023d0,  21.85517853633162d0,  &
                   6686.409456851624d0,  666686.4060346312d0  /)
    rfermi5o4 = (/ 5.068992073170225d0,  5.167362268139454d0,  &
                   83.47034204824101d0,  833.3470413445044d0  /)

    do i = 1, 4
        err = abs(fermi0(eta(i)) - rfermi0(i))/rfermi0(i)
        if(err.gt.eps) then
            write(6,*) 'fermi0: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi1(eta(i)) - rfermi1(i))/rfermi1(i)
        if(err.gt.eps) then
            write(6,*) 'fermi1: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi2(eta(i)) - rfermi2(i))/rfermi2(i)
        if(err.gt.eps) then
            write(6,*) 'fermi2: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi3(eta(i)) - rfermi3(i))/rfermi3(i)
        if(err.gt.eps) then
            write(6,*) 'fermi3: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi4(eta(i)) - rfermi4(i))/rfermi4(i)
        if(err.gt.eps) then
            write(6,*) 'fermi4: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi5(eta(i)) - rfermi5(i))/rfermi5(i)
        if(err.gt.eps) then
            write(6,*) 'fermi5: eta = ', eta(i), 'rel.err = ', err
        end if
    end do

    do i = 1, 4
        err = abs(fermi4o2(eta(i)) - rfermi4o2(i))/rfermi4o2(i)
        if(err.gt.eps) then
            write(6,*) 'fermi4o2: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi4o3(eta(i)) - rfermi4o3(i))/rfermi4o3(i)
        if(err.gt.eps) then
            write(6,*) 'fermi4o3: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi5o3(eta(i)) - rfermi5o3(i))/rfermi5o3(i)
        if(err.gt.eps) then
            write(6,*) 'fermi5o3: eta = ', eta(i), 'rel.err = ', err
        end if

        err = abs(fermi5o4(eta(i)) - rfermi5o4(i))/rfermi5o4(i)
        if(err.gt.eps) then
            write(6,*) 'fermi5o4: eta = ', eta(i), 'rel.err = ', err
        end if
    end do
end program

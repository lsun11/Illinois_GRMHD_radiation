module weak_constants

  implicit none

  !UNITS CONVERSION
  real(kind=8) :: mev_to_erg = 1.60217733d-6
  real(kind=8) :: erg_to_mev = 6.24150636d5
  real(kind=8) :: amu_cgs = 1.66053873d-24 !Atomic mass unit in g
  real(kind=8) :: amu_mev = 931.49432d0 !Atomic mass unit in MeV
  real(kind=8) :: kb_erg = 1.380658d-16 !Boltzmann constant in erg
  real(kind=8) :: kb_mev = 8.61738568d-11 !Boltzmann constant in MeV
  real(kind=8) :: temp_mev_to_kelvin = 1.1604447522806d10
  real(kind=8) :: fermi_cubed_cgs = 1.d-39
  real(kind=8) :: clite = 2.99792458d10 !speed of light in cm s^-1
  real(kind=8) :: me_mev = 0.510998910d0 !mass of the electron in MeV
  real(kind=8) :: me_erg = 8.187108692567103d-07 !mass of the electron in erg
  real(kind=8) :: sigma_0 = 1.76d-44 !cross section in unit of cm^2
  real(kind=8) :: alpha = 1.23d0 !dimensionless
  real(kind=8) :: multipl_nuea = 1 !multiplicity factor for nu_e and anti nu_e
  real(kind=8) :: multipl_nux = 2 !multiplicity factor for nu_x
  real(kind=8) :: Qnp = 1.293333d0 !neutron-proton mass difference in MeV
  real(kind=8) :: hc_mevcm = 1.23984172d-10 !hc in units of MeV*cm
  real(kind=8) :: hc_ergvcm = 1.986445683269303d-16 !hc in units of erg*cm/s
  real(kind=8) :: Cv = 0.5d0 + 2.0d0*0.23d0 !vector  const. dimensionless
  real(kind=8) :: Ca = 0.5d0 !axial const. dimensionless
  real(kind=8) :: gamma_0 = 5.565d-2 ! dimensionless
  real(kind=8) :: fsc = 1.0d0/137.036d0 !fine structure constant, dimensionless
  real(kind=8) :: planck = 6.626176d-27 !Planck constant erg s
  real(kind=8) :: avo = 6.0221367d23 !Avogadro's number mol^-1
  real(kind=8) :: Ggrav = 6.6742d-8
  real(kind=8) :: solarMass = 1.9891d33

end module weak_constants

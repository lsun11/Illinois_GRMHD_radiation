#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

CCTK_REAL FUNCTION AtomicMassImpl()

    use table3d_mod
    use weak_constants

    implicit none

    DECLARE_CCTK_PARAMETERS

    !mass_fact is in MeV
    AtomicMassImpl = normfact * cgs2cactusMass * mass_fact * mev_to_erg / (clite*clite)

END FUNCTION AtomicMassImpl


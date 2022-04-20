//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2012, David Radice <david.radice@aei.mpg.de>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void THC_LK_MoLRegister(CCTK_ARGUMENTS) {
    int ierr = 0;

    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::alp"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::betax"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::betay"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::betaz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gxx"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gxy"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gxz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gyy"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gyz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gzz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kxx"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kxy"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kxz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kyy"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kyz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kzz"));

    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::rho"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[0]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[1]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::vel[2]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::eps"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::press"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::w_lorentz"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::entropy"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::temperature"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("HydroBase::Y_e"));

    ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::zvec[0]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::zvec[1]"));
    ierr |= MoLRegisterConstrained(CCTK_VarIndex("THC_Core::zvec[2]"));

    if(ierr) {
        CCTK_ERROR("Could not register with MoL");
    }
}

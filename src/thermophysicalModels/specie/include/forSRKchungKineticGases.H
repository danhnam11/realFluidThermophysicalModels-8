/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef forSRKchungKineticGases_H
#define forSRKchungKineticGases_H

#include "rfSpecie.H"
#include "soaveRedlichKwong.H"
#include "rfJanafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "chungKineticTransport.H"
#include "thermo.H"
#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forSRKchungKineticGasEqns(Mu, He, Cp, Macro, Args...)                    \
    forThermo(Mu, He, Cp, soaveRedlichKwong, rfSpecie, Macro, Args)

#define forSRKchungKineticGasEnergiesAndThermos(Mu, Macro, Args...)               \
    forSRKchungKineticGasEqns(Mu, sensibleEnthalpy, rfJanafThermo, Macro, Args); \
    forSRKchungKineticGasEqns(Mu, sensibleInternalEnergy, rfJanafThermo, Macro, Args)

#define forSRKchungKineticGasTransports(Macro, Args...)                              \
    forSRKchungKineticGasEnergiesAndThermos(chungKineticTransport, Macro, Args)

#define forSRKchungKineticGases(Macro, Args...)                                   \
    forSRKchungKineticGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

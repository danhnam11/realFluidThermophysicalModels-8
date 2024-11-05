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

#include "ode.H"

//#include "kineticStandardChemistryModel.H"  //
#include "SRKchungTakaStandardChemistryModel.H"  //
#include "PRchungTakaStandardChemistryModel.H"  //
#include "SRKchungKineticStandardChemistryModel.H"  //
#include "PRchungKineticStandardChemistryModel.H"  //
#include "idKineticStandardChemistryModel.H"  //
#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

//#include "forKineticGases.H"  //
#include "forSRKchungTakaGases.H"  //
#include "forPRchungTakaGases.H"  //
#include "forSRKchungKineticGases.H"  //
#include "forPRchungKineticGases.H"  //
#include "forIdKineticGases.H"  //
#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeChemistrySolver.H"
#include "makeKineticChemistrySolver.H" //
#include "makeRealFluidChemistrySolver.H" //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //forKineticGases(makeKineticChemistrySolvers, ode, psiReactionThermo); //
    //forKineticGases(makeKineticChemistrySolvers, ode, rhoReactionThermo); //

    forSRKchungTakaGases(makeSRKchungTakaChemistrySolvers, ode, psiReactionThermo); //
    forSRKchungTakaGases(makeSRKchungTakaChemistrySolvers, ode, rhoReactionThermo); //

    forPRchungTakaGases(makePRchungTakaChemistrySolvers, ode, psiReactionThermo); //
    forPRchungTakaGases(makePRchungTakaChemistrySolvers, ode, rhoReactionThermo); //

    forSRKchungKineticGases(makeSRKchungKineticChemistrySolvers, ode, psiReactionThermo); //
    forSRKchungKineticGases(makeSRKchungKineticChemistrySolvers, ode, rhoReactionThermo); //
    forPRchungKineticGases(makePRchungKineticChemistrySolvers, ode, psiReactionThermo); //
    forPRchungKineticGases(makePRchungKineticChemistrySolvers, ode, rhoReactionThermo); //

    forIdKineticGases(makeIdKineticChemistrySolvers, ode, psiReactionThermo); //
    forIdKineticGases(makeIdKineticChemistrySolvers, ode, rhoReactionThermo); //


    forCommonGases(makeChemistrySolvers, ode, psiReactionThermo);
    forCommonGases(makeChemistrySolvers, ode, rhoReactionThermo);

    forCommonLiquids(makeChemistrySolvers, ode, rhoReactionThermo);

    forPolynomials(makeChemistrySolvers, ode, rhoReactionThermo);
}


// ************************************************************************* //

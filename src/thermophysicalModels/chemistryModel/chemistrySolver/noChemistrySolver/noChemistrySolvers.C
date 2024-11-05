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

#include "noChemistrySolver.H"

//#include "kineticStandardChemistryModel.H" //
#include "SRKchungTakaStandardChemistryModel.H" //
#include "PRchungTakaStandardChemistryModel.H" //
#include "SRKchungKineticStandardChemistryModel.H" //
#include "PRchungKineticStandardChemistryModel.H" //
#include "idKineticStandardChemistryModel.H" //
#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

//#include "forKineticGases.H" //
#include "forSRKchungTakaGases.H" //
#include "forPRchungTakaGases.H" //
#include "forSRKchungKineticGases.H" //
#include "forPRchungKineticGases.H" //
#include "forIdKineticGases.H" //
#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeKineticChemistrySolver.H" //
#include "makeRealFluidChemistrySolver.H" //
#include "makeChemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //forKineticGases(makeKineticChemistrySolvers, noChemistrySolver, psiReactionThermo); //
    //forKineticGases(makeKineticChemistrySolvers, noChemistrySolver, rhoReactionThermo); //

    forSRKchungTakaGases(makeSRKchungTakaChemistrySolvers, noChemistrySolver, psiReactionThermo); //
    forSRKchungTakaGases(makeSRKchungTakaChemistrySolvers, noChemistrySolver, rhoReactionThermo); //
    forPRchungTakaGases(makePRchungTakaChemistrySolvers, noChemistrySolver, psiReactionThermo); //
    forPRchungTakaGases(makePRchungTakaChemistrySolvers, noChemistrySolver, rhoReactionThermo); //

    forSRKchungKineticGases(makeSRKchungKineticChemistrySolvers, noChemistrySolver, psiReactionThermo); //
    forSRKchungKineticGases(makeSRKchungKineticChemistrySolvers, noChemistrySolver, rhoReactionThermo); //
    forPRchungKineticGases(makePRchungKineticChemistrySolvers, noChemistrySolver, psiReactionThermo); //
    forPRchungKineticGases(makePRchungKineticChemistrySolvers, noChemistrySolver, rhoReactionThermo); //

    forIdKineticGases(makeIdKineticChemistrySolvers, noChemistrySolver, psiReactionThermo); //
    forIdKineticGases(makeIdKineticChemistrySolvers, noChemistrySolver, rhoReactionThermo); //


    forCommonGases(makeChemistrySolvers, noChemistrySolver, psiReactionThermo);
    forCommonGases(makeChemistrySolvers, noChemistrySolver, rhoReactionThermo);

    forCommonLiquids
    (
        makeChemistrySolvers,
        noChemistrySolver,
        rhoReactionThermo
    );

    forPolynomials(makeChemistrySolvers, noChemistrySolver, rhoReactionThermo);
}


// ************************************************************************* //

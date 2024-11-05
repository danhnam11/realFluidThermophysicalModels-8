/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "infiniteReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "FallOffReactionRate.H"

#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "LangmuirHinshelwoodReactionRate.H"

#include "MichaelisMentenReactionRate.H"

#include "forCommonGases.H"
//#include "forKineticGases.H" //
#include "forSRKchungTakaGases.H" //
#include "forPRchungTakaGases.H" //
#include "forSRKchungKineticGases.H" //
#include "forPRchungKineticGases.H" //
#include "forIdKineticGases.H" //
#include "forCommonLiquids.H"
#include "forPolynomials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* const Foam::Tuple2<Foam::word, Foam::scalar>::typeName
(
    "Tuple2<word,scalar>"
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //forKineticGases(defineReaction, nullArg); //
    forSRKchungTakaGases(defineReaction, nullArg); //
    forPRchungTakaGases(defineReaction, nullArg); //
    forSRKchungKineticGases(defineReaction, nullArg); //
    forPRchungKineticGases(defineReaction, nullArg); //
    forIdKineticGases(defineReaction, nullArg); //
    forCommonGases(defineReaction, nullArg);
    forCommonLiquids(defineReaction, nullArg);
    forPolynomials(defineReaction, nullArg);


    // Irreversible/reversible/non-equilibrium-reversible reactions

    //forKineticGases(makeIRNReactions, ArrheniusReactionRate); //
    forSRKchungTakaGases(makeIRNReactions, ArrheniusReactionRate); //
    forPRchungTakaGases(makeIRNReactions, ArrheniusReactionRate); //
    forSRKchungKineticGases(makeIRNReactions, ArrheniusReactionRate); //
    forPRchungKineticGases(makeIRNReactions, ArrheniusReactionRate); //
    forIdKineticGases(makeIRNReactions, ArrheniusReactionRate); //
    forCommonGases(makeIRNReactions, ArrheniusReactionRate);
    forCommonLiquids(makeIRNReactions, ArrheniusReactionRate);
    forPolynomials(makeIRNReactions, ArrheniusReactionRate);

    //forKineticGases(makeIRNReactions, infiniteReactionRate); //
    forSRKchungTakaGases(makeIRNReactions, infiniteReactionRate); //
    forPRchungTakaGases(makeIRNReactions, infiniteReactionRate); //
    forSRKchungKineticGases(makeIRNReactions, infiniteReactionRate); //
    forPRchungKineticGases(makeIRNReactions, infiniteReactionRate); //
    forIdKineticGases(makeIRNReactions, infiniteReactionRate); //
    forCommonGases(makeIRNReactions, infiniteReactionRate);
    forCommonLiquids(makeIRNReactions, infiniteReactionRate);
    forPolynomials(makeIRNReactions, infiniteReactionRate);

    //forKineticGases(makeIRNReactions, LandauTellerReactionRate); //
    forSRKchungTakaGases(makeIRNReactions, LandauTellerReactionRate); //
    forPRchungTakaGases(makeIRNReactions, LandauTellerReactionRate); //
    forSRKchungKineticGases(makeIRNReactions, LandauTellerReactionRate); //
    forPRchungKineticGases(makeIRNReactions, LandauTellerReactionRate); //
    forIdKineticGases(makeIRNReactions, LandauTellerReactionRate); //
    forCommonGases(makeIRNReactions, LandauTellerReactionRate);
    forCommonLiquids(makeIRNReactions, LandauTellerReactionRate);
    forPolynomials(makeIRNReactions, LandauTellerReactionRate);

    //forKineticGases(makeIRNReactions, thirdBodyArrheniusReactionRate); //
    forSRKchungTakaGases(makeIRNReactions, thirdBodyArrheniusReactionRate); //
    forPRchungTakaGases(makeIRNReactions, thirdBodyArrheniusReactionRate); //
    forSRKchungKineticGases(makeIRNReactions, thirdBodyArrheniusReactionRate); //
    forPRchungKineticGases(makeIRNReactions, thirdBodyArrheniusReactionRate); //
    forIdKineticGases(makeIRNReactions, thirdBodyArrheniusReactionRate); //
    forCommonGases(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forCommonLiquids(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forPolynomials(makeIRNReactions, thirdBodyArrheniusReactionRate);


    // Irreversible/reversible reactions

    //forKineticGases(makeIRReactions, JanevReactionRate); //
    forSRKchungTakaGases(makeIRReactions, JanevReactionRate); //
    forPRchungTakaGases(makeIRReactions, JanevReactionRate); //
    forSRKchungKineticGases(makeIRReactions, JanevReactionRate); //
    forPRchungKineticGases(makeIRReactions, JanevReactionRate); //
    forIdKineticGases(makeIRReactions, JanevReactionRate); //
    forCommonGases(makeIRReactions, JanevReactionRate);
    forCommonLiquids(makeIRReactions, JanevReactionRate);
    forPolynomials(makeIRReactions, JanevReactionRate);

    //forKineticGases(makeIRReactions, powerSeriesReactionRate); //
    forSRKchungTakaGases(makeIRReactions, powerSeriesReactionRate); //
    forPRchungTakaGases(makeIRReactions, powerSeriesReactionRate); //
    forSRKchungKineticGases(makeIRReactions, powerSeriesReactionRate); //
    forPRchungKineticGases(makeIRReactions, powerSeriesReactionRate); //
    forIdKineticGases(makeIRReactions, powerSeriesReactionRate); //
    forCommonGases(makeIRReactions, powerSeriesReactionRate);
    forCommonLiquids(makeIRReactions, powerSeriesReactionRate);
    forPolynomials(makeIRReactions, powerSeriesReactionRate);


    // Pressure dependent reactions
    /*
    forKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    */
    forSRKchungTakaGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );//
    forPRchungTakaGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forSRKchungKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forPRchungKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forIdKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //


    forCommonGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    /*
    forKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    */
    forSRKchungTakaGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forPRchungTakaGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forSRKchungKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forPRchungKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forIdKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    /*
    forKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    */
    forSRKchungTakaGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forPRchungTakaGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forSRKchungKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forPRchungKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forIdKineticGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //


    forCommonGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );

    /*
    forKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    */
    forSRKchungTakaGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forPRchungTakaGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forSRKchungKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forPRchungKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //
    forIdKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    ); //


    forCommonGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    /*
    forKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    */
    forSRKchungTakaGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forPRchungTakaGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forSRKchungKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forPRchungKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //
    forIdKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    ); //


    forCommonGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    /*
    forKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    */
    forSRKchungTakaGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //

    forPRchungTakaGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forSRKchungKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forPRchungKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //
    forIdKineticGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    ); //


    forCommonGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
}

// ************************************************************************* //

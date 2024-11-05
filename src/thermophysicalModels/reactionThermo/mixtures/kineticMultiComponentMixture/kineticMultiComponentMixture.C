/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "kineticMultiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::PtrList<ThermoType>
Foam::kineticMultiComponentMixture<ThermoType>::readSpeciesData
(
    const dictionary& thermoDict
) const
{
    PtrList<ThermoType> specieThermos(species_.size());

    forAll(species_, i)
    {
        specieThermos.set
        (
            i,
            new ThermoType(thermoDict.subDict(species_[i]))
        );
    }

    return specieThermos;
}


template<class ThermoType>
typename Foam::kineticMultiComponentMixture<ThermoType>::speciesCompositionTable
Foam::kineticMultiComponentMixture<ThermoType>::readSpeciesComposition
(
    const dictionary& thermoDict,
    const speciesTable& species
) const
{
    speciesCompositionTable speciesComposition_;

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(species, si)
    {
        if (thermoDict.subDict(species[si]).isDict("elements"))
        {
            dictionary currentElements
            (
                thermoDict.subDict(species[si]).subDict("elements")
            );

            wordList currentElementsName(currentElements.toc());
            List<specieElement> currentComposition(currentElementsName.size());

            forAll(currentElementsName, eni)
            {
                currentComposition[eni].name() = currentElementsName[eni];

                currentComposition[eni].nAtoms() =
                    currentElements.lookupOrDefault
                    (
                        currentElementsName[eni],
                        0
                    );
            }

            // Add current specie composition to the hash table
            speciesCompositionTable::iterator specieCompositionIter
            (
                speciesComposition_.find(species[si])
            );

            if (specieCompositionIter != speciesComposition_.end())
            {
                speciesComposition_.erase(specieCompositionIter);
            }

            speciesComposition_.insert(species[si], currentComposition);
        }
    }

    return speciesComposition_;
}


template<class ThermoType>
void Foam::kineticMultiComponentMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(max(Yt).value()) < rootVSmall)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << exit(FatalError);
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::kineticMultiComponentMixture<ThermoType>::kineticMultiComponentMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        mesh,
        phaseName
    ),
    specieThermos_(readSpeciesData(thermoDict)),
    speciesComposition_(readSpeciesComposition(thermoDict, species())),
    mixture_("mixture", specieThermos_[0]),
    mixtureVol_("volMixture", specieThermos_[0]),
    numberOfSpecies_(species_.size()),  //
    // for diffusivity
    ListW_(species_.size()),
    linearityMk_(species_.size()),
    epsilonOverKbMk_(species_.size()),
    sigmaMk_(species_.size()),
    miuiMk_(species_.size()),
    polarMk_(species_.size()),
    ZrotMk_(species_.size()),
    CpCoeffTableMk_(species_.size()),

    EPSILONijOVERKB_(species_.size()),
    DELTAij_(species_.size()),
    Mij_(species_.size()),
    SIGMAij_(species_.size())
    //
{
    correctMassFractions();

    // precalculation for kinetic theory model
    //for Kinetic model
    forAll(ListW_, i)
    { 
        ListW_[i]           = specieThermos_[i].W();
        linearityMk_[i]     = specieThermos_[i].linearity();
        epsilonOverKbMk_[i] = specieThermos_[i].epsilonOverKb();
        sigmaMk_[i]         = specieThermos_[i].sigma();
        miuiMk_[i]          = specieThermos_[i].miui();
        polarMk_[i]         = specieThermos_[i].polar();
        ZrotMk_[i]          = specieThermos_[i].Zrot();        
        CpCoeffTableMk_[i]  = specieThermos_[i].CpCoeffTable();
    }

 
    List<scalar> nEPSILONijOVERKB(species_.size());
    List<scalar> nDELTAij(species_.size());
    List<scalar> nMij(species_.size());
    List<scalar> nSIGMAij(species_.size());    

    // Dip = dipole moment 
    const scalar DipMin = 1e-20;

    forAll(EPSILONijOVERKB_, i)
    {
        forAll(nEPSILONijOVERKB, j)
        {
            nMij[j] = 
               1/(1/specieThermos_[i].W() + 1/specieThermos_[j].W());

            nEPSILONijOVERKB[j] = 
               sqrt(specieThermos_[i].epsilonOverKb()*specieThermos_[j].epsilonOverKb());

            nSIGMAij[j] = 
               0.5*(specieThermos_[i].sigma() + specieThermos_[j].sigma());

            // calculate coeficient xi
            scalar xi = 1.0; 
            if ((specieThermos_[i].miui() < DipMin) && (specieThermos_[j].miui() > DipMin)) 
            {
                // miui_j > DipMin > miui_i --> j is polar, i is nonpolar              
                nDELTAij[j] = 0;
                xi = 1.0 + 
                     specieThermos_[i].polar()*pow(specieThermos_[j].miui(), 2)* 
                     sqrt(specieThermos_[j].epsilonOverKb()/specieThermos_[i].epsilonOverKb())/
                     (4*pow(specieThermos_[i].sigma(), 3)*
                            (
                             1e+19*specieThermos_[j].Kb()*specieThermos_[j].epsilonOverKb()*
                             pow(specieThermos_[j].sigma(), 3)
                            )
                     );
            }
            else if ((specieThermos_[i].miui() > DipMin) && (specieThermos_[j].miui() < DipMin))
            {
                // miui_j < DipMin < miui_i --> i is polar, j is nonpolar
                nDELTAij[j] = 0;
                xi = 1.0 +
                     specieThermos_[j].polar()*pow(specieThermos_[i].miui(), 2)*
                     sqrt(specieThermos_[i].epsilonOverKb()/specieThermos_[j].epsilonOverKb())/
                     (4*pow(specieThermos_[j].sigma(), 3)*   
                            (
                             1e+19*specieThermos_[i].Kb()*specieThermos_[i].epsilonOverKb()*
                             pow(specieThermos_[i].sigma(), 3)
                            )
                     ); 
            }
            else 
            {            
                xi = 1.0;
                nDELTAij[j] =
                   0.5*(specieThermos_[i].miui()*specieThermos_[j].miui())/
                   (nEPSILONijOVERKB[j]*1e+19*specieThermos_[j].Kb()*pow(nSIGMAij[j], 3));
            }
           
            nEPSILONijOVERKB[j] = nEPSILONijOVERKB[j]*pow(xi, 2);
            nSIGMAij[j] = nSIGMAij[j]*pow(xi, -1/6);
        }

        EPSILONijOVERKB_[i] = nEPSILONijOVERKB;
        Mij_[i]             = nMij;
        SIGMAij_[i]         = nSIGMAij;
        DELTAij_[i]         = nDELTAij;
    }
    // End of pre-calculation for kinetic model

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::kineticMultiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*specieThermos_[n];
    }

   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(Y_.size()); 
    List<scalar> Y(Y_.size()); 
    scalar sumXb = 0.0;  
    forAll(X, i)
    {
        sumXb = sumXb + Y_[i][celli]/ListW_[i]; 
    }  
    if (sumXb == 0){ sumXb = 1e-30;} 

    forAll(X, i)
    {
        X[i] = (Y_[i][celli]/ListW_[i])/sumXb;
        Y[i] = Y_[i][celli];
        if(X[i] <= 0) { X[i] = 0; }
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }
    
    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }

    
    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS
    (
        Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
        

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::kineticMultiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*specieThermos_[n];
    }

   //- Update coefficients for diffusivity mixture
    //- List of secies mole and mass fraction 
    List<scalar> X(Y_.size());
    List<scalar> Y(Y_.size());
    scalar sumXb = 0.0;
    forAll(X, i)
    {
        sumXb = sumXb + Y_[i].boundaryField()[patchi][facei]/ListW_[i];
    }
    if (sumXb == 0){ sumXb = 1e-30;}

    forAll(X, i)
    {
        X[i] = (Y_[i].boundaryField()[patchi][facei]/ListW_[i])/sumXb;
        Y[i] = Y_[i].boundaryField()[patchi][facei];
        if(X[i] <= 0) { X[i] = 0; } 
        if(Y[i] <= 0) { Y[i] = 0; }
    }

    scalar WmixCorrect = 0.0, sumXcorrected = 0.0;
    forAll(X, i)
    {
        X[i] = X[i] + 1e-40;
        sumXcorrected = sumXcorrected + X[i];
    }

    forAll(X, i)
    {
        X[i] = X[i]/sumXcorrected;
        WmixCorrect = WmixCorrect + X[i]*ListW_[i];
    }

    forAll(Y, i)
    {
        Y[i] = X[i]*ListW_[i]/WmixCorrect;
    }

   
    // Update coefficients for mixture of Kinetic model
    mixture_.updateTRANS
    (
        Y, X, EPSILONijOVERKB_, DELTAij_, Mij_, SIGMAij_,
        linearityMk_, epsilonOverKbMk_, sigmaMk_, miuiMk_, polarMk_, ZrotMk_, ListW_,
        CpCoeffTableMk_
    );
   

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::kineticMultiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(specieThermos_, i)
    {
        rhoInv += Y_[i][celli]/specieThermos_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/specieThermos_[0].rho(p, T)/rhoInv*specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/specieThermos_[n].rho(p, T)/rhoInv*specieThermos_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::kineticMultiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(specieThermos_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/specieThermos_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/specieThermos_[0].rho(p, T)/rhoInv
      * specieThermos_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/specieThermos_[n].rho(p,T)
          / rhoInv*specieThermos_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
void Foam::kineticMultiComponentMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        specieThermos_[i] = ThermoType(thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //

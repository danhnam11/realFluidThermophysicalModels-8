#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species() << exit(FatalIOError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// Nam

//- Creates Molecular weight field [kg/kmol] 
Info<< "Creating field molecular weight - W \n" << endl;
volScalarField Wmix
(
    IOobject
    (
        "Wmix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, 
        IOobject::NO_WRITE 
    ),
    thermo.W()
);

//- Creates Dynamic viscosity field [kg/m.s] 
Info<< "Creating field dynamic viscosity - mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

//- Creates Heat Capacity at constant pressure "Cp" [J/kgK] 
Info<< "Creating field heat capacity Cp\n" << endl;
volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);

//- Creates Heat Capacity at constant volume "Cv" [J/kgK] 
Info<< "Creating field heat capacity Cv\n" << endl;
volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cv()
);

//- Create Thermal Conductivity field"kappa" [W/mK] 
Info<< "Creating field thermal conductivity - kappa \n" << endl;
volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.kappa()
);

//- Create the total enthalpy field "he" [J/kg]-for generating properties only
Info<< "Creating field total enthaly - he\n" << endl;
volScalarField HE
(
    IOobject
    (
        "HE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),  
    thermo.he()
); 

//- Creates List of species name

   forAll(Y, i)
   {
       Info << " =========== print speciesName ==============  " << endl;
       Info << " ith = " << i << ";--> speciesName =  " << composition.species()[i] << endl;
   }

// correction terms for diffusion velocity
volScalarField sumDiffusionCorrections1
(
    IOobject
    (
        "gas_sumCorrections1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gas_sumCorrections1", dimensionSet(1, -3, -1, 0, 0), 0.)
);

volScalarField sumDiffusionCorrections2
(
    IOobject
    (
        "gas_sumCorrections2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gas_sumCorrections2", dimensionSet(1, -3, -1, 0, 0), 0.)
);

volScalarField sumDiffusionCorrections3
(
    IOobject
    (
        "gas_sumCorrections3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("gas_sumCorrections3", dimensionSet(1, -3, -1, 0, 0), 0.)
);
// --------

// - for Energy Equation
volScalarField sumHeatDiffusion
(
    IOobject
    (
        "sumHeatDiffusion",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sumHeatDiffusion", dimensionSet(1, -1, -3, 0, 0), 0.)
);

volScalarField sumHeatDiffusion2
(
    IOobject
    (
        "sumHeatDiffusion2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("sumHeatDiffusion2", dimensionSet(1, -1, -3, 0, 0), 0.)
);
// --------

PtrList<volVectorField> YVi(Y.size());
forAll(YVi, i) 
{  
    YVi.set
    (
        i, 
         new volVectorField 
         ( 
             IOobject
             (
               "YVi_"+composition.species()[i],
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::NO_WRITE 
             ), 
             mesh,
             dimensionedVector("YVi_"+composition.species()[i], dimensionSet(0,1,-1, 0, 0), vector(0.,0.,0.))
          )
    ); 
}


//- Creates List of mass diffusivity field Dimix_i[m^2/s] of species ith in the mixture
Info<< "Creating field mass diffusivity Dimix_i\n" << endl;
PtrList<volScalarField> Dimix(Y.size());
forAll(Dimix, i) 
{  
    Dimix.set
    (
        i, 
         new volScalarField 
         ( 
             IOobject
             (
               "Dimix_"+composition.species()[i],
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE 
             ), 
             thermo.Dimix(i)
         )
    ); 
}


//- Creates List of Enthalpy/Internal Energy field he_i[J/kg] of individual species
Info<< "\nCreating field individual enthalpy he_i \n" << endl;
PtrList<volScalarField> hei(Y.size());
forAll(hei, i) 
{  
    hei.set
    (
        i, 
         new volScalarField 
         ( 
             IOobject
             (
               "hei_"+composition.species()[i],
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE 
             ), 
             thermo.hei(i)
         )
    ); 
}

// 



Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& p = thermo.p();

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<psiReactionThermophysicalTransportModel> thermophysicalTransport
(
    psiReactionThermophysicalTransportModel::New(turbulence(), thermo)
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> reaction
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvOptions.H"

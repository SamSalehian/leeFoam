/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    leeFoam createFields

Description
	this file contains necessary field definitions for leeFoam solver

\*---------------------------------------------------------------------------*/

// creating/reading base flow fields

// mean pressure
Info << nl << "Reading base-flow pressure field" << endl;
volScalarField p0
(
 	IOobject
	(
	 	"p0",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
    mesh
);

// mean velocity vector
Info << nl << "Reading bese-flow velocity field" << endl;
volVectorField U0
(
 	IOobject
	(
	 	"U0",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
    mesh
);

// mean density 
Info << nl << "Reading base-flow density field" << endl;
volScalarField rho0
(
 	IOobject
	(
	 	"rho0",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
    mesh
);

// p' pressure pertibation
Info << nl << "Reading perturbation pressure field" << endl;
volScalarField p
(
 	IOobject
	(
	 	"p",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

// U' velocity pertibation vector
Info << nl << "Reading perturbation velocity field" << endl;
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

// rho' desnsity pertibation
Info << nl << "Reading perturbation density field" << endl;
volScalarField rho
(
 	IOobject
	(
	 	"rho",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
    mesh
);

// damping coef ; D = w*frequency*blendFactor
Info << nl << "Reading dampining coef field" << endl;
volScalarField D
(
 	IOobject
	(
	 	"D",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
    mesh
);

// Acoustic Relaxation Term (ART)
Info << nl << "Reading dampining coef field" << endl;
volScalarField D2
(
 	IOobject
	(
	 	"D2",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	),
    mesh
);

// creating a flux field
// this field is used in the matrix operations like laplace and divergence
// it is computed by interpolating velocity (volume field) to faces (surface field)
// and multiplying with face area.
Info << nl << "Creating flux field" << endl;
surfaceScalarField phi
(
 	IOobject
	(
	 	"phi",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::interpolate(U) & mesh.Sf()
);
Info << nl << "Creating flux field" << endl;
surfaceScalarField phi0
(
 	IOobject
	(
	 	"phi0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	fvc::interpolate(U0) & mesh.Sf()
);

Info << nl << "Creating flux field" << endl;
surfaceScalarField rhoPhi0_rho0
(
 	IOobject
	(
	 	"rhoPhi0_rho0",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	 linearInterpolate(rho*U0/rho0) & mesh.Sf()
);

// This part reads the constants defined in transportProperties
Info << nl << "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
// Read gamma from transportProperties
dimensionedScalar gamma("gamma", dimensionSet(0,0,0,0,0,0,0), transportProperties);
dimensionedScalar c0("c0", dimVelocity, transportProperties);

// declaring fvSolution dictionary entry
IOdictionary fvSolution
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
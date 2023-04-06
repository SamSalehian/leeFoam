/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Global
    LEECourantNo

Description
    Calculates and outputs the perturbation and base-flow Courant Numbers.

\*---------------------------------------------------------------------------*/

scalar Co0 = 0.0;
scalar meanCo0 = 0.0;
scalar Co = 0.0;
scalar meanCo = 0.0;

{
    scalarField sumPhi0
    (
        fvc::surfaceSum(mag(phi0))().primitiveField()
    );

    Co0 = 0.5*gMax(sumPhi0/mesh.V().field())*runTime.deltaTValue();

    meanCo0 =
        0.5*(gSum(sumPhi0)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi))().primitiveField()
    );

    Co = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanCo =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Base-Flow Courant Number mean: " << meanCo0
    << " max: " << Co0 << endl;

Info<< "Perturbation Courant Number mean: " << meanCo
    << " max: " << Co << endl;
	

// ************************************************************************* //

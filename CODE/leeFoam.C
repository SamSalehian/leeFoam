/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | https://urldefense.com/v3/__http://www.openfoam.com__;!!Pe9-MuB9P6Cr!BfDls6hCL059-0XJXkWiQuXDqA7XAjhR3P1uL2HtVDkPOHGlR1Qw3XRhk2AI4P17auM$ 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Patrick Good, Embry-Riddle Aeronautical University 
				  & Seyyed Salehian, Tuskegee University
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
    along with OpenFOAM.  If not, see <https://urldefense.com/v3/__http://www.gnu.org/licenses/__;!!Pe9-MuB9P6Cr!BfDls6hCL059-0XJXkWiQuXDqA7XAjhR3P1uL2HtVDkPOHGlR1Qw3XRhk2AIupn31DA$ >.

Application
    leeFoam

Description
    Solves the Linearized Euler Equations (LEE) for a aeroacoustics applications


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

	//Calculation maximun Acoustic Courant Number
	scalar  MaxCoA =
			max(mesh.surfaceInterpolation::deltaCoeffs()*c0).value()
			*runTime.deltaT().value();
	Info<< "Max Acoustic Courant Number = " << MaxCoA << endl;
	
	Info << nl << "Starting the time loop" << endl;
    while (runTime.loop()) 
    {
		Info << nl  << "Time = " << runTime.timeName() << endl;
		Info << "dt = " << runTime.deltaT().value() << endl;
		// Base-flow & Perturbation Cournant No are calculated 
		// and plotted in #include "leeCourantNo.H"
		#include "leeCourantNo.H"

		// fluxes defined in createFields.h
		// rhoPhi0_rho0 -> rho*U0/rho0
		// phi -> U
		// phi0 -> U0 
		
		// fluxes defined within the time loop
		// rhoPhi0 -> rho*U0
		// rho0Phi -> rho0*U
		surfaceScalarField rho0Phi
		(
			"rho0Phi",
			linearInterpolate(rho0*U) & mesh.Sf()
		);
		surfaceScalarField rPhi0
		(
			"rPhi0",
			linearInterpolate(U0) & mesh.Sf()
		);

		//Continuity eqn
		fvScalarMatrix contEqn
        ( 
           fvm::ddt(rho) 
		   + fvm::div(rPhi0, rho)
		   + fvc::div(rho0Phi)
        );
		contEqn.solve();
		
				
		//Momentum eqn
		fvVectorMatrix momEqn
		(
		   fvm::ddt(U)
		   + fvm::div(phi0,U)
		   + rho/rho0*fvc::div(phi0,U0)
		   + fvc::div(phi,U0)
		   + (1/rho0)*fvc::grad(p)
		   ==
		     -(D*U) + fvc::grad(D2*fvc::div(U)) 
		);		
		momEqn.solve();

		//Pressure eqn
		fvScalarMatrix PEqn //scalor matrix eqn???
		(
			fvm::ddt(p)
			+ (U & fvc::grad(p0)) 
			+ (gamma*p*fvc::div(U0))
			+ (U0 & fvc::grad(p))
			+ (gamma*p0*fvc::div(U))
		);	
		PEqn.solve();

		runTime.write();

	    // print execution time for preformance calculations  
	    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;
    }

    // printing execution time information at the end of simulation
    Info<< nl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

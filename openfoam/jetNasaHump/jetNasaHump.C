/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    NASAHumpBCs

Group
    no group

Description
    Preprocessor to create Dirichlet BCs for the NASA hump case    
\*---------------------------------------------------------------------------*/

//some OF forks do not fvCFD.H predefined, I have included an fvCFD_alt.H just in case
//comment out #include "fvCFD.H" and uncomment "include fvCFD_alt.H" 
#include "fvCFD.H" 
#include "simpleControl.H"
#include "fixedGradientFvPatchField.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting preprocessor time loop\n" << endl;
    
    // Update of the boundary patches - - - - - - - - - - - - - - - - - -// 
    const scalar Pi=constant::mathematical::pi;
    Info<<"Setting up Dirichlet BCs"<<endl<<endl;
    {
	// Remove singularities at the inlet (velocity and nuTilda)
	Info<<"Adding parabolic layer at the inlet extension"<<endl;
        
	// Find patch id for inlet
        label patch1=mesh.boundaryMesh().findPatchID("inlet");
 
        // Refcasting constant boundary fields
        fixedValueFvPatchVectorField& UPatch=
            refCast<fixedValueFvPatchVectorField>(U.boundaryFieldRef()[patch1]);
        vectorField& UField = UPatch;

        fixedValueFvPatchScalarField& nuTildaPatch=
            refCast<fixedValueFvPatchScalarField>(nuTilda.boundaryFieldRef()[patch1]);
        scalarField& nuTildaField = nuTildaPatch;

	// Check mass flow after profile is changed
	scalar massCheckInlet(0.0);

	// Loop over face centres to change BCs
        forAll(mesh.boundaryMesh()[patch1].faceCentres(),faceI)
        {
            scalar y=mesh.boundaryMesh()[patch1].faceCentres()[faceI].y();

            if (y<=0.02)
            {
                UField[faceI].component(0)=inletMag*2500*(0.04-y)*y;
                UField[faceI].component(1)=0.0;
                UField[faceI].component(2)=0.0;
            }
            else 
            {
                UField[faceI].component(0)=inletMag;
                UField[faceI].component(1)=0.0;
                UField[faceI].component(2)=0.0;
            }
            if (y<=0.002)
            {
                nuTildaField[faceI]=4.5822e-5*(0.004-y)*y/0.000004;
            }
            else 
            {
                nuTildaField[faceI]=4.5822e-5;
            }
            massCheckInlet += (UField[faceI] & mesh.Sf().boundaryField()[patch1][faceI]);
        }
        
	// Same process for jet, profile is harmonic (as in the paper)
        // Check mass flow for jet
	scalar massCheckJet(0.0);   

	// Jet patch id
        label patch2=mesh.boundaryMesh().findPatchID(jetName);
 
        // Refcasting velocity boundary fields
        fixedValueFvPatchVectorField& UPatch2=
            refCast<fixedValueFvPatchVectorField>(U.boundaryFieldRef()[patch2]);
        vectorField& UField2 = UPatch2;

        forAll(mesh.boundaryMesh()[patch2].faceCentres(),faceI)
        {
            scalar x=mesh.boundaryMesh()[patch2].faceCentres()[faceI].x();
            // Normalise for unknown jet width, symmetry needed
            scalar maxX=gMax(mesh.boundaryMesh()[patch2].faceCentres().component(0));
            scalar minX=gMin(mesh.boundaryMesh()[patch2].faceCentres().component(0));
            
	    x=(x-minX)/(maxX-minX);
            
	    // Find velocity magnitude
            scalar uMag(jetMag*(1.0-Foam::cos(2*x*Pi))/2);  
	    // Project on surface normal value (worked better with PGD)
            UField2[faceI] = uMag * mesh.Sf().boundaryField()[patch2][faceI]/mesh.magSf().boundaryField()[patch2][faceI];
	    // Zero out z-component (might be important for quasi-3D, can be removed alltogether)
            UField2[faceI].component(2)=0.0;
            massCheckJet = massCheckJet + (UField2[faceI] & mesh.Sf().boundaryField()[patch2][faceI]);
        }
  
	// Write fields and print mass flow rate
        U.write();
        nuTilda.write();
        Info<< "Inlet volume is "<<massCheckInlet<<"m^3/s"<<endl;   
        Info<< "Jet suction volume is "<<massCheckJet<<"m^3/s"<<endl;   

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
       
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

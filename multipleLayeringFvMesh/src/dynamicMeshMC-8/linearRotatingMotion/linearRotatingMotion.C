/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "linearRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearRotatingMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearRotatingMotion::linearRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
	velocity_(SBMFCoeffs_.lookup("velocity")),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    omega_(DataEntry<scalar>::New("omega", SBMFCoeffs_))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearRotatingMotion::~linearRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::linearRotatingMotion::transformation() const
{
    scalar t = time_.value();
	scalar dt = time_.deltaT().value();
			
    // Rotation around axis
    scalar angle = omega_->integrate(0, dt);
	//Reescribir el código con arreglos "old--" de forma de poder realizar movimientos relativos. Ahora solo se pude correr teniendo en cuenta de arrancar siempre en t=0.
	vector displacement;
	// Translation of centre of gravity with constant velocity
	if (t > 0)
	{    
		displacement = velocity_*(t-dt);
	}
	else
	{
		displacement = 0*velocity_;
	}

    quaternion R(axis_, angle);
	
	Info << "Tiempo: " << t << "origen: " << origin_ + displacement << endl;
    septernion TR(septernion(origin_ + displacement)*R*septernion(-origin_ - displacement));
		
	if(t>0)
	{
		const vector deltaD = velocity_*dt;
		Info << "Desplazamiento: " << deltaD << endl;
		TR.operator+=(deltaD);
	}
	

    Info<< "solidBodyMotionFunctions::linearRotatingMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::linearRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    omega_.reset
    (
        DataEntry<scalar>::New("omega", SBMFCoeffs_).ptr()
    );

    return true;
}

// ************************************************************************* //

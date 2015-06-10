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

#include "circularLinearRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(circularLinearRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        circularLinearRotatingMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::circularLinearRotatingMotion::circularLinearRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
	circularOmega_(readScalar(SBMFCoeffs_.lookup("circularOmega"))),
	radius_(readScalar(SBMFCoeffs_.lookup("radius"))),
    origin_(SBMFCoeffs_.lookup("origin")),
    axis_(SBMFCoeffs_.lookup("axis")),
    omega_(DataEntry<scalar>::New("omega", SBMFCoeffs_))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::circularLinearRotatingMotion::~circularLinearRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::circularLinearRotatingMotion::transformation() const
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
		displacement.x() = radius_ * cos(circularOmega_*(t-dt)) - radius_;
		displacement.y() = radius_ * sin(circularOmega_*(t-dt));
		displacement.z() = 0;
 	}
	else
	{
		displacement.x() = 0;
		displacement.y() = 0;
		displacement.z() = 0;
	}

    quaternion R(axis_, angle);
	
	Info <<"Tiempo: " << t << "origen: " << origin_ + displacement << endl;
    septernion TR(septernion(origin_ + displacement)*R*septernion(-origin_ - displacement));
	vector deltaD;
		
	if(t>0)
	{
		deltaD.x() = radius_ * cos(circularOmega_*(t)) - radius_ * cos(circularOmega_*(t-dt));
		deltaD.y() = radius_ * sin(circularOmega_*(t)) - radius_ * sin(circularOmega_*(t-dt));		
		deltaD.z() = 0;

		Info << "Desplazamiento: " << deltaD << endl;
		TR.operator+=(deltaD);
	}
	

    Info<< "solidBodyMotionFunctions::circularLinearRotatingMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::circularLinearRotatingMotion::read
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

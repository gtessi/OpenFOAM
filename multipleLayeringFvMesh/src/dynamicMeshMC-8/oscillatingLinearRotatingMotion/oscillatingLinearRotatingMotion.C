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

#include "oscillatingLinearRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"


using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(oscillatingLinearRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        oscillatingLinearRotatingMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingLinearRotatingMotion::oscillatingLinearRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
	linearAmplitude_(SBMFCoeffs_.lookup("linearAmplitude")),
	linearOmega_(readScalar(SBMFCoeffs_.lookup("linearOmega"))),
    origin_(SBMFCoeffs_.lookup("origin")),
    rotatingAmplitude_(SBMFCoeffs_.lookup("rotatingAmplitude")),
	rotatingOmega_(DataEntry<scalar>::New("rotatingOmega", SBMFCoeffs_))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::oscillatingLinearRotatingMotion::~oscillatingLinearRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::oscillatingLinearRotatingMotion::transformation() const
{
    scalar t = time_.value();
	scalar dt = time_.deltaT().value();
	scalar timeIndex = time_.timeIndex();
	scalar startTimeIndex = time_.startTimeIndex();
	scalar relativeTimeIndex = timeIndex -startTimeIndex;
			
	// Ojo al construir la malla ya esta girando!!!. 
//////////////////////////////////////////////////////////
    // Rotation around axis
	scalar angle;
	scalar oldAngle;
	
	if (relativeTimeIndex > 0)
	{
    	angle = rotatingOmega_->integrate(0, t);
		oldAngle = rotatingOmega_->integrate(0,(t-dt)); 
	}
	else
	{
		angle = 0;	
		oldAngle = 0;
	}
		
	vector eulerAngles = rotatingAmplitude_*(sin(angle)-sin(oldAngle));

	// Convert the rotational motion from deg to rad
    eulerAngles *= pi/180.0;
///////////////////////////////////////////////////////////
	//Reescribir el código con arreglos "old--" de forma de poder realizar movimientos relativos. Ahora solo se pude correr teniendo en cuenta de arrancar siempre en t=0.
	vector displacement;
	
	// Translation of centre of gravity with constant velocity
	if (relativeTimeIndex > 0)
	{    
		displacement = linearAmplitude_*sin(linearOmega_*(t-dt));
	}
	else
	{
		displacement = 0*linearAmplitude_;
	}

    //quaternion R(axis_, angle);
    quaternion R(eulerAngles.x(), eulerAngles.y(), eulerAngles.z());
	
	Info << "Tiempo: " << t << "origen: " << origin_ + displacement << endl;

    septernion TR(septernion(origin_ + displacement)*R*septernion(-origin_ - displacement));
		
	if(relativeTimeIndex > 0)
	{
		const vector deltaD = linearAmplitude_*sin(linearOmega_*(t)) -linearAmplitude_*sin(linearOmega_*(t-dt));
		TR.operator+=(deltaD);
	}
	

    Info<< "solidBodyMotionFunctions::oscillatingLinearRotatingMotion::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::oscillatingLinearRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    rotatingOmega_.reset
    (
        DataEntry<scalar>::New("rotatingOmega", SBMFCoeffs_).ptr()
    );

    return true;
}

// ************************************************************************* //

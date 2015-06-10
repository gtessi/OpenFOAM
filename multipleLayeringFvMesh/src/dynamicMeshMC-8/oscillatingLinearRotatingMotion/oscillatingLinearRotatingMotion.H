/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::solidBodyMotionFunctions::oscillatingLinearRotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The rotation is defined by an origin and axis of rotation and an angular
    speed.

SourceFiles
    oscillatingLinearRotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingLinearRotatingMotion_H
#define oscillatingLinearRotatingMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "DataEntry.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class oscillatingLinearRotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class oscillatingLinearRotatingMotion
:
    public solidBodyMotionFunction
{
    // Private data

		//- Amplitude of linear motion
        vector linearAmplitude_;

        //- Omega of linear motion
        scalar linearOmega_;

        //- Origin of the axis
        mutable vector origin_;
		
        //- Amplitude of rotating motion
        vector rotatingAmplitude_;

        //- Angular velocty (rad/sec)
        autoPtr<DataEntry<scalar> > rotatingOmega_;
          
	
    // Private Member Functions

        //- Disallow copy construct
        oscillatingLinearRotatingMotion(const oscillatingLinearRotatingMotion&);

        //- Disallow default bitwise assignment
        void operator=(const oscillatingLinearRotatingMotion&);


public:

    //- Runtime type information
    TypeName("oscillatingLinearRotatingMotion");


    // Constructors

        //- Construct from components
        oscillatingLinearRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new oscillatingLinearRotatingMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~oscillatingLinearRotatingMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

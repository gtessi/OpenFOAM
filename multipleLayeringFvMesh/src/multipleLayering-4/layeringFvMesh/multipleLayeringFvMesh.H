﻿/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::multipleLayeringFvMesh

Description
    Sample topoChangerFvMesh that moves a cell zone in a given direction,
    introducing/removing layers around it.

    Cells in the movingCells zone shall be moved given the prescribed velocity
    and will be bounded in "front" and "back" by other cell zones.
    Layer addition/removal interfaces are inserted at boundaries between the
    moving zone and front and back, pointing outside of the moving cell zone

Author and rewrite
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    multipleLayeringFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef multipleLayeringFvMesh_H
#define multipleLayeringFvMesh_H

#include "topoChangerFvMesh.H"
#include "solidBodyMotionFunction.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multipleLayeringFvMesh Declaration
\*---------------------------------------------------------------------------*/

class multipleLayeringFvMesh
:
    public topoChangerFvMesh
{
    // Private data
		

	//- Dynamic mesh dictionary
        dictionary dict_;

        //- Moving cell zone name
        List<word> movingCellsName_;

        //- Front face zone name
        List<word> frontFacesName_;

        //- Back face zone name
        List<word> backFacesName_;

        //- Moving point zone name
        List<word> movingPointsName_;

        //- Motion control function
        List<autoPtr<solidBodyMotionFunction> >  SBMFPtr_;

        //- Vertex motion mask
        scalarField motionMask_;

        PtrList<dictionary> groups_;
        
        //last displacement
        List<vector> xn_;

        int i;

        //Pair of modifiers counter
        int modK;


        //- Processors modifiers name list
        List<List<label> > procsModifiersList_;

        //- Processors groups name list
        List<List<label> > procsGroupsList_;

        //- Group name list
        List<label> groupsNameList_;

        //- Modifiers processors list
        List<List<label> > subRankList_;

        //- Communicator index list
        List<label> commIndexList_;

        //-----------------------------------------


    // Private Member Functions

        //- Disallow default bitwise copy construct
        multipleLayeringFvMesh(const multipleLayeringFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const multipleLayeringFvMesh&);


        //- Add mixer zones and modifiers
        void addZonesAndModifiers();

        //- Mark motion vertices
        tmp<scalarField> calcMotionMask(int) const;

        //- Select a pair of opposite layering modifiers and make it live, else -> disable.
        void selectPairLayersLive();


        //- Unique, sorted element list from list of lists
        List<label> elementsList(List<List<label> > list);

        //- Belong to list
        bool belongsTo(label val, List<label> list);

        //- Invert boolean list
        void invBoolList(List<bool>& list);



public:

    //- Runtime type information
    TypeName("multipleLayeringFvMesh");


    // Constructors

        //- Construct from database
        explicit multipleLayeringFvMesh(const IOobject& io);


    // Destructor

        virtual ~multipleLayeringFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

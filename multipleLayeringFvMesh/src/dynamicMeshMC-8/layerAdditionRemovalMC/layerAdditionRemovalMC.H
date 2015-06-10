/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::layerAdditionRemovalMC

Description
    Cell layer addition mesh modifier

SourceFiles
    layerAdditionRemovalMC.C
    addCellLayer.C
    removeCellLayer.C

\*---------------------------------------------------------------------------*/

#ifndef layerAdditionRemovalMC_H
#define layerAdditionRemovalMC_H

#include "polyMeshModifier.H"
#include "primitiveFacePatch.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class layerAdditionRemovalMC Declaration
\*---------------------------------------------------------------------------*/

class layerAdditionRemovalMC
:
    public polyMeshModifier
{
    // Private data

        //- Master face zone ID
        faceZoneID faceZoneID_;

	//-Remove fraction of extrusion layer.  Triggers layer removal
        mutable scalar removeFactor_;

        //- Add fraction of extrusion layer.  Triggers layer addition
        mutable scalar addFactor_;
		
	//- Fraction of the actual layer height for points offset.  Triggers layer addition
        mutable scalar addLayerFactor_;

        //- Layer thickness from previous step
        //  Used to decide whether to add or remove layers
        mutable scalar oldLayerThickness_;

        //- Point pairing
        mutable labelList* pointsPairingPtr_;

        //- Face pairing
        mutable labelList* facesPairingPtr_;

        //- Layer removal trigger time index
        mutable label triggerRemoval_;

        //- Layer addition trigger  time index
        mutable label triggerAddition_;


        //- Modifier name
        mutable label modifierName_;

        //- Group name
        mutable label groupName_;

        //- Communicator index
        mutable label commIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        layerAdditionRemovalMC(const layerAdditionRemovalMC&);

        //- Disallow default bitwise assignment
        void operator=(const layerAdditionRemovalMC&);

        //- Check validity of construction data
        void checkDefinition();


        // Topological changes

            //- Check for valid layer
            bool validCollapse() const;

            //- Set layer pairing.  Return true if a valid layer exists
            bool setLayerPairing() const;
           
            //- Return points pairing in a layer (not automatic!)
            const labelList& pointsPairing() const;

            //- Return faces pairing in a layer (not automatic!)
            const labelList& facesPairing() const;

            //- Calculate the offset to the next layer
            tmp<vectorField> extrusionDir() const;

            //- Add a layer of cells
            void addCellLayer(polyTopoChange&) const;

            //- Remove a layer of cells
            void removeCellLayer(polyTopoChange&) const;

            //- Clear addressing
            void clearAddressing() const;

        // Helpers

            //- Optionally read old thickness
            static scalar readOldThickness(const dictionary&);


    // Static data members

        //- Thickness insertion fraction for the pre-motion
        static const scalar addDelta_;

        //- Thickness removal fraction for the cell collapse
        //  Note: the cell will be collapsed to this relative
        //  thickness before the layer is removed.
        static const scalar removeDelta_;

public:

    //- Runtime type information
    TypeName("layerAdditionRemovalMC");


    // Constructors

        //- Construct from components
        layerAdditionRemovalMC
        (
            const word& name,
            const label index,
            const polyTopoChanger& mme,
            const word& zoneName,
            const scalar removeFactor,
            const scalar addFactor,           
            const label modifierName, // Modifier name
            const label groupName // Group name
        );

        //- Construct from dictionary
        layerAdditionRemovalMC
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    //- Destructor
    virtual ~layerAdditionRemovalMC();


    // Member Functions

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);


        // Edit

            //- Return min layer thickness which triggers removal
            scalar removeFactor() const
            {
                return removeFactor_;
            }

            //- Set min layer thickness which triggers removal
            void setRemoveFactor(const scalar t) const;

            //- Return max layer thickness which triggers removal
            scalar addFactor() const
            {
                return addFactor_;
            }

            //- Set max layer thickness which triggers removal
            void setAddFactor(const scalar t) const;


            //- Return modifier name
            label modifierName() const
            {
               return modifierName_;
            }

            //- Return group name
            label groupName() const
            {
               return groupName_;
            }

            //- Load communicators index list
            void loadCommunicators(label commIndex) const
            {
                commIndex_ = commIndex;
            }

            //- Return master face zone ID
            faceZoneID getFaceZoneID() const
            {
                return faceZoneID_;
            }


        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

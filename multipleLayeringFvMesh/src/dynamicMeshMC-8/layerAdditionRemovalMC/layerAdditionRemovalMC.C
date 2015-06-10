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

Description
    Cell layer addition/removal mesh modifier

\*---------------------------------------------------------------------------*/

#include "layerAdditionRemovalMC.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "Time.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "oppositeFace.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(layerAdditionRemovalMC, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        layerAdditionRemovalMC,
        dictionary
    );
}


const Foam::scalar Foam::layerAdditionRemovalMC::addDelta_ = 0.5;
const Foam::scalar Foam::layerAdditionRemovalMC::removeDelta_ = 0.1;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layerAdditionRemovalMC::checkDefinition()
{
    /*if (!faceZoneID_.active())
    {
        FatalErrorIn
        (
            "void Foam::layerAdditionRemovalMC::checkDefinition()"
        )   << "Master face zone named " << faceZoneID_.name()
            << " cannot be found."
            << abort(FatalError);
    }

    if
    (
        minLayerThickness_ < VSMALL
     || maxLayerThickness_ < minLayerThickness_
    )
    {
        FatalErrorIn
        (
            "void Foam::layerAdditionRemovalMC::checkDefinition()"
        )   << "Incorrect layer thickness definition."
            << abort(FatalError);
    }*/

  /*  if (topoChanger().mesh().faceZones()[faceZoneID_.index()].empty())
    {
        FatalErrorIn
        (
            "void Foam::layerAdditionRemovalMC::checkDefinition()"
        )   << "Face extrusion zone contains no faces. "
            << " Please check your mesh definition."
            << abort(FatalError);
    }*/

    if (debug)
    {
        Pout<< "Cell layer addition/removal object " << name() << " :" << nl
            << "    faceZoneID: " << faceZoneID_ << endl;
    }
}

Foam::scalar Foam::layerAdditionRemovalMC::readOldThickness
(
    const dictionary& dict
)
{
    return dict.lookupOrDefault("oldLayerThickness", -1.0);
}


void Foam::layerAdditionRemovalMC::clearAddressing() const
{
    // Layer removal data
    deleteDemandDrivenData(pointsPairingPtr_);
    deleteDemandDrivenData(facesPairingPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::layerAdditionRemovalMC::layerAdditionRemovalMC
(const word& name,
    const label index,
    const polyTopoChanger& mme,
    const word& zoneName,
    const scalar removeFactor,
    const scalar addFactor,

    const label modifierName, // Modifier name
    const label groupName // Group name
)
:
    polyMeshModifier(name, index, mme, true),
    faceZoneID_(zoneName, mme.mesh().faceZones()),
    removeFactor_(removeFactor),
    addFactor_(addFactor),
    addLayerFactor_(),

    modifierName_(modifierName), // Modifier name
    groupName_(groupName), // Group name

    oldLayerThickness_(-1.0),
    pointsPairingPtr_(NULL),
    facesPairingPtr_(NULL),
    triggerRemoval_(-1),
    triggerAddition_(-1)
{
    checkDefinition();
}


// Construct from dictionary
Foam::layerAdditionRemovalMC::layerAdditionRemovalMC
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, Switch(dict.lookup("active"))),
    faceZoneID_(dict.lookup("faceZoneName"), mme.mesh().faceZones()),
    removeFactor_(readScalar(dict.lookup("removeFactor"))),
    addFactor_(readScalar(dict.lookup("addFactor"))),
    addLayerFactor_(),
    oldLayerThickness_(readOldThickness(dict)),

    modifierName_(readScalar(dict.lookup("modifierName"))), // Modifier name
    groupName_(readScalar(dict.lookup("groupName"))), // Group name

    pointsPairingPtr_(NULL),
    facesPairingPtr_(NULL),
    triggerRemoval_(-1),
    triggerAddition_(-1)
{
    checkDefinition();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::layerAdditionRemovalMC::~layerAdditionRemovalMC()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::layerAdditionRemovalMC::changeTopology() const
{
    // Protect from multiple calculation in the same time-step
    if (triggerRemoval_ > -1 || triggerAddition_ > -1)
    {
        return true;
    }

    // Go through all the cells in the master layer and calculate
    // approximate layer thickness as the ratio of the cell volume and
    // face area in the face zone.
    // Layer addition:
    //     When the max thickness exceeds the threshold, trigger refinement.
    // Layer removal:
    //     When the min thickness falls below the threshold, trigger removal.

    const faceZone& fz = topoChanger().mesh().faceZones()[faceZoneID_.index()];
    const labelList& mc = fz.masterCells();

    const scalarField& V = topoChanger().mesh().cellVolumes();
    const vectorField& S = topoChanger().mesh().faceAreas();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Nuevo, para Layering-MM, segundo layer
    const polyMesh& mesh = topoChanger().mesh();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    labelList secondLayerFaceLabels(mc.size());
    labelList secondLayerCellLabels(mc.size());

    
    forAll(mc, faceI)
    {
        secondLayerFaceLabels[faceI] = cells[mc[faceI]].opposingFaceLabel(fz[faceI], faces);
        secondLayerCellLabels[faceI] = (mesh.faceOwner()[secondLayerFaceLabels[faceI]]==mc[faceI]) ? mesh.faceNeighbour()[secondLayerFaceLabels[faceI]] : mesh.faceOwner()[secondLayerFaceLabels[faceI]];
    }


    scalar avgDeltaSecond = 0;
    scalar minDeltaSecond = GREAT;
    scalar maxDeltaSecond = 0;

    forAll(secondLayerFaceLabels, faceI)
    {
        scalar curDeltaSecond = V[secondLayerCellLabels[faceI]]/mag(S[secondLayerFaceLabels[faceI]]);
        avgDeltaSecond += curDeltaSecond;
        minDeltaSecond = min(minDeltaSecond, curDeltaSecond);
        maxDeltaSecond = max(maxDeltaSecond, curDeltaSecond);
    }

    avgDeltaSecond /= secondLayerFaceLabels.size();
	Info << "Second Layer height " << avgDeltaSecond << endl;

	//cellList secondLayerCells = secondLayerFaceList.faceOwner();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (min(V) < -VSMALL)
    {
        FatalErrorIn("bool layerAdditionRemovalMC::changeTopology() const")
            << "negative cell volume. Error in mesh motion before "
            << "topological change.\n V: " << V
            << abort(FatalError);
    }

    scalar avgDelta = 0;
    scalar minDelta = GREAT;
    scalar maxDelta = 0;
	scalar avgPosFz = 0;	
	

    forAll(fz, faceI)
    {
        scalar curDelta = V[mc[faceI]]/mag(S[fz[faceI]]);
        avgDelta += curDelta;
        minDelta = min(minDelta, curDelta);
        maxDelta = max(maxDelta, curDelta);
    }

    avgDelta /= fz.size();
    avgPosFz /= fz.size();

    addLayerFactor_ = (avgDelta - avgDeltaSecond)/avgDelta;
    Info << "addLayerFactor_ value: " << addLayerFactor_ << endl;
  
    if (debug)
    {
        Pout<< "bool layerAdditionRemovalMC::changeTopology() const "
            << " for object " << name() << " : " << nl
            << "Layer thickness: min: " << minDelta
            << " max: " << maxDelta << " avg: " << avgDelta
            << " old thickness: " << oldLayerThickness_ << nl
            << "Removal threshold: " << removeFactor_
            << " addition threshold: " << addFactor_ << endl;
    }

    bool topologicalChange = false;

    // If the thickness is decreasing and crosses the min thickness,
    // trigger removal
    if (oldLayerThickness_ < 0)
    {
        if (debug)
        {
            Pout<< "First step. No addition/removal" << endl;
        }

        // No topological changes allowed before first mesh motion
        //
        oldLayerThickness_ = avgDelta;

        topologicalChange = false;
    }
    else if (avgDelta < oldLayerThickness_) // REMOVAL
    {
        /**
        // *** debug info
        Pout << "** REMOVAL  - modifier Name: " << modifierName_ << " - commIndex_: " << commIndex_ << nl;
        **/

        bool localTopologicalChange;

        if (minDelta < removeFactor_*avgDeltaSecond)
        {
            localTopologicalChange = true;
        }
        else
        {
            localTopologicalChange = false;
        }

        bool globalTopologicalChange = localTopologicalChange;

        /**
        // *** debug info
        Pout << "*** commIndex_ " << commIndex_ << nl;
        **/

        // Applies reduce operation only if running in parallel
        if (Pstream::parRun())
        {
            reduce(globalTopologicalChange, andOp<bool>(), 0, commIndex_);
        }

        // Layers moving towards addition
        if (globalTopologicalChange)
        {
            // Check layer pairing
            if (setLayerPairing())
            {
                // A mesh layer detected.  Check that collapse is valid
                if (validCollapse())
                {
                    // At this point, info about moving the old mesh
                    // in a way to collapse the cells in the removed
                    // layer is available.  Not sure what to do with
                    // it.

                    if (debug)
                    {
                        Pout<< "bool layerAdditionRemovalMC::changeTopology() "
                            << " const for object " << name() << " : "
                            << "Triggering layer removal" << endl;
                    }

                    triggerRemoval_ = topoChanger().mesh().time().timeIndex();

                    // Old thickness looses meaning.
                    // Set it up to indicate layer removal
                    oldLayerThickness_ = GREAT;

                    topologicalChange = true;
                }
                else
                {
                    // No removal, clear addressing
                    clearAddressing();
                }
            }
        }
        else
        {
            oldLayerThickness_ = avgDelta;
        }

        
    }
    else // ADDITION
    {
        /**
        // *** debug info
        Pout << "** ADDITION - modifier Name: " << modifierName_ << " - commIndex_: " << commIndex_ << nl;
        **/

        bool localTopologicalChange;

        if (maxDelta > addFactor_*avgDeltaSecond)
		{
			localTopologicalChange = true;
		}
		else
		{
			localTopologicalChange = false;
		}

		bool globalTopologicalChange = localTopologicalChange;

        /**
        // *** debug info
        Pout << "*** commIndex_ " << commIndex_ << nl;
        **/

        // Applies reduce operation only if running in parallel
        if (Pstream::parRun())
        {
            reduce(globalTopologicalChange, andOp<bool>(), 0, commIndex_);
        }

		// Layers moving towards addition
        if (globalTopologicalChange)
        {
            if (debug)
            {
                Pout<< "bool layerAdditionRemovalMC::changeTopology() const "
                    << " for object " << name() << " : "
                    << "Triggering layer addition" << endl;
            }

            triggerAddition_ = topoChanger().mesh().time().timeIndex();

            // Old thickness looses meaning.
            // Set it up to indicate layer removal
            oldLayerThickness_ = 0;

            topologicalChange = true;
        }
        else
        {
            oldLayerThickness_ = avgDelta;
        }
    }

    return topologicalChange;
}


void Foam::layerAdditionRemovalMC::setRefinement(polyTopoChange& ref) const
{
    // Insert the layer addition/removal instructions
    // into the topological change

    if (triggerRemoval_ == topoChanger().mesh().time().timeIndex())
    {
        removeCellLayer(ref);

        // Clear addressing.  This also resets the addition/removal data
        if (debug)
        {
            Pout<< "layerAdditionRemovalMC::setRefinement(polyTopoChangeMC& ref) "
                << " for object " << name() << " : "
                << "Clearing addressing after layer removal. " << endl;
        }

        triggerRemoval_ = -1;
        clearAddressing();
    }

    if (triggerAddition_ == topoChanger().mesh().time().timeIndex())
    {
		//Pout << "Estoy por entrar a addCellLayer" << endl;
        addCellLayer(ref);

        // Clear addressing.  This also resets the addition/removal data
        if (debug)
        {
            Pout<< "layerAdditionRemovalMC::setRefinement(polyTopoChangeMC& ref) "
                << " for object " << name() << " : "
                << "Clearing addressing after layer addition. " << endl;
        }

        triggerAddition_ = -1;
        clearAddressing();
    }
}


void Foam::layerAdditionRemovalMC::updateMesh(const mapPolyMesh&)
{
    if (debug)
    {
        Pout<< "layerAdditionRemovalMC::updateMesh(const mapPolyMesh&) "
            << " for object " << name() << " : "
            << "Clearing addressing on external request. ";

        if (pointsPairingPtr_ || facesPairingPtr_)
        {
            Pout<< "Pointers set." << endl;
        }
        else
        {
            Pout<< "Pointers not set." << endl;
        }
    }

    // Mesh has changed topologically.  Update local topological data
    faceZoneID_.update(topoChanger().mesh().faceZones());

    clearAddressing();
}


void Foam::layerAdditionRemovalMC::setRemoveFactor(const scalar t) const
{
    if
    (
        t < VSMALL
     || addFactor_ < t
    )
    {
        FatalErrorIn
        (
            "void layerAdditionRemovalMC::setMinLayerThickness("
            "const scalar t) const"
        )   << "Incorrect layer thickness definition."
            << abort(FatalError);
    }

    removeFactor_ = t;
}


void Foam::layerAdditionRemovalMC::setAddFactor(const scalar t) const
{
    if
    (
        t < removeFactor_
    )
    {
        FatalErrorIn
        (
            "void layerAdditionRemovalMC::setMaxLayerThickness("
            "const scalar t) const"
        )   << "Incorrect layer thickness definition."
            << abort(FatalError);
    }

    addFactor_ = t;
}



void Foam::layerAdditionRemovalMC::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name()<< nl
        << faceZoneID_ << nl
        << oldLayerThickness_ << nl
        << addLayerFactor_ << nl //<< endl;
        << modifierName_ << nl	// Print modifier name
        << groupName_ << endl;	// Print group name
}


void Foam::layerAdditionRemovalMC::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << "    faceZoneName " << faceZoneID_.name()
        << token::END_STATEMENT << nl
        << "    removeFactor " << removeFactor_
        << token::END_STATEMENT << nl
        << "    addFactor " << addFactor_
        << token::END_STATEMENT << nl
        << "    oldLayerThickness " << oldLayerThickness_
        << token::END_STATEMENT << nl
        << "    addFactor " << addLayerFactor_
        << token::END_STATEMENT << nl
        << "    modifierName " << modifierName_	// Print modifier name
        << token::END_STATEMENT << nl
        << "    groupName " << groupName_	// Print group name
        << token::END_STATEMENT << nl
        << "    active " << active()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// ************************************************************************* //

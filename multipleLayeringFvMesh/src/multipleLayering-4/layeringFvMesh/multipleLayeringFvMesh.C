/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

#include "multipleLayeringFvMesh.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "layerAdditionRemovalMC.H"
#include "volMesh.H"
#include "transformField.H"
#include "addToRunTimeSelectionTable.H"
#include "instantList.H"
#include "linearRotatingMotion.H"
#include "oscillatingLinearRotatingMotion.H"
#include "circularLinearRotatingMotion.H"

#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multipleLayeringFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        multipleLayeringFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::tmp<Foam::scalarField>
Foam::multipleLayeringFvMesh::calcMotionMask(int i) const
{
    Info << "Updating vertex markup" << endl;

    //tmp<scalarField> tvertexMarkup(new scalarField(allPoints().size(), 0));
    tmp<scalarField> tvertexMarkup(new scalarField(points().size(), 0));
    scalarField& vertexMarkup = tvertexMarkup();
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    pointZoneID movingPointsID(movingPointsName_[i], pointZones());

    if (movingPointsID.active())
    {
        // Get labels of all moving points
        const labelList& movingPoints = pointZones()[movingPointsID.index()];
        
        forAll (movingPoints, pointI)
        {
            vertexMarkup[movingPoints[pointI]] = 1;
        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cellZoneID movingCellsID(movingCellsName_[i], cellZones());

    // In order to do a correct update on a mask on processor boundaries,
    // Detection of moving cells should use patchNeighbourField for
    // processor (not coupled!) boundaries.  This is done by expanding
    // a moving cell set into a field and making sure that processor patch
    // points move in sync.  Not done at the moment, probably best to do
    // using parallel update of pointFields.  HJ, 19/Feb/2011

    // If moving cells are found, perform mark-up
    if (movingCellsID.active())
    {
        // Get cell-point addressing
        const labelListList& cp = cellPoints();

        // Get labels of all moving cells
        const labelList& movingCells = cellZones()[movingCellsID.index()];

        forAll (movingCells, cellI)
        {
            const labelList& curCp = cp[movingCells[cellI]];

            forAll (curCp, pointI)
            {
                vertexMarkup[curCp[pointI]] = 1;
            }
        }
    }

    faceZoneID frontFacesID(frontFacesName_[i], faceZones());

    if (frontFacesID.active())
    {
        const faceZone& frontFaces = faceZones()[frontFacesID.index()];

        const labelList& mp = frontFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }

    faceZoneID backFacesID(backFacesName_[i], faceZones());

    if (backFacesID.active())
    {
        const faceZone& backFaces = faceZones()[backFacesID.index()];

        const labelList& mp = backFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }

    return tvertexMarkup;
}


void Foam::multipleLayeringFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

   /* if (topoChanger_.size() > 0)
    {
        Info<< "void multipleLayeringFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }*/

    // Add layer addition/removal interfaces
    topoChanger_.setSize(2*groups_.size());
    label nMods = 0;

    forAll(groups_,j)
    {
            faceZoneID frontFacesID(frontFacesName_[j], faceZones());
            faceZoneID backFacesID(backFacesName_[j], faceZones());

            if (frontFacesID.active())
            {
                const faceZone& frontFaces = faceZones()[frontFacesID.index()];

                if (!frontFaces.empty())
                {
                    topoChanger_.set
                    (
                        nMods,
                        new layerAdditionRemovalMC
                        (
                            frontFacesName_[j] + "Layer",
                            nMods,
                            topoChanger_,
                            frontFacesName_[j],
                            readScalar
				            (
				                groups_[j].subDict("front").lookup("removeFactor")
				            ),
				            readScalar
				            (
				                groups_[j].subDict("front").lookup("addFactor")
				            ),
                            // Front modifier value
                            readScalar
                            (
                                groups_[j].subDict("front").lookup("frontModifierValue")
                            ),
                            // Group modifier value
                            readScalar
                            (
                                groups_[j].subDict("front").lookup("groupValue")
                            )
                        )
                    );

                    nMods++;
                }
            }

            if (backFacesID.active())
            {
                const faceZone& backFaces = faceZones()[backFacesID.index()];

                if (!backFaces.empty())
                {
                    topoChanger_.set
                    (
                        nMods,
                        new layerAdditionRemovalMC
                        (
                            backFacesName_[j] + "Layer",
                            nMods,
                            topoChanger_,
                            backFacesName_[j],
                            readScalar
				            (
				                groups_[j].subDict("back").lookup("removeFactor")
				            ),
				            readScalar
				            (
				                groups_[j].subDict("back").lookup("addFactor")
				            ),
                            readScalar
                            (
                                // Front modifier value
                                groups_[j].subDict("back").lookup("frontModifierValue")
                            ),
                            // Group modifier value
                            readScalar
                            (
                                groups_[j].subDict("back").lookup("groupValue")
                            )
                        )
                    );

                    nMods++;
                }
            }
    }

    topoChanger_.setSize(nMods);

    reduce(nMods, sumOp<label>());

    Info << "Adding " << nMods << " mesh modifiers" << endl;

    // Write mesh and modifiers
    topoChanger_.write();

    // No need to write the mesh - only modifiers are added.
    // HJ, 18/Feb/2011
//     write();frontModifierValue
}

void Foam::multipleLayeringFvMesh::selectPairLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    List<label> localGroupList = procsGroupsList_[Pstream::myProcNo()];

    /**
    // *** debug info
    Pout << "** localGroupList: " << localGroupList << nl;
    **/

    // Enable modifiers for layering
    forAll (localGroupList, groupI)
    {
        /**
        // *** debug info
        Pout << "** Local modifier: " << procsModifiersList_[Pstream::myProcNo()] << ", group: " << localGroupList[groupI] << nl;
        **/

        if (localGroupList[groupI] == groupsNameList_[modK])
        {
            topoChanges[groupI].enable();
        }
        else
        {
            topoChanges[groupI].disable();
        }
    }

    // Update group current index
    int size = groupsNameList_.size()-1;

    if (modK == size)
    {
        modK = 0;
    }
    else
    {
        ++modK;
    }

}



Foam::List<Foam::label> Foam::multipleLayeringFvMesh::elementsList(List<List<label> > list)
{
    SortableList<label> tmp;

    // concatenates and flats the list
    forAll(list, listI)
    {
        forAll(list[listI], listJ)
        {
            tmp.append(list[listI][listJ]);
        }
    }

    // sorts into ascending order
    tmp.sort();

    // gets begin and end iterators
    UList<label>::iterator pFirst = tmp.begin();
    UList<label>::iterator pLast = tmp.end();

    // removes repeated elements
    UList<label>::iterator new_pLast = std::unique(pFirst,pLast);

    // resizes the list
    tmp.resize(std::distance(pFirst,new_pLast));

    return tmp;
}



bool Foam::multipleLayeringFvMesh::belongsTo(label val, List<label> list)
{
    forAll(list, listI)
    {
        if (val == list[listI]) {
            return true;
        }
    }

    return false;
}



void Foam::multipleLayeringFvMesh::invBoolList(List<bool>& list)
{
    forAll(list, listI)
    {
        list[listI] = !list[listI];
    }
}







// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::multipleLayeringFvMesh::multipleLayeringFvMesh(const IOobject& io)
: 
    topoChangerFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    
    motionMask_(),
    groups_(dict_.lookup("groups"))
{
    movingCellsName_.resize(groups_.size());
    frontFacesName_.resize(groups_.size());
    backFacesName_.resize(groups_.size());
    movingPointsName_.resize(groups_.size());
    SBMFPtr_.resize(groups_.size());
    xn_.resize(groups_.size());

    forAll(groups_, j)
    {
        const dictionary& dict = groups_[j];
        const word& setMovingCells(dict.lookup("movingCells"));
        const word& setFrontFaces(dict.lookup("frontFaces"));
        const word& setBackFaces(dict.lookup("backFaces"));
        const word& setMovingPoints(dict.lookup("movingPoints"));
        movingCellsName_[j] = setMovingCells;
        frontFacesName_[j] = setFrontFaces;
        backFacesName_[j] = setBackFaces;     
        movingPointsName_[j] = setMovingPoints;
        SBMFPtr_[j] = solidBodyMotionFunction::New(dict, time());
        xn_[j] = SBMFPtr_[j]().transformation().t();
    }
    motionMask_ = calcMotionMask(0);
    addZonesAndModifiers();
    modK = 0;


    // *** layering running in parallel with modifier's name method ***

    // Resize
    procsModifiersList_.resize(Pstream::nProcs()); // Modifiers
    procsGroupsList_.resize(Pstream::nProcs()); // Groups

    // Modifier's name and group in each processor
    const PtrList<polyMeshModifier>& topoChanges = topoChanger_;

    forAll(topoChanges, modI)
    {
        const polyMeshModifier& currentPMM = topoChanges[modI];

        const layerAdditionRemovalMC& currentLAR = refCast<const layerAdditionRemovalMC>(currentPMM);


        procsModifiersList_[Pstream::myProcNo()].append(currentLAR.modifierName()); // Modifiers
        procsGroupsList_[Pstream::myProcNo()].append(currentLAR.groupName()); // Groups

        /**
        // *** debug info
        Pout << currentLAR.modifierName() << nl;
        Pout << currentLAR.groupName() << nl;
        **/
    }

    Pstream::gatherList(procsModifiersList_);
    Pstream::gatherList(procsGroupsList_);

    /**
    // *** debug info
    Info << " *** entrando en los procsModifiersList_ " << nl << nl;

    Pout << procsModifiersList_ << nl;

    Info << " -------------------- " << nl;

    Info << " *** entrando en los procsGroupsList_ " << nl << nl;

    Pout << procsGroupsList_ << nl;

    Info << " -------------------- " << nl;
    **/


    // Processors in each modifier
    if (Pstream::master()) {
        // Make flat, unique, sorted list
        List<label> modifiersNameList = elementsList(procsModifiersList_);
        groupsNameList_ = elementsList(procsGroupsList_);

        /**
        // *** debug info
        Pout << "** Modifiers name list: " << modifiersNameList << nl;
        Pout << "** Groups name list" << groupsNameList_ << nl;
        Info << " -------------------- " << nl;
        **/

        // Resize
        subRankList_.resize(modifiersNameList.size());

        // Go through each modifier name
        forAll(modifiersNameList, modI)
        {
            // Go through each processor
            forAll(procsModifiersList_, procModI)
            {
                if (belongsTo(modifiersNameList[modI], procsModifiersList_[procModI])) {
                    subRankList_[modI].append(procModI);
                }
            }
        }

        /**
        // *** debug info
        Pout << subRankList_ << nl;
        Info << " -------------------- " << nl;
        **/
    }


    // Scatter subRankList_ to all processors
    Pstream::scatter(subRankList_);

    // Allocate communicators for each modifier
    forAll(subRankList_, rankI)
    {
        commIndexList_.append(Pstream::allocateCommunicator(0, subRankList_[rankI]));
    }

    /**
    // *** debug info
    Pout << commIndexList_ << nl;
    Info << " -------------------- " << nl;
    **/

    // Send the communicators list
    forAll(topoChanges, modI)
    {
        const polyMeshModifier& currentPMM = topoChanges[modI];

        const layerAdditionRemovalMC& currentLAR = refCast<const layerAdditionRemovalMC>(currentPMM);

        currentLAR.loadCommunicators(commIndexList_[currentLAR.modifierName()]);
    }


    // Scatter groupsNameList_ to all processors
    Pstream::scatter(groupsNameList_);


}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multipleLayeringFvMesh::~multipleLayeringFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multipleLayeringFvMesh::update()
{
    pointField oldPointsNew = points();
    pointField newPoints = points();

    bool localMeshChanged;

    // Loop over all the groups during the current timestep
    forAll(groupsNameList_, groupI) {

        //Select correct pair of modifiers in order to do update mesh.
        selectPairLayersLive();

        /**
        // *** debug info
        const PtrList<polyMeshModifier>& topoChanges = topoChanger_;

        forAll (topoChanges, morphI)
        {
           Pout << "** Modifier named " << topoChanges[morphI].name() << " - status " << topoChanges[morphI].active() << nl;
        }
        **/

        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh(true); //(true, false);

        //bool localMeshChanged = topoChangeMap.valid();
        localMeshChanged = topoChangeMap.valid();

        bool globalMeshChanged = localMeshChanged;

        reduce(globalMeshChanged, orOp<bool>());


        if (globalMeshChanged)
        {
    //         // Map old points onto the new mesh
    //         pointField mappedOldPointsNew(allPoints().size());
    //         mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());

    //         movePoints(mappedOldPointsNew);
    //         resetMotion();
    //         setV0();

            // Get new points from preMotion
            Info << "Local Topology change. Calculating premotion points" << endl;
            newPoints = topoChangeMap().preMotionPoints();
        }

        forAll(groups_,j)
        {
            if(isA<Foam::solidBodyMotionFunctions::linearRotatingMotion>(SBMFPtr_[j]()) || isA<Foam::solidBodyMotionFunctions::oscillatingLinearRotatingMotion>(SBMFPtr_[j]()) || isA<Foam::solidBodyMotionFunctions::circularLinearRotatingMotion>(SBMFPtr_[j]()) )
            {

                ////////////////////////////////////////////////////////////////////////////////
                pointZoneID movingPointsID(movingPointsName_[j], pointZones());

                // Get labels of all moving (rotating) points
                const labelList& pointIDs_ = pointZones()[movingPointsID.index()];
                ////////////////////////////////////////////////////////////////////////////////
                 /*cellZoneID movingCellsID(movingCellsName_[i], cellZones());

                // Get cell-point addressing
                const labelListList& cp = cellPoints();

                // Get labels of all moving cells
                const labelList& movingCells = cellZones()[movingCellsID.index()];
                labelList pointIDs_;

                forAll (movingCells, cellI)
                {
                    //const labelList& curCp = cp[movingCells[cellI]];
                    const UList<label>& curCp = cp[movingCells[cellI]];
                    pointIDs_.append(curCp);
                }*/
                ////////////////////////////////////////////////////////////////////////////////


                //Info << "Por rotar" << endl;
                UIndirectList<point>(newPoints, pointIDs_) =
                transform
                (
                    SBMFPtr_[j]().transformation(),
                    pointField(newPoints, pointIDs_)
                );


            }

            else
            {
                motionMask_ = calcMotionMask(j);
                scalar tn = time().value();
                scalar dt = time().deltaT().value();
                Foam::vector xn1 = SBMFPtr_[j]().transformation().t();
                Foam::septernion velocity( (xn1-this->xn_[j])/dt, quaternion::zero );
                Info << "Executing mesh motion with velocity " << velocity << " and dt: " << dt << endl;
                newPoints += motionMask_*transform(velocity, newPoints)*dt;
                this->xn_[j] = xn1;
            }


        }

        movePoints(newPoints);
    }
    
    return localMeshChanged;      
    
}


// ************************************************************************* //


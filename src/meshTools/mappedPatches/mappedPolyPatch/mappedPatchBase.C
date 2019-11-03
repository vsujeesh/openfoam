/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearchMeshObject.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "treeDataFace.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "Time.H"
#include "mapDistribute.H"
#include "SubField.H"
#include "triPointRef.H"
#include "syncTools.H"
#include "treeDataCell.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBase, 0);
}


const Foam::Enum
<
    Foam::mappedPatchBase::sampleMode
>
Foam::mappedPatchBase::sampleModeNames_
({
    { sampleMode::NEARESTCELL, "nearestCell" },
    { sampleMode::NEARESTPATCHFACE, "nearestPatchFace" },
    { sampleMode::NEARESTPATCHFACEAMI, "nearestPatchFaceAMI" },
    { sampleMode::NEARESTPATCHPOINT, "nearestPatchPoint" },
    { sampleMode::NEARESTFACE, "nearestFace" },
    { sampleMode::NEARESTONLYCELL, "nearestOnlyCell" },
});


const Foam::Enum
<
    Foam::mappedPatchBase::offsetMode
>
Foam::mappedPatchBase::offsetModeNames_
({
    { offsetMode::UNIFORM, "uniform" },
    { offsetMode::NONUNIFORM, "nonuniform" },
    { offsetMode::NORMAL, "normal" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::mappedPatchBase::communicator
(
    const word& myWorld,
    const word& sampleWorld
)
{
    // Start off with local world
    label comm = UPstream::worldComm;

    if (!sampleWorld.empty() && !myWorld.empty())
    {
        const wordList& procWorlds = UPstream::worlds();

        if (!procWorlds.found(sampleWorld))
        {
            FatalErrorInFunction << "Cannot find sampleWorld " << sampleWorld
                << " in set of worlds " << procWorlds
                << exit(FatalError);
        }

        DynamicList<label> subRanks(procWorlds.size());
        forAll(procWorlds, proci)
        {
            const word& world = procWorlds[proci];
            if (world == myWorld || world == sampleWorld)
            {
                subRanks.append(proci);
            }
        }

        // Allocate new communicator with parent 0 (= world)
        comm = UPstream::allocateCommunicator(0, subRanks, true);

        Pout<< "*** myWorld:" << myWorld << " sampleWorld:" << sampleWorld
            << " using subRanks:" << subRanks << " new comm:" << comm << endl;
    }

    return comm;
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::facePoints
(
    const polyPatch& pp
) const
{
    const polyMesh& mesh = pp.boundaryMesh().mesh();

    // Force construction of min-tet decomp
    (void)mesh.tetBasePtIs();

    // Initialise to face-centre
    tmp<pointField> tfacePoints(new pointField(patch_.size()));
    pointField& facePoints = tfacePoints.ref();

    forAll(pp, facei)
    {
        facePoints[facei] = facePoint
        (
            mesh,
            pp.start()+facei,
            polyMesh::FACE_DIAG_TRIS
        ).rawPoint();
    }

    return tfacePoints;
}


void Foam::mappedPatchBase::collectSamples
(
    const label mySampleWorld,      // Wanted world
    const pointField& facePoints,
    pointField& samples,            // All samples
    labelList& patchFaceWorlds,     // Per sample: wanted world
    labelList& patchFaceProcs,      // Per sample: originating processor 
    labelList& patchFaces,          // Per sample: originating patchFace index
    pointField& patchFc             // Per sample: originating centre
) const
{
    const label oldComm(Pstream::warnComm);
    Pstream::warnComm = comm_;

    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    // Collect all sample points and the faces they come from.
    {
        List<pointField> globalFc(nProcs);
        globalFc[myRank] = facePoints;
        Pstream::gatherList(globalFc, Pstream::msgType(), comm_);
        Pstream::scatterList(globalFc, Pstream::msgType(), comm_);
        // Rework into straight list
        patchFc = ListListOps::combine<pointField>
        (
            globalFc,
            accessOp<pointField>()
        );
    }

    {
        List<pointField> globalSamples(nProcs);
        globalSamples[myRank] = samplePoints(facePoints);
        Pstream::gatherList(globalSamples, Pstream::msgType(), comm_);
        Pstream::scatterList(globalSamples, Pstream::msgType(), comm_);
        // Rework into straight list
        samples = ListListOps::combine<pointField>
        (
            globalSamples,
            accessOp<pointField>()
        );
    }

    {
        labelListList globalFaces(nProcs);
        globalFaces[myRank] = identity(patch_.size());
        // Distribute to all processors
        Pstream::gatherList(globalFaces, Pstream::msgType(), comm_);
        Pstream::scatterList(globalFaces, Pstream::msgType(), comm_);

        patchFaces = ListListOps::combine<labelList>
        (
            globalFaces,
            accessOp<labelList>()
        );
    }

    {
        labelList procToWorldIndex(nProcs);
        procToWorldIndex[myRank] = mySampleWorld;
        Pstream::gatherList(procToWorldIndex, Pstream::msgType(), comm_);
        Pstream::scatterList(procToWorldIndex, Pstream::msgType(), comm_);

        labelList nPerProc(nProcs);
        nPerProc[myRank] = patch_.size();
        Pstream::gatherList(nPerProc, Pstream::msgType(), comm_);
        Pstream::scatterList(nPerProc, Pstream::msgType(), comm_);

        patchFaceWorlds.setSize(patchFaces.size());
        patchFaceProcs.setSize(patchFaces.size());

        label sampleI = 0;
        forAll(nPerProc, proci)
        {
            for (label i = 0; i < nPerProc[proci]; i++)
            {
                patchFaceWorlds[sampleI] = procToWorldIndex[proci];
                patchFaceProcs[sampleI] = proci;
                sampleI++;
            }
        }
    }
    Pstream::warnComm = oldComm;
}


void Foam::mappedPatchBase::findLocalSamples
(
    const sampleMode mode,

    const label mySampleWorld,  // local world to sample == my own world
    const word& sampleRegion,   // local region to sample
    const word& samplePatch,    // local patch to sample

    const pointField& samples,
    List<nearInfoWorld>& nearest
) const
{
    // Find the local cell containing the samples

    const label myRank = Pstream::myProcNo(comm_);

    // Lookup the correct region
    const polyMesh& mesh = lookupMesh(sampleRegion);

    // All the info for nearest. Construct to miss
    nearest.setSize(samples.size());
    nearInfoWorld miss;
    {
        miss.first().second() = Tuple2<scalar, label>(Foam::sqr(GREAT), -1);
        miss.second() = -1; // set world to be ignored
    }
    nearest = miss;

    switch (mode)
    {
        case NEARESTCELL:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];
                nearInfoWorld& near = nearest[sampleI];

                label celli = tree.findInside(sample);

                if (celli == -1)
                {
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
                else
                {
                    const point& cc = mesh.cellCentres()[celli];

                    near.first().first() = pointIndexHit
                    (
                        true,
                        cc,
                        celli
                    );
                    near.first().second().first() = magSqr(cc-sample);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            break;
        }

        case NEARESTONLYCELL:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];
                nearInfoWorld& near = nearest[sampleI];

                near.first().first() = tree.findNearest(sample, sqr(GREAT));
                near.first().second().first() = magSqr
                (
                    near.first().first().hitPoint()
                   -sample
                );
                near.first().second().second() = myRank;
                near.second() = mySampleWorld;
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            Random rndGen(123456);

            const polyPatch& pp = lookupPatch(sampleRegion, samplePatch);

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearInfoWorld& near = nearest[sampleI];
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            else
            {
                // patch faces
                const labelList patchFaces(identity(pp.size(), pp.start()));

                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1e-4
                    )
                );

                patchBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
                patchBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

                indexedOctree<treeDataFace> boundaryTree
                (
                    treeDataFace    // all information needed to search faces
                    (
                        false,      // do not cache bb
                        mesh,
                        patchFaces  // boundary faces only
                    ),
                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    nearInfoWorld& near = nearest[sampleI];
                    pointIndexHit& nearInfo = near.first().first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        magSqr(patchBb.span())
                    );

                    if (!nearInfo.hit())
                    {
                        near.first().second().first() = Foam::sqr(GREAT);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                    else
                    {
                        point fc(pp[nearInfo.index()].centre(pp.points()));
                        nearInfo.setPoint(fc);
                        near.first().second().first() = magSqr(fc-sample);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                }
            }
            break;
        }

        case NEARESTPATCHPOINT:
        {
            Random rndGen(123456);

            const polyPatch& pp = lookupPatch(sampleRegion, samplePatch);

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearInfoWorld& near = nearest[sampleI];
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            else
            {
                // patch (local) points
                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1e-4
                    )
                );
                patchBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
                patchBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

                indexedOctree<treeDataPoint> boundaryTree
                (
                    treeDataPoint   // all information needed to search faces
                    (
                        mesh.points(),
                        pp.meshPoints() // selection of points to search on
                    ),
                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    nearInfoWorld& near = nearest[sampleI];
                    pointIndexHit& nearInfo = near.first().first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        magSqr(patchBb.span())
                    );

                    if (!nearInfo.hit())
                    {
                        near.first().second().first() = Foam::sqr(GREAT);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                    else
                    {
                        const point& pt = nearInfo.hitPoint();

                        near.first().second().first() = magSqr(pt-sample);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                }
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const meshSearchMeshObject& meshSearchEngine =
                meshSearchMeshObject::New(mesh);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];
                nearInfoWorld& near = nearest[sampleI];

                label facei = meshSearchEngine.findNearestFace(sample);

                if (facei == -1)
                {
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
                else
                {
                    const point& fc = mesh.faceCentres()[facei];

                    near.first().first() = pointIndexHit(true, fc, facei);
                    near.first().second().first() = magSqr(fc-sample);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            break;
        }

        case NEARESTPATCHFACEAMI:
        {
            // nothing to do here
            return;
        }

        default:
        {
            FatalErrorInFunction
                << "problem." << abort(FatalError);
        }
    }
}


// Find the processor/cell containing the samples. Does not account
// for samples being found in two processors.
void Foam::mappedPatchBase::findSamples
(
    const sampleMode mode,
    const label myWorld,
    const pointField& samples,
    const labelList& wantedWorlds,
    const labelList& origProcs,

    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    const label myRank = Pstream::myProcNo(comm_);
    const label nProcs = Pstream::nProcs(comm_);

    wordList samplePatches(nProcs);
    {
        const label oldComm(Pstream::warnComm);
        Pstream::warnComm = comm_;
        samplePatches[myRank] = samplePatch_;
        Pstream::gatherList(samplePatches, Pstream::msgType(), comm_);
        Pstream::scatterList(samplePatches, Pstream::msgType(), comm_);
        Pstream::warnComm = oldComm;
    }
    wordList sampleRegions(nProcs);
    {
        const label oldComm(Pstream::warnComm);
        Pstream::warnComm = comm_;
        sampleRegions[myRank] = sampleRegion_;
        Pstream::gatherList(sampleRegions, Pstream::msgType(), comm_);
        Pstream::scatterList(sampleRegions, Pstream::msgType(), comm_);
        Pstream::warnComm = oldComm;
    }

    // Find all the info for nearest
    List<nearInfoWorld> nearest(samples.size());
    forAll(nearest, samplei)
    {
        nearest[samplei].first() = nearInfo
        (
            pointIndexHit(),
            Tuple2<scalar, label>(Foam::sqr(GREAT), -1)
        );
        nearest[samplei].second() = wantedWorlds[samplei];
    }


    // Extract samples to search for locally
    {
        DynamicList<label> localMap(samples.size());
        forAll(wantedWorlds, samplei)
        {
            if (wantedWorlds[samplei] == myWorld)
            {
                localMap.append(samplei);
            }
        }
        DebugVar(localMap.size());

        if (localMap.size())
        {
            pointField localSamples(samples, localMap);
            labelList localOrigProcs(origProcs, localMap);

            //Assume single patch to sample for now
            const word localOrigPatch(samplePatches[localOrigProcs[0]]);
            const word localOrigRegion(sampleRegions[localOrigProcs[0]]);
            List<nearInfoWorld> localNearest(localSamples.size());
            Pout<< "*** Searching locally for " << localSamples.size()
                << " samples on region:" << localOrigRegion
                << " on patch:" << localOrigPatch << endl;
            findLocalSamples
            (
                mode,
                myWorld,
                localOrigRegion,
                localOrigPatch,
                localSamples,
                localNearest
            );
            UIndirectList<nearInfoWorld>(nearest, localMap) = localNearest;
        }
    }

    const label oldComm(Pstream::warnComm);
    Pstream::warnComm = comm_;

    // Find nearest. Combine on master.
    Pstream::listCombineGather
    (
        nearest,
        nearestWorldEqOp(),
        Pstream::msgType(),
        comm_
    );
    Pstream::listCombineScatter(nearest, Pstream::msgType(), comm_);

    if (debug)
    {
        Pout<< "** AFter combining:" << endl;
        forAll(nearest, samplei)
        {
            Pout<< "  sample:" << samples[samplei]
                << " origating from proc:" << origProcs[samplei]
                << " to be found on world:" << wantedWorlds[samplei] << nl
                << "    found on world:" << nearest[samplei].second() << nl
                << "    found on proc:"
                << nearest[samplei].first().second().second() << nl
                << "    found on patchfacei:"
                << nearest[samplei].first().first().index() << nl
                << "    found at location:"
                << nearest[samplei].first().first().rawPoint() << nl;
        }
        Pout<< endl;
    }

    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());
    sampleLocations.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        const nearInfo& ni = nearest[sampleI].first();

        if (!ni.first().hit())
        {
            sampleProcs[sampleI] = -1;
            sampleIndices[sampleI] = -1;
            sampleLocations[sampleI] = vector::max;
        }
        else
        {
            sampleProcs[sampleI] = ni.second().second();
            sampleIndices[sampleI] = ni.first().index();
            sampleLocations[sampleI] = ni.first().hitPoint();
        }
    }

    Pstream::warnComm = oldComm;
}


void Foam::mappedPatchBase::calcMapping() const
{
    static bool hasWarned = false;
    if (mapPtr_.valid())
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    // Get points on face (since cannot use face-centres - might be off
    // face-diagonal decomposed tets.
    tmp<pointField> patchPoints(facePoints(patch_));

    // Get offsetted points
    const pointField offsettedPoints(samplePoints(patchPoints()));

    // Do a sanity check - am I sampling my own patch?
    // This only makes sense for a non-zero offset.
    bool sampleMyself =
    (
        mode_ == NEARESTPATCHFACE
     && sampleWorld() == UPstream::myWorld()
     && sampleRegion() == patch_.boundaryMesh().mesh().name()
     && samplePatch() == patch_.name()
    );

    if (sampleMyself)
    {
        // Check offset
        vectorField d(offsettedPoints-patchPoints());
        bool coincident = (gAverage(mag(d)) <= ROOTVSMALL);

        if (sampleMyself && coincident)
        {
            WarningInFunction
                << "Invalid offset " << d << endl
                << "Offset is the vector added to the patch face centres to"
                << " find the patch face supplying the data." << endl
                << "Setting it to " << d
                << " on the same patch, on the same region"
                << " will find the faces themselves which does not make sense"
                << " for anything but testing." << endl
                << "patch_:" << patch_.name() << endl
                << "sampleRegion_:" << sampleRegion() << endl
                << "mode_:" << sampleModeNames_[mode_] << endl
                << "samplePatch_:" << samplePatch() << endl
                << "offsetMode_:" << offsetModeNames_[offsetMode_] << endl;
        }
    }

    // Collect per processor the world
    const label nProcs = Pstream::nProcs(comm_);

    DynamicList<word> allSampleWorlds(nProcs);
    for (const auto& sampleWorld : UPstream::worlds())
    {
        if (!allSampleWorlds.found(sampleWorld))
        {
            allSampleWorlds.append(sampleWorld);
        }
    }

    // Get local world and world-to-sample in index form
    const label myWorld = allSampleWorlds.find(UPstream::myWorld());
    const label mySampleWorld = allSampleWorlds.find(sampleWorld_);


    // Get global list of all samples and the processor and face they come from.
    pointField samples;         // coordinates
    labelList patchFaceWorlds;  // world to sample
    labelList patchFaceProcs;   // originating processor
    labelList patchFaces;       // originating face on processor patch
    pointField patchFc;         // originating face centre
    collectSamples
    (
        mySampleWorld,          // world I want to sample
        patchPoints,

        samples,
        patchFaceWorlds,
        patchFaceProcs,
        patchFaces,
        patchFc
    );

    if (debug)
    {
        forAll(samples, samplei)
        {
            Pout<< "    sample:" << samples[samplei]
                << " origating from proc:" << patchFaceProcs[samplei]
                << "  face:" << patchFaces[samplei]
                << " to be found on world:" << patchFaceWorlds[samplei] << nl;
        }
    }

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples
    (
        mode_,
        myWorld,        // my world (in index form)
        samples,
        patchFaceWorlds,
        patchFaceProcs,

        sampleProcs,
        sampleIndices,
        sampleLocations
    );

    // Check for samples that were not found. This will only happen for
    // NEARESTCELL since finds cell containing a location
    if (mode_ == NEARESTCELL)
    {
        label nNotFound = 0;
        forAll(sampleProcs, sampleI)
        {
            if (sampleProcs[sampleI] == -1)
            {
                nNotFound++;
            }
        }
        reduce(nNotFound, sumOp<label>(), Pstream::msgType(), comm_);

        if (nNotFound > 0)
        {
            if (!hasWarned)
            {
                WarningInFunction
                    << "Did not find " << nNotFound
                    << " out of " << sampleProcs.size() << " total samples."
                    << " Sampling these on owner cell centre instead." << endl
                    << "On patch " << patch_.name()
                    << " on region " << sampleRegion()
                    << " in mode " << sampleModeNames_[mode_] << endl
                    << "with offset mode " << offsetModeNames_[offsetMode_]
                    << ". Suppressing further warnings from " << type() << endl;

                hasWarned = true;
            }

            // Collect the samples that cannot be found
            DynamicList<label> subMap;
            forAll(sampleProcs, sampleI)
            {
                if (sampleProcs[sampleI] == -1)
                {
                    subMap.append(sampleI);
                }
            }

            // And re-search for pure nearest (should not fail)
            labelList subSampleProcs;
            labelList subSampleIndices;
            pointField subSampleLocations;
            findSamples
            (
                mode_,
                myWorld,        // my world (in index form)

                pointField(samples, subMap),
                UIndirectList<label>(patchFaceWorlds, subMap)(),
                UIndirectList<label>(patchFaceProcs, subMap)(),

                subSampleProcs,
                subSampleIndices,
                subSampleLocations
            );

            // Insert
            labelUIndList(sampleProcs, subMap) = subSampleProcs;
            labelUIndList(sampleIndices, subMap) = subSampleIndices;
            UIndirectList<point>(sampleLocations, subMap) = subSampleLocations;
        }
    }

    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell/face sample is in (so source when mapping)
    //   sampleIndices, sampleProcs.

    if (debug && Pstream::master(comm_))
    {
        forAll(samples, i)
        {
            Pout<< i << " need data in region "
                << patch_.boundaryMesh().mesh().name()
                << " for proc:" << patchFaceProcs[i]
                << " face:" << patchFaces[i]
                << " at:" << patchFc[i] << endl
                << "Found data in region " << sampleRegion()
                << " at proc:" << sampleProcs[i]
                << " face:" << sampleIndices[i]
                << " at:" << sampleLocations[i]
                << nl << endl;
        }

        OFstream str
        (
            patch_.boundaryMesh().mesh().time().path()
          / patch_.name()
          + "_mapped.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres/points to file " << str.name()
            << endl;

        label vertI = 0;

        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }

    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchFaceProcs, comm_));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.

    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();

    forAll(subMap, proci)
    {
        subMap[proci] = labelUIndList(sampleIndices, subMap[proci]);
        constructMap[proci] = labelUIndList(patchFaces, constructMap[proci]);

        if (debug)
        {
            Pout<< "To proc:" << proci << " sending values of cells/faces:"
                << subMap[proci] << endl;
            Pout<< "From proc:" << proci
                << " receiving values of patch faces:"
                << constructMap[proci] << endl;
        }
    }


    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        bitSet used(patch_.size());
        forAll(constructMap, proci)
        {
            const labelList& map = constructMap[proci];

            forAll(map, i)
            {
                label facei = map[i];

                if (used.test(facei))
                {
                    FatalErrorInFunction
                        << "On patch " << patch_.name()
                        << " patchface " << facei
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
                else
                {
                    used.set(facei);
                }
            }
        }
        forAll(used, facei)
        {
            if (!used.test(facei))
            {
                FatalErrorInFunction
                    << "On patch " << patch_.name()
                    << " patchface " << facei
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }
}


const Foam::autoPtr<Foam::searchableSurface>& Foam::mappedPatchBase::surfPtr()
const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


void Foam::mappedPatchBase::calcAMI() const
{
    if (AMIPtr_.valid())
    {
        FatalErrorInFunction
            << "AMI already calculated" << exit(FatalError);
    }

    AMIPtr_.clear();


    const polyPatch& nbr = samplePolyPatch();

    // Transform neighbour patch to local system
    pointField nbrPoints(samplePoints(nbr.localPoints()));

    primitivePatch nbrPatch0
    (
        SubList<face>
        (
            nbr.localFaces(),
            nbr.size()
        ),
        nbrPoints
    );


    if (debug)
    {
        OFstream os(patch_.name() + "_neighbourPatch-org.obj");
        meshTools::writeOBJ(os, samplePolyPatch().localFaces(), nbrPoints);

        OFstream osN(patch_.name() + "_neighbourPatch-trans.obj");
        meshTools::writeOBJ(osN, nbrPatch0, nbrPoints);

        OFstream osO(patch_.name() + "_ownerPatch.obj");
        meshTools::writeOBJ(osO, patch_.localFaces(), patch_.localPoints());
    }

    // Construct/apply AMI interpolation to determine addressing and weights
    AMIPtr_.reset
    (
        new AMIPatchToPatchInterpolation
        (
            patch_,
            nbrPatch0,
            surfPtr(),
            faceAreaIntersect::tmMesh,
            true,
            AMIPatchToPatchInterpolation::imFaceAreaWeight,
            -1,
            AMIReverse_
        )
    );
}


// Hack to read old (List-based) format. See Field.C. The difference
// is only that in case of falling back to old format it expects a non-uniform
// list instead of a single vector.
Foam::tmp<Foam::pointField> Foam::mappedPatchBase::readListOrField
(
    const word& keyword,
    const dictionary& dict,
    const label size
)
{
    tmp<pointField> tfld(new pointField());
    pointField& fld = tfld.ref();

    if (size)
    {
        ITstream& is = dict.lookup(keyword);

        // Read first token
        token firstToken(is);

        if (firstToken.isWord())
        {
            if (firstToken.wordToken() == "uniform")
            {
                fld.setSize(size);
                fld = pTraits<vector>(is);
            }
            else if (firstToken.wordToken() == "nonuniform")
            {
                is >> static_cast<List<vector>&>(fld);
                if (fld.size() != size)
                {
                    FatalIOErrorInFunction(dict)
                        << "size " << fld.size()
                        << " is not equal to the given value of " << size
                        << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Expected keyword 'uniform' or 'nonuniform', found "
                    << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else if (is.version() == IOstream::versionNumber(2,0))
        {
            IOWarningInFunction(dict)
                << "Expected keyword 'uniform' or 'nonuniform', "
                   "assuming List format for backwards compatibility."
                   "Foam version 2.0." << endl;

            is.putBack(firstToken);
            is >> static_cast<List<vector>&>(fld);
        }
    }
    return tfld;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::mappedPatchBase(const polyPatch& pp)
:
    patch_(pp),
    sampleWorld_(""),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    mode_(NEARESTPATCHFACE),
    samplePatch_(""),
    coupleGroup_(),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(pp.size(), offset_),
    distance_(0),
    comm_(communicator(UPstream::myWorld(), sampleWorld_)),
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vectorField& offsets
)
:
    patch_(pp),
    sampleWorld_(""),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(NONUNIFORM),
    offset_(Zero),
    offsets_(offsets),
    distance_(0),
    comm_(communicator(UPstream::myWorld(), sampleWorld_)),
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vector& offset
)
:
    patch_(pp),
    sampleWorld_(""),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(UNIFORM),
    offset_(offset),
    offsets_(0),
    distance_(0),
    comm_(communicator(UPstream::myWorld(), sampleWorld_)),
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const scalar distance
)
:
    patch_(pp),
    sampleWorld_(""),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(NORMAL),
    offset_(Zero),
    offsets_(0),
    distance_(distance),
    comm_(communicator(UPstream::myWorld(), sampleWorld_)),
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleWorld_(dict.lookupOrDefault<word>("sampleWorld", "")),
    sampleRegion_(dict.lookupOrDefault<word>("sampleRegion", "")),
    mode_(sampleModeNames_.get("sampleMode", dict)),
    samplePatch_(dict.lookupOrDefault<word>("samplePatch", "")),
    coupleGroup_(dict),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0.0),
    comm_(communicator(UPstream::myWorld(), sampleWorld_)),
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(dict.lookupOrDefault("flipNormals", false)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (!coupleGroup_.valid())
    {
        if (sampleWorld_.empty() && sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }

    if (offsetModeNames_.readIfPresent("offsetMode", dict, offsetMode_))
    {
        switch (offsetMode_)
        {
            case UNIFORM:
            {
                dict.readEntry("offset", offset_);
            }
            break;

            case NONUNIFORM:
            {
                //offsets_ = pointField(dict.lookup("offsets"));
                offsets_ = readListOrField("offsets", dict, patch_.size());
            }
            break;

            case NORMAL:
            {
                dict.readEntry("distance", distance_);
            }
            break;
        }
    }
    else if (dict.readIfPresent("offset", offset_))
    {
        offsetMode_ = UNIFORM;
    }
    else if (dict.found("offsets"))
    {
        offsetMode_ = NONUNIFORM;
        //offsets_ = pointField(dict.lookup("offsets"));
        offsets_ = readListOrField("offsets", dict, patch_.size());
    }
    else if (mode_ != NEARESTPATCHFACE && mode_ != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply the offsetMode as one of "
            << offsetModeNames_
            << exit(FatalIOError);
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const sampleMode mode,
    const dictionary& dict
)
:
    patch_(pp),
    sampleWorld_(dict.lookupOrDefault<word>("sampleWorld", "")),
    sampleRegion_(dict.lookupOrDefault<word>("sampleRegion", "")),
    mode_(mode),
    samplePatch_(dict.lookupOrDefault<word>("samplePatch", "")),
    coupleGroup_(dict), //dict.lookupOrDefault<word>("coupleGroup", "")),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0.0),
    comm_(communicator(UPstream::myWorld(), sampleWorld_)),
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(dict.lookupOrDefault("flipNormals", false)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (mode != NEARESTPATCHFACE && mode != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Construct from sampleMode and dictionary only applicable for "
            << " collocated patches in modes "
            << sampleModeNames_[NEARESTPATCHFACE] << ','
            << sampleModeNames_[NEARESTPATCHFACEAMI]
            << exit(FatalIOError);
    }


    if (!coupleGroup_.valid())
    {
        if (sampleWorld_.empty() && sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb
)
:
    patch_(pp),
    sampleWorld_(mpb.sampleWorld_),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_(mpb.offsets_),
    distance_(mpb.distance_),
    comm_(mpb.comm_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb,
    const labelUList& mapAddressing
)
:
    patch_(pp),
    sampleWorld_(mpb.sampleWorld_),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_
    (
        offsetMode_ == NONUNIFORM
      ? vectorField(mpb.offsets_, mapAddressing)
      : vectorField(0)
    ),
    distance_(mpb.distance_),
    comm_(mpb.comm_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::~mappedPatchBase()
{
    clearOut();
}


void Foam::mappedPatchBase::clearOut()
{
    mapPtr_.clear();
    AMIPtr_.clear();
    surfPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedPatchBase::lookupMesh
(
    const word& sampleRegion
) const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();
    return
    (
        sampleRegion.empty() || sampleRegion == thisMesh.name()
      ? thisMesh
      : thisMesh.time().lookupObject<polyMesh>(sampleRegion)
    );
}


const Foam::polyPatch& Foam::mappedPatchBase::lookupPatch
(
    const word& sampleRegion,
    const word& samplePatch
) const
{
    const polyMesh& nbrMesh = lookupMesh(sampleRegion);

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch
            << " in region " << sampleRegion_ << endl
            << exit(FatalError);
    }
    return nbrMesh.boundaryMesh()[patchi];
}


const Foam::polyMesh& Foam::mappedPatchBase::sampleMesh() const
{
    if (UPstream::myWorld() != sampleWorld_)
    {
        //FatalErrorInFunction
        WarningInFunction
            << "sampleWorld : " << sampleWorld_
            << " is not the current world : " << UPstream::myWorld()
            //<< exit(FatalError);
            << endl;
    }
    return lookupMesh(sampleRegion());
}


const Foam::polyPatch& Foam::mappedPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch()
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints
(
    const pointField& fc
) const
{
    tmp<pointField> tfld(new pointField(fc));
    pointField& fld = tfld.ref();

    switch (offsetMode_)
    {
        case UNIFORM:
        {
            fld += offset_;
            break;
        }
        case NONUNIFORM:
        {
            fld += offsets_;
            break;
        }
        case NORMAL:
        {
            // Get outwards pointing normal
            vectorField n(patch_.faceAreas());
            n /= mag(n);

            fld += distance_*n;
            break;
        }
    }

    return tfld;
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints() const
{
    return samplePoints(facePoints(patch_));
}


Foam::pointIndexHit Foam::mappedPatchBase::facePoint
(
    const polyMesh& mesh,
    const label facei,
    const polyMesh::cellDecomposition decompMode
)
{
    const point& fc = mesh.faceCentres()[facei];

    switch (decompMode)
    {
        case polyMesh::FACE_PLANES:
        case polyMesh::FACE_CENTRE_TRIS:
        {
            // For both decompositions the face centre is guaranteed to be
            // on the face
            return pointIndexHit(true, fc, facei);
        }
        break;

        case polyMesh::FACE_DIAG_TRIS:
        case polyMesh::CELL_TETS:
        {
            // Find the intersection of a ray from face centre to cell centre
            // Find intersection of (face-centre-decomposition) centre to
            // cell-centre with face-diagonal-decomposition triangles.

            const pointField& p = mesh.points();
            const face& f = mesh.faces()[facei];

            if (f.size() <= 3)
            {
                // Return centre of triangle.
                return pointIndexHit(true, fc, 0);
            }

            label celli = mesh.faceOwner()[facei];
            const point& cc = mesh.cellCentres()[celli];
            vector d = fc-cc;

            const label fp0 = mesh.tetBasePtIs()[facei];
            const point& basePoint = p[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = p[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = p[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);
                pointHit hitInfo = tri.intersection
                (
                    cc,
                    d,
                    intersection::HALF_RAY
                );

                if (hitInfo.hit() && hitInfo.distance() > 0)
                {
                    return pointIndexHit(true, hitInfo.hitPoint(), i-2);
                }

                fp = nextFp;
            }

            // Fall-back
            return pointIndexHit(false, fc, -1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "problem" << abort(FatalError);
            return pointIndexHit();
        }
    }
}


void Foam::mappedPatchBase::write(Ostream& os) const
{
    os.writeEntry("sampleMode", sampleModeNames_[mode_]);
    if (!sampleWorld_.empty())
    {
        os.writeEntry("sampleWorld", sampleWorld_);
    }
    if (!sampleRegion_.empty())
    {
        os.writeEntry("sampleRegion", sampleRegion_);
    }
    if (!samplePatch_.empty())
    {
        os.writeEntry("samplePatch", samplePatch_);
    }
//    os.writeEntryIfDifferent<word>("sampleWorld", sampleWorld_, "");
//    os.writeEntryIfDifferent<word>("sampleRegion", sampleRegion_, "");
//    os.writeEntryIfDifferent<word>("samplePatch", samplePatch_, "");
    coupleGroup_.write(os);

    if
    (
        offsetMode_ == UNIFORM
     && offset_ == vector::zero
     && (mode_ == NEARESTPATCHFACE || mode_ == NEARESTPATCHFACEAMI)
    )
    {
        // Collocated mode. No need to write offset data
    }
    else
    {
        os.writeEntry("offsetMode", offsetModeNames_[offsetMode_]);

        switch (offsetMode_)
        {
            case UNIFORM:
            {
                os.writeEntry("offset", offset_);
                break;
            }
            case NONUNIFORM:
            {
                offsets_.writeEntry("offsets", os);
                break;
            }
            case NORMAL:
            {
                os.writeEntry("distance", distance_);
                break;
            }
        }

        if (mode_ == NEARESTPATCHFACEAMI)
        {
            if (AMIReverse_)
            {
                os.writeEntry("flipNormals", AMIReverse_);
            }

            if (!surfDict_.empty())
            {
                surfDict_.writeEntry(surfDict_.dictName(), os);
            }
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "ensightMesh.H"
#include "fvMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "mapDistribute.H"

#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "ensightOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::ensightMesh::internalZone = -1;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Find matching ids based on whitelist, blacklist
//
// An empty whitelist accepts everything that is not blacklisted.
// A regex match is trumped by a literal match.
// 
// Eg,
//     input:  ( abc apple wall wall1 wall2 )
//     whitelist:  ( abc  def  "wall.*" )
//     blacklist:  ( "[ab].*"  wall )
//
//     result:  (abc wall1 wall2)
//
static labelList getSelected
(
    const UList<word>& input,
    const wordRes& whitelist,
    const wordRes& blacklist
)
{
    const label len = input.size();

    if (whitelist.empty() && blacklist.empty())
    {
        return identity(len);
    }

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        const auto& text = input[i];

        bool accept = false;

        if (whitelist.size())
        {
            const auto result = whitelist.matched(text);

            accept =
            (
                result == wordRe::LITERAL
              ? true
              : (result == wordRe::REGEX && !blacklist.match(text))
            );
        }
        else
        {
            accept = !blacklist.match(text);
        }

        if (accept)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightMesh::clear()
{
    cellZoneParts_.clear();
    faceZoneParts_.clear();
    boundaryParts_.clear();
}


void Foam::ensightMesh::renumber()
{
    label partNo = 0;

    for (const label zoneId : cellZoneParts_.sortedToc())
    {
        cellZoneParts_[zoneId].index() = partNo++;
    }

    for (const label patchId : boundaryParts_.sortedToc())
    {
        boundaryParts_[patchId].index() = partNo++;
    }

    for (const label zoneId : faceZoneParts_.sortedToc())
    {
        faceZoneParts_[zoneId].index() = partNo++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::ensightMesh
(
    const fvMesh& mesh,
    const ensightMesh::options& opts
)
:
    options_(new options(opts)),
    mesh_(mesh),
    needsUpdate_(true)
{
    if (!option().lazy())
    {
        correct();
    }
}


Foam::ensightMesh::ensightMesh(const fvMesh& mesh)
:
    ensightMesh(mesh, ensightMesh::options(IOstream::streamFormat::BINARY))
{}


Foam::ensightMesh::ensightMesh
(
    const fvMesh& mesh,
    const IOstream::streamFormat format
)
:
    ensightMesh(mesh, ensightMesh::options(format))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightMesh::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::ensightMesh::expire()
{
    clear();

    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


void Foam::ensightMesh::correct()
{
    clear();

    // Track which cells are in a zone or not
    bitSet cellSelection;

    // Boundary faces to be excluded from export
    bitSet excludeFace;


    // All or specified cellZones first
    const wordRes& czMatcher = option().cellZoneSelection();

    if
    (
        option().useCellZones()
     && (option().useInternalMesh() || !czMatcher.empty())
    )
    {
        const wordList zoneNames(mesh_.cellZones().names());

        const labelList zoneIds =
        (
            czMatcher.empty()
          ? identity(zoneNames.size())     // Use all
          : czMatcher.matching(zoneNames)  // Use selected names
        );

        for (const label zoneId : zoneIds)
        {
            const word& zoneName = zoneNames[zoneId];
            const cellZone& zn = mesh_.cellZones()[zoneId];

            if (returnReduce(!zn.empty(), orOp<bool>()))
            {
                cellSelection.set(zn);

                ensightCells& part = cellZoneParts_(zoneId);

                part.clear();
                part.identifier() = zoneId;
                part.rename(zoneName);

                part.classify(mesh_, zn);

                // Finalize
                part.reduce();
            }
        }
    }

    if (option().useInternalMesh())
    {
        // The internal mesh

        if (!cellZoneParts_.empty())
        {
            // The unzoned cells - flip selection from zoned to unzoned
            cellSelection.flip();

            if (returnReduce(cellSelection.any(), orOp<bool>()))
            {
                ensightCells& part = cellZoneParts_(internalZone);

                part.clear();
                part.identifier() = internalZone;
                part.rename("internalMesh");

                part.classify(mesh_, cellSelection);

                // Finalize
                part.reduce();
            }

            // Handled all cells
            cellSelection.clearStorage();
        }
        else if (czMatcher.empty())
        {
            // The entire mesh, but only if there were no zones.
            // Skip this is any zones were specified, regardless of their existence

            ensightCells& part = cellZoneParts_(internalZone);

            part.clear();
            part.identifier() = internalZone;
            part.rename("internalMesh");

            part.classify(mesh_);

            // Finalize
            part.reduce();

            // Handled all cells
            cellSelection.clearStorage();
        }
    }


    if (option().useBoundaryMesh())
    {
        // Patches are output. Check that they are synced.
        mesh_.boundaryMesh().checkParallelSync(true);

        wordList patchNames = mesh_.boundaryMesh().names();
        if (Pstream::parRun())
        {
            // Do not include processor patches in matching
            patchNames.resize(mesh_.boundaryMesh().nNonProcessor());
        }

        const labelList patchIds
        (
            getSelected
            (
                patchNames,
                option().patchSelection(),
                option().patchExclude()
            )
        );

        for (const label patchId : patchIds)
        {
            const word& patchName = patchNames[patchId];

            // Use fvPatch (not polyPatch) to automatically remove empty patches
            const fvPatch& p = mesh_.boundary()[patchId];

            ensightFaces& part = boundaryParts_(patchId);

            part.clear();
            part.identifier() = patchId;
            part.rename(patchName);

            if (p.size())
            {
                // Local face addressing (offset = 0),
                // - this is what we'll need later when writing fields
                part.classify(p.patch());
            }
            else
            {
                // The patch is empty (on this processor)
                // or the patch is 'empty' (as fvPatch type)
                part.clear();
            }

            // Finalize
            part.reduce();

            if (!part.total())
            {
                boundaryParts_.erase(patchId);
            }
        }
    }


    if (option().useFaceZones())
    {
        const wordRes& fzMatcher = option().faceZoneSelection();

        const wordList zoneNames(mesh_.faceZones().names());

        const labelList zoneIds =
        (
            fzMatcher.empty()
          ? identity(zoneNames.size())      // Use all
          : fzMatcher.matching(zoneNames)   // Use selected names
        );


        if (zoneIds.size())
        {
            excludeFace.resize(mesh_.nFaces());

            for (const polyPatch& pp : mesh_.boundaryMesh())
            {
                const auto* procPatch = isA<processorPolyPatch>(pp);

                if (isA<emptyPolyPatch>(pp))
                {
                    excludeFace.set(pp.range());
                }
                else if (procPatch && !procPatch->owner())
                {
                    // Exclude neighbour-side, retain owner-side only
                    excludeFace.set(pp.range());
                }
            }
        }


        for (const label zoneId : zoneIds)
        {
            const word& zoneName = zoneNames[zoneId];
            const faceZone& zn = mesh_.faceZones()[zoneId];

            ensightFaces& part = faceZoneParts_(zoneId);

            part.clear();
            part.identifier() = zoneId;
            part.rename(zoneName);

            if (zn.size())
            {
                part.classify
                (
                    mesh_.faces(),
                    zn,
                    zn.flipMap(),
                    excludeFace
                );
            }

            // Finalize
            part.reduce();

            if (!part.total())
            {
                faceZoneParts_.erase(zoneId);
            }
        }
    }

    renumber();

    needsUpdate_ = false;
}


void Foam::ensightMesh::write(ensightGeoFile& os) const
{
    //
    // Write cellZones / internalMesh
    //
    for (const label zoneId : cellZoneParts_.sortedToc())
    {
        const ensightCells& part = cellZoneParts_[zoneId];

        // Renumber the points/faces into unique points
        autoPtr<globalIndex> globalPointsPtr;
        labelList pointToGlobal;  // local point to unique global index
        labelList uniqueMeshPointLabels;  // unique global points

        const bool usesAllCells =
            returnReduce((part.size() == mesh_.nCells()), andOp<bool>());

        if (usesAllCells)
        {
            // All cells used, and thus all points

            globalPointsPtr =
                mesh_.globalData().mergePoints
                (
                    pointToGlobal,
                    uniqueMeshPointLabels
                );
        }
        else
        {
            // Map mesh point index to local (compact) point index
            Map<label> meshPointMap(part.meshPointMap(mesh_));

            labelList meshPoints(meshPointMap.sortedToc());

            globalPointsPtr =
                mesh_.globalData().mergePoints
                (
                    meshPoints,
                    meshPointMap,
                    pointToGlobal,
                    uniqueMeshPointLabels
                );

            meshPointMap.clear();

            // The mergePoints returns pointToGlobal assuming local addressing
            // (eg, patch localFaces).
            // Recast as original mesh points to new global points

            labelList oldToNew(mesh_.nPoints(), -1);

            forAll(meshPoints, i)
            {
                const label orig = meshPoints[i];
                const label glob = pointToGlobal[i];

                oldToNew[orig] = glob;
            }

            pointToGlobal.transfer(oldToNew);
        }

        ensightOutput::Detail::writeCoordinates
        (
            os,
            part.index(),
            part.name(),
            globalPointsPtr().size(),   // nPoints (global)
            UIndirectList<point>(mesh_.points(), uniqueMeshPointLabels),
            Pstream::parRun()           //!< Collective write?
        );

        writeCellConnectivity(os, part, pointToGlobal);
    }


    //
    // Write patches - sorted by index
    //
    for (const label patchId : boundaryParts_.sortedToc())
    {
        const ensightFaces& part = boundaryParts_[patchId];

        const polyPatch& pp = mesh_.boundaryMesh()[patchId];

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;  // local point to unique global index
        labelList uniqueMeshPointLabels;  // unique global points

        autoPtr<globalIndex> globalPointsPtr =
            mesh_.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );


        ensightOutput::Detail::writeCoordinates
        (
            os,
            part.index(),
            part.name(),
            globalPointsPtr().size(),  // nPoints (global)
            UIndirectList<point>(mesh_.points(), uniqueMeshPointLabels),
            Pstream::parRun() //!< Collective write?
        );

        // Renumber the patch faces,
        // from local patch indexing to unique global index
        faceList patchFaces(pp.localFaces());
        ListListOps::inplaceRenumber(pointToGlobal, patchFaces);

        ensightOutput::writeFaceConnectivity
        (
            os,
            part,
            patchFaces,
            Pstream::parRun()           //!< Collective write?
        );
    }


    //
    // Write requested faceZones - sorted by index
    //
    for (const label zoneId : faceZoneParts_.sortedToc())
    {
        const ensightFaces& part = faceZoneParts_[zoneId];

        // Use the properly sorted faceIds (ensightFaces) and do NOT use the
        // faceZone directly, otherwise the point-maps will not correspond.
        // - perform face-flipping later

        indirectPrimitivePatch pp
        (
            IndirectList<face>(mesh_.faces(), part.faceIds()),
            mesh_.points()
        );

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;  // local point to unique global index
        labelList uniqueMeshPointLabels;  // unique global points

        autoPtr<globalIndex> globalPointsPtr =
            mesh_.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        ensightOutput::Detail::writeCoordinates
        (
            os,
            part.index(),
            part.name(),
            globalPointsPtr().size(),  // nPoints (global)
            UIndirectList<point>(mesh_.points(), uniqueMeshPointLabels),
            Pstream::parRun() //!< Collective write?
        );

        // Renumber the faces belonging to the faceZone,
        // from local numbering to unique global index.

        faceList patchFaces(pp.localFaces());
        ListListOps::inplaceRenumber(pointToGlobal, patchFaces);

        // Also a good place to perform face flipping
        const boolList& flip = part.flipMap();

        forAll(patchFaces, facei)
        {
            face& f = patchFaces[facei];

            if (flip[facei])
            {
                f.flip();
            }
        }

        ensightOutput::writeFaceConnectivityPresorted
        (
            os,
            part,
            patchFaces,
            Pstream::parRun()           //!< Collective write?
        );
    }
}


// ************************************************************************* //

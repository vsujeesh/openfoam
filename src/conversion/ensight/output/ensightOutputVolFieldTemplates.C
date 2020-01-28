/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "ensightOutputVolField.H"
#include "ensightMesh.H"

#include "fvMesh.H"
#include "globalIndex.H"
#include "linear.H"
#include "volPointInterpolation.H"
#include "interpolation.H"
#include "processorFvPatch.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * *  Detail * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::Detail::writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh
)
{
    constexpr bool parallel = true;

    const fvMesh& mesh = ensMesh.mesh();

    const Map<ensightCells>& cellZoneParts = ensMesh.cellZoneParts();
    const Map<ensightFaces>& faceZoneParts = ensMesh.faceZoneParts();
    const Map<ensightFaces>& boundaryParts = ensMesh.boundaryParts();

    //
    // Write internalMesh and cellZones - sorted by index
    //
    for (const label zoneId : cellZoneParts.sortedToc())
    {
        const ensightCells& part = cellZoneParts[zoneId];

        Detail::writeField(os, vf, part, parallel);
    }

    //
    // Write patches - sorted by index
    //
    for (const label patchId : boundaryParts.sortedToc())
    {
        const ensightFaces& part = boundaryParts[patchId];

        writeField
        (
            os,
            vf.boundaryField()[patchId],
            part,
            parallel,
            true // localField
        );
    }

    //
    // Write requested faceZones - sorted by index
    //
    if (!faceZoneParts.empty())
    {
        // Interpolates cell values to faces
        // - only needed when exporting faceZones...
        GeometricField<Type, fvsPatchField, surfaceMesh> sf
        (
            Foam::linearInterpolate(vf)
        );

        // Flat boundary field
        // similar to volPointInterpolation::flatBoundaryField()

        Field<Type> flat(mesh.nBoundaryFaces(), Zero);

        const fvBoundaryMesh& bm = mesh.boundary();
        forAll(vf.boundaryField(), patchi)
        {
            const polyPatch& pp = bm[patchi].patch();
            const auto& bf = vf.boundaryField()[patchi];

            if (isA<processorFvPatch>(bm[patchi]))
            {
                // Use average value for processor faces
                // own cell value = patchInternalField
                // nei cell value = evaluated boundary values
                SubList<Type>
                (
                    flat,
                    bf.size(),
                    pp.offset()
                ) = (0.5 * (bf.patchInternalField() + bf));
            }
            else if (!isA<emptyFvPatch>(bm[patchi]))
            {
                SubList<Type>
                (
                    flat,
                    bf.size(),
                    pp.offset()
                ) = bf;
            }
        }

        for (const label zoneId : faceZoneParts.sortedToc())
        {
            const ensightFaces& part = faceZoneParts[zoneId];

            // Loop over face ids to store the needed field values
            // - internal faces use linear interpolation
            // - boundary faces use the corresponding patch value

            // Field (local size)
            Field<Type> values(part.size());
            auto valIter = values.begin();

            for (const label faceId : part.faceIds())
            {
                *valIter =
                (
                    mesh.isInternalFace(faceId)
                  ? sf[faceId]
                  : flat[faceId - mesh.nInternalFaces()]
                );

                ++valIter;
            }

            // The field is already in the proper element order
            // - just need its corresponding sub-fields
            Detail::writeFaceSubField(os, values, part, parallel);
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writePointField
(
    ensightFile& os,
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& ensMesh
)
{
    const char* coordKeyword = "coordinates";

    constexpr bool parallel = true;

    const fvMesh& mesh = ensMesh.mesh();

    const Map<ensightCells>& cellZoneParts = ensMesh.cellZoneParts();
    const Map<ensightFaces>& faceZoneParts = ensMesh.faceZoneParts();
    const Map<ensightFaces>& boundaryParts = ensMesh.boundaryParts();

    //
    // Write internalMesh and cellZones - sorted by index
    //
    for (const label zoneId : cellZoneParts.sortedToc())
    {
        const ensightCells& part = cellZoneParts[zoneId];

        const bool usesAllCells =
            returnReduce((part.size() == mesh.nCells()), andOp<bool>());

        // Renumber the points/faces into unique points
        autoPtr<globalIndex> globalPointsPtr;
        labelList pointToGlobal;  // local point to unique global index
        labelList uniqueMeshPointLabels;  // unique global points

        if (usesAllCells)
        {
            // All cells used, and thus all points

            mesh.globalData().mergePoints
            (
                pointToGlobal,
                uniqueMeshPointLabels
            );
        }
        else
        {
            // Map mesh point index to local (compact) point index
            Map<label> meshPointMap(part.meshPointMap(mesh));

            globalPointsPtr =
                mesh.globalData().mergePoints
                (
                    meshPointMap.sortedToc(),
                    meshPointMap,
                    pointToGlobal,
                    uniqueMeshPointLabels
                );
        }

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        Detail::writeFieldComponents
        (
            os,
            coordKeyword,
            UIndirectList<Type>(pf.internalField(), uniqueMeshPointLabels),
            parallel
        );
    }


    //
    // Write patches - sorted by index
    //
    for (const label patchId : boundaryParts.sortedToc())
    {
        const ensightFaces& part = boundaryParts[patchId];

        const fvPatch& p = mesh.boundary()[patchId];

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                p.patch().meshPoints(),
                p.patch().meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        Detail::writeFieldComponents
        (
            os,
            coordKeyword,
            UIndirectList<Type>(pf.internalField(), uniqueMeshPointLabels),
            parallel
        );
    }

    //
    // Write requested faceZones - sorted by index
    //
    for (const label zoneId : faceZoneParts.sortedToc())
    {
        const ensightFaces& part = faceZoneParts[zoneId];

        uindirectPrimitivePatch p
        (
            UIndirectList<face>
            (
                mesh.faces(),
                part.faceIds()
            ),
            mesh.points()
        );

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                p.meshPoints(),
                p.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        Detail::writeFieldComponents
        (
            os,
            coordKeyword,
            UIndirectList<Type>(pf.internalField(), uniqueMeshPointLabels),
            parallel
        );
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
    const bool nodeValues
)
{
    if (nodeValues)
    {
        tmp<GeometricField<Type, pointPatchField, pointMesh>> pfld
        (
            volPointInterpolation::New(vf.mesh()).interpolate(vf)
        );
        pfld.ref().checkOut();
        pfld.ref().rename(vf.name());

        return Detail::writePointField<Type>(os, pfld, ensMesh);
    }

    return Detail::writeVolField<Type>(os, vf, ensMesh);
}


// * * * * * * * * * * * * * * * *  Serial * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::Serial::writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightParts& list
)
{
    // serial only until the geometry is also in parallel
    constexpr bool parallel = false;

    for (const ensightPart& item : list)
    {
        const ensightCells* cptr = isA<ensightCells>(item);

        if (cptr)
        {
            const ensightCells& part = *cptr;

            ensightOutput::Detail::writeField
            (
                os,
                vf.internalField(),
                part,
                parallel
            );
            continue;
        }


        const ensightFaces* fptr = isA<ensightFaces>(item);

        if (fptr)
        {
            const ensightFaces& part = *fptr;

            const label patchi = part.identifier();

            if (patchi >= 0 && patchi < vf.boundaryField().size())
            {
                ensightOutput::Detail::writeField
                (
                    os,
                    vf.boundaryField()[patchi],
                    part,
                    parallel,
                    true // localField
                );
            }
            continue;
        }
    }

    return true;
}


// ************************************************************************* //

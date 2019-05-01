/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "dynamicMotionSolverFvMeshAMI.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cyclicAMIPolyPatch.H"
#include "polyTopoChange.H"
#include "MeshObject.H"
#include "lduMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMotionSolverFvMeshAMI, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMotionSolverFvMeshAMI,
        IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMotionSolverFvMeshAMI::dynamicMotionSolverFvMeshAMI
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    motionPtr_(motionSolver::New(*this))
{}


Foam::dynamicMotionSolverFvMeshAMI::dynamicMotionSolverFvMeshAMI
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    dynamicFvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    motionPtr_(motionSolver::New(*this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::motionSolver& Foam::dynamicMotionSolverFvMeshAMI::motion() const
{
    return *motionPtr_;
}


bool Foam::dynamicMotionSolverFvMeshAMI::update()
{
    // Mesh not moved/changed yet
    moving(false);
    topoChanging(false);

    pointField newPoints(motionPtr_->curPoints());

    polyBoundaryMesh& pbm = const_cast<polyBoundaryMesh&>(boundaryMesh());

    // Scan all patches and see if we want to apply a mesh topology  update
    bool changeRequired = false;
    for (polyPatch& pp : pbm)
    {
        DebugInfo
            << "pre- topology change: patch " << pp.name()
            << " size:" << returnReduce(pp.size(), sumOp<label>())
            << " mag(faceAreas):" << gSum(mag(pp.faceAreas())) << endl;

        //changeRequired = pp.changeTopology(newPoints) || changeRequired;
        changeRequired = pp.changeTopology() || changeRequired;
    }

    reduce(changeRequired, orOp<bool>());

    if (changeRequired)
    {
        polyTopoChange polyTopo(*this);

        // Set new point positions in polyTopo object
        polyTopo.movePoints(newPoints);

        // Accumulate the patch-based mesh changes
        // Note:
        // - updates the AMIs using the new points
        // - creates a topo change object that removes old added faces and
        //   adds the new faces
        for (polyPatch& pp : pbm)
        {
            pp.setTopology(polyTopo);
        }

        // Update geometry
        // Note
        // - changeMesh leads to polyMesh::resetPrimitives which will also
        //   trigger polyBoundaryMesh::updateMesh (init and update) and
        //   ::calcGeometry
        // - BUT: mesh still corresponds to original (non-extended mesh) so
        //   we want to bypass these calls...
        autoPtr<mapPolyMesh> map =
            polyTopo.changeMesh
            (
                *this,
                true       // We will be calling movePoints after this update
            );

        // Apply topology change - update fv geometry and map fields
        // - polyMesh::updateMesh -> (fires initUpdateMesh and updateMesh in AMI BCs) called before mapFields
        // - AMI addressing must be up-to-date - used by, e.g. FaceCellWave
        // - will trigger (again) polyBoundaryMesh::updateMesh (init and update)
        updateMesh(map());

        // Move points and update derived properties
        // Note: resets face areas based on raw point locations!!!polyBoundaryMesh::updateMesh (init and update)polyBoundaryMesh::updateMesh (init and update)polyBoundaryMesh::updateMesh (init and update)
        // Note: processorPolyPatches will trigger calculation of faceCentres
        // (and therefore cell volumes), so need to update faceAreas in
        // initMovePoints since proc patches will be evaluated later than
        // AMI patches
        if (map().hasMotionPoints())
        {
            movePoints(map().preMotionPoints());
        }

        surfaceScalarField& meshPhi = setPhi();
        surfaceScalarField::Boundary& meshPhiBf = meshPhi.boundaryFieldRef();
        for (polyPatch& pp : pbm)
        {
            if (isA<cyclicAMIPolyPatch>(pp))
            {
// Might be the nbr side that was moving!!!

                const cyclicAMIPolyPatch& cami =
                    dynamic_cast<const cyclicAMIPolyPatch&>(pp);

                if (cami.owner())
                {
                    scalarField& phip = meshPhiBf[cami.index()];
                    forAll(phip, facei)
                    {
                        label meshFacei = cami.start() + facei;
                        const face& f = faces()[meshFacei];
                        scalar geomArea = f.mag(points());
                        scalar scaledArea = mag(cami.faceAreas()[facei]);
                        phip[facei] *= scaledArea/geomArea;
                    }

                    scalarField srcMeshPhi(phip);
                    if (Pstream::parRun())
                    {
                        cami.AMI().srcMap().distribute(srcMeshPhi);
                    }

                    const labelListList& tgtToSrcAddr = cami.AMI().tgtAddress();
                    const cyclicAMIPolyPatch& nbr = cami.neighbPatch();
                    scalarField& nbrPhip = meshPhiBf[nbr.index()];

                    forAll(nbr, tgtFacei)
                    {
                        label srcFacei = tgtToSrcAddr[tgtFacei][0];
                        nbrPhip[tgtFacei] = -srcMeshPhi[srcFacei];
                    }

                    DebugInfo
                        << "patch:" << cami.name()
                        << " sum(area):" << gSum(mag(cami.faceAreas()))
                        << " min(mag(faceAreas):" << gMin(mag(cami.faceAreas()))
                        <<  " sum(meshPhi):" << gSum(phip) << nl
                        << "patch:" << nbr.name()
                        << " sum(area):" << gSum(mag(nbr.faceAreas()))
                        << " min(mag(faceAreas):" << gMin(mag(nbr.faceAreas()))
                        << " sum(meshPhi):" << gSum(nbrPhip)
                        << endl;
                }
            }
        }
    }
    else
    {
        fvMesh::movePoints(newPoints);
    }

    volVectorField* Uptr = getObjectPtr<volVectorField>("U");

    if (Uptr)
    {
        Uptr->correctBoundaryConditions();
    }

    return true;
}


// ************************************************************************* //

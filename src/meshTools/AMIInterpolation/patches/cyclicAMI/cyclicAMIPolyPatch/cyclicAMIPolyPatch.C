/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "cyclicAMIPolyPatch.H"
#include "transformField.H"
#include "SubField.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaIntersect.H"
#include "ops.H"
#include "polyTopoChange.H"
#include "OBJstream.H"
#include "vectorList.H"

#define DEBUG(msg){Pout<< "[" << __FILE__ << ":" << __LINE__ << "] " << name() << ": " << msg << endl;}
#define DEBUG2(msg){Pout<< "[" << __FILE__ << ":" << __LINE__ << "] " << name() << ": " << #msg << "=" << msg << endl;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::vector Foam::cyclicAMIPolyPatch::findFaceNormalMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const vectorField n((faceCentres - rotationCentre_) ^ rotationAxis_);

    const scalarField magRadSqr(magSqr(n));

    label facei = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&) : patch: " << name() << nl
            << "    rotFace  = " << facei << nl
            << "    point    = " << faceCentres[facei] << nl
            << "    distance = " << Foam::sqrt(magRadSqr[facei])
            << endl;
    }

    return n[facei];
}


void Foam::cyclicAMIPolyPatch::calcTransforms
(
    const primitivePatch& half0,
    const pointField& half0Ctrs,
    const vectorField& half0Areas,
    const pointField& half1Ctrs,
    const vectorField& half1Areas
)
{
    if (transform() != neighbPatch().transform())
    {
        FatalErrorInFunction
            << "Patch " << name()
            << " has transform type " << transformTypeNames[transform()]
            << ", neighbour patch " << neighbPatchName()
            << " has transform type "
            << neighbPatch().transformTypeNames[neighbPatch().transform()]
            << exit(FatalError);
    }


    // Calculate transformation tensors

    switch (transform())
    {
        case ROTATIONAL:
        {
            tensor revT = Zero;

            if (rotationAngleDefined_)
            {
                const tensor T(rotationAxis_*rotationAxis_);

                const tensor S
                (
                    0, -rotationAxis_.z(), rotationAxis_.y(),
                    rotationAxis_.z(), 0, -rotationAxis_.x(),
                    -rotationAxis_.y(), rotationAxis_.x(), 0
                );

                const tensor revTPos
                (
                    T
                  + cos(rotationAngle_)*(tensor::I - T)
                  + sin(rotationAngle_)*S
                );

                const tensor revTNeg
                (
                    T
                  + cos(-rotationAngle_)*(tensor::I - T)
                  + sin(-rotationAngle_)*S
                );

                // Check - assume correct angle when difference in face areas
                // is the smallest
                const vector transformedAreaPos = gSum(half1Areas & revTPos);
                const vector transformedAreaNeg = gSum(half1Areas & revTNeg);
                const vector area0 = gSum(half0Areas);
                const scalar magArea0 = mag(area0) + ROOTVSMALL;

                // Areas have opposite sign, so sum should be zero when correct
                // rotation applied
                const scalar errorPos = mag(transformedAreaPos + area0);
                const scalar errorNeg = mag(transformedAreaNeg + area0);

                const scalar normErrorPos = errorPos/magArea0;
                const scalar normErrorNeg = errorNeg/magArea0;

                if (errorPos > errorNeg && normErrorNeg < matchTolerance())
                {
                    revT = revTNeg;
                    rotationAngle_ *= -1;
                }
                else
                {
                    revT = revTPos;
                }

                const scalar areaError = min(normErrorPos, normErrorNeg);

                if (areaError > matchTolerance())
                {
                    WarningInFunction
                        << "Patch areas are not consistent within "
                        << 100*matchTolerance()
                        << " % indicating a possible error in the specified "
                        << "angle of rotation" << nl
                        << "    owner patch     : " << name() << nl
                        << "    neighbour patch : " << neighbPatch().name()
                        << nl
                        << "    angle           : "
                        << radToDeg(rotationAngle_) << " deg" << nl
                        << "    area error      : " << 100*areaError << " %"
                        << "    match tolerance : " <<  matchTolerance()
                        << endl;
                }

                if (debug)
                {
                    scalar theta = radToDeg(rotationAngle_);

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }
            else
            {
                point n0 = Zero;
                point n1 = Zero;
                if (half0Ctrs.size())
                {
                    n0 = findFaceNormalMaxRadius(half0Ctrs);
                }
                if (half1Ctrs.size())
                {
                    n1 = -findFaceNormalMaxRadius(half1Ctrs);
                }

                reduce(n0, maxMagSqrOp<point>());
                reduce(n1, maxMagSqrOp<point>());

                n0.normalise();
                n1.normalise();

                // Extended tensor from two local coordinate systems calculated
                // using normal and rotation axis
                const tensor E0
                (
                    rotationAxis_,
                    (n0 ^ rotationAxis_),
                    n0
                );
                const tensor E1
                (
                    rotationAxis_,
                    (-n1 ^ rotationAxis_),
                    -n1
                );
                revT = E1.T() & E0;

                if (debug)
                {
                    scalar theta = radToDeg(acos(-(n0 & n1)));

                    Pout<< "cyclicAMIPolyPatch::calcTransforms: patch:"
                        << name()
                        << " Specified rotation:"
                        << " n0:" << n0 << " n1:" << n1
                        << " swept angle: " << theta << " [deg]"
                        << " reverse transform: " << revT
                        << endl;
                }
            }

            const_cast<tensorField&>(forwardT()) = tensorField(1, revT.T());
            const_cast<tensorField&>(reverseT()) = tensorField(1, revT);
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        case TRANSLATIONAL:
        {
            if (debug)
            {
                Pout<< "cyclicAMIPolyPatch::calcTransforms : patch:" << name()
                    << " Specified translation : " << separationVector_
                    << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()) = vectorField
            (
                1,
                separationVector_
            );
            const_cast<boolList&>(collocated()) = boolList(1, false);

            break;
        }
        default:
        {
            if (debug)
            {
                Pout<< "patch:" << name()
                    << " Assuming cyclic AMI pairs are colocated" << endl;
            }

            const_cast<tensorField&>(forwardT()).clear();
            const_cast<tensorField&>(reverseT()).clear();
            const_cast<vectorField&>(separation()).setSize(0);
            const_cast<boolList&>(collocated()) = boolList(1, true);

            break;
        }
    }

    if (debug)
    {
        Pout<< "patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


void Foam::cyclicAMIPolyPatch::restoreScaledGeometry()
{
    DebugInFunction << endl;

    if (boundaryMesh().mesh().hasCellVolumes())
    {
        WarningInFunction
            << "Mesh already has volumes set!"
            << endl;
    }

    const polyPatch& nbr = neighbPatch();
    vectorField::subField srcFaceAreas = faceAreas();
    vectorField::subField tgtFaceAreas = nbr.faceAreas();
    vectorField::subField srcFaceCentres = faceCentres();
    vectorField::subField tgtFaceCentres = nbr.faceCentres();

    if (debug)
    {
        Info<< "before: sum(mag(srcFaceAreas)):"
            << gSum(mag(srcFaceAreas)) << endl;
        Info<< "before: sum(mag(faceAreas0)):"
            << gSum(mag(faceAreas0_)) << endl;
        Info<< "before: sum(mag(tgtFaceAreas)):"
            << gSum(mag(tgtFaceAreas)) << endl;
        Info<< "before: sum(mag(nbrFaceAreas0)):"
            << gSum(mag(nbrFaceAreas0_)) << endl;
    }

    srcFaceAreas = faceAreas0_;
    srcFaceCentres = faceCentres0_;
    faceAreas0_.clear();
    faceCentres0_.clear();

    tgtFaceAreas = nbrFaceAreas0_;
    tgtFaceCentres = nbrFaceCentres0_;
    nbrFaceAreas0_.clear();
    nbrFaceCentres0_.clear();

    if (debug)
    {
        Info<< "after: sum(mag(srcFaceAreas)):"
            << gSum(mag(srcFaceAreas)) << endl;
        Info<< "after: sum(mag(faceAreas0)):"
            << gSum(mag(faceAreas0_)) << endl;
        Info<< "after: sum(mag(tgtFaceAreas)):"
            << gSum(mag(tgtFaceAreas)) << endl;
        Info<< "after: sum(mag(nbrFaceAreas0)):"
            << gSum(mag(nbrFaceAreas0_)) << endl;
    }
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

bool Foam::cyclicAMIPolyPatch::removeAMIFaces(polyTopoChange& topoChange)
{
    DebugInFunction << endl;

    if (!owner())
    {
        return false;
    }

    bool changeRequired = false;

    // Remove any faces that we inserted...

    const cyclicAMIPolyPatch& nbr = neighbPatch();

    const label newSrcFaceStart = srcFaceIDs_.size();

    if (newSrcFaceStart != 0)
    {
        for (label facei = newSrcFaceStart; facei < size(); ++facei)
        {
            changeRequired = true;
            label meshFacei = start() + facei;
            topoChange.removeFace(meshFacei, -1);
        }
    }

    const label newTgtFaceStart = tgtFaceIDs_.size();

    if (newTgtFaceStart != 0)
    {
        for (label facei = newTgtFaceStart; facei < nbr.size(); ++facei)
        {
            changeRequired = true;
            label meshFacei = nbr.start() + facei;
            topoChange.removeFace(meshFacei, -1);
        }
    }

    srcFaceIDs_.clear();
    tgtFaceIDs_.clear();

    return changeRequired;
}


bool Foam::cyclicAMIPolyPatch::addAMIFaces(polyTopoChange& topoChange)
{
    DebugInFunction << endl;

    bool changedFaces = false;
    const cyclicAMIPolyPatch& nbr = neighbPatch();

    polyMesh& mesh = const_cast<polyMesh&>(boundaryMesh().mesh());
    const faceZoneMesh& faceZones = mesh.faceZones();

    // First face address and weight are used to manipulate the
    // original face - all other addresses and weights are used to
    // create additional faces
    const labelListList& srcToTgtAddr = AMI().srcAddress();
    const labelListList& tgtToSrcAddr = AMI().tgtAddress();
    srcFaceIDs_.setSize(srcToTgtAddr.size());
    tgtFaceIDs_.setSize(tgtToSrcAddr.size());

    label nNewSrcFaces = 0;
    forAll(srcToTgtAddr, srcFacei)
    {
        const labelList& tgtAddr = srcToTgtAddr[srcFacei];

        // No tgt faces linked to srcFacei
        if (tgtAddr.empty()) continue;

        srcFaceIDs_[srcFacei].setSize(tgtAddr.size());
        srcFaceIDs_[srcFacei][0] = srcFacei;

        label meshFacei = start() + srcFacei;
        for (label addri = 1; addri < tgtAddr.size(); ++addri)
        {
            changedFaces = true;

            // Note: new faces reuse originating face points
            // - but areas are scaled by the weights (later)

            // New source face for each target face address
            srcFaceIDs_[srcFacei][addri] = nNewSrcFaces + srcToTgtAddr.size();
            ++nNewSrcFaces;
            (void)topoChange.addFace
            (
                mesh.faces()[meshFacei],        // modified face
                mesh.faceOwner()[meshFacei],    // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                meshFacei,                      // master face
                false,                          // face flip
                index(),                        // patch for face
                faceZones.whichZone(meshFacei), // zone for original face
                false                           // face flip in zone
            );
        }
    }

    label nNewTgtFaces = 0;
    forAll(tgtToSrcAddr, tgtFacei)
    {
        const labelList& srcAddr = tgtToSrcAddr[tgtFacei];

        // No src faces linked to tgtFacei
        if (srcAddr.empty()) continue;

        tgtFaceIDs_[tgtFacei].setSize(srcAddr.size());
        tgtFaceIDs_[tgtFacei][0] = tgtFacei;

        label meshFacei = nbr.start() + tgtFacei;
        for (label addri = 1; addri < srcAddr.size(); ++addri)
        {
            changedFaces = true;

            // Note: new faces reuse originating face points
            // - but areas are scaled by the weights (later)

            // New target face for each source face address
            tgtFaceIDs_[tgtFacei][addri] = nNewTgtFaces + tgtToSrcAddr.size();
            ++nNewTgtFaces;

            (void)topoChange.addFace
            (
                mesh.faces()[meshFacei],        // modified face
                mesh.faceOwner()[meshFacei],    // owner
                -1,                             // neighbour
                -1,                             // master point
                -1,                             // master edge
                meshFacei,                      // master face
                false,                          // face flip
                nbr.index(),                    // patch for face
                faceZones.whichZone(meshFacei), // zone for original face
                false                           // face flip in zone
            );
        }
    }

    Info<< "New faces - " << name() << ": "
        << returnReduce(nNewSrcFaces, sumOp<label>())
        << " "  << nbr.name() << ": "
        << returnReduce(nNewTgtFaces, sumOp<label>())
        << endl;

    if (debug)
    {
        Pout<< "New faces - " << name() << ": " << nNewSrcFaces
            << " "  << nbr.name() << ": " << nNewTgtFaces << endl;
    }

    return returnReduce(changedFaces, orOp<bool>());
}


void Foam::cyclicAMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    resetAMI(boundaryMesh().mesh().points(), AMIMethod);
}


void Foam::cyclicAMIPolyPatch::resetAMI
(
    const UList<point>& points,
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    DebugInFunction << endl;

    if (!owner())
    {
        return;
    }

    AMIPtr_.clear();

    const cyclicAMIPolyPatch& nbr = neighbPatch();
    pointField srcPoints(points, meshPoints());
    pointField nbrPoints(points, neighbPatch().meshPoints());

    if (debug)
    {
        const Time& t = boundaryMesh().mesh().time();
        OFstream os(t.path()/name() + "_neighbourPatch-org.obj");
        meshTools::writeOBJ(os, neighbPatch().localFaces(), nbrPoints);
    }

    label patchSize0 = size();
    label nbrPatchSize0 = nbr.size();

    if (createAMIFaces_)
    {
        // AMI is created based on the original patch faces (non-extended patch)
        if (srcFaceIDs_.size())
        {
            patchSize0 = srcFaceIDs_.size();
        }
        if (tgtFaceIDs_.size())
        {
            nbrPatchSize0 = tgtFaceIDs_.size();
        }
    }

    // Transform neighbour patch to local system
    transformPosition(nbrPoints);
    primitivePatch nbrPatch0
    (
        SubList<face>
        (
            nbr.localFaces(),
            nbrPatchSize0
        ),
        nbrPoints
    );
    primitivePatch patch0
    (
        SubList<face>
        (
            localFaces(),
            patchSize0
        ),
        srcPoints
    );


    if (debug)
    {
        const Time& t = boundaryMesh().mesh().time();
        OFstream osN(t.path()/name() + "_neighbourPatch-trans.obj");
        meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

        OFstream osO(t.path()/name() + "_ownerPatch.obj");
        meshTools::writeOBJ(osO, this->localFaces(), localPoints());
    }

    // Construct/apply AMI interpolation to determine addressing and weights
    AMIPtr_.reset
    (
        new AMIPatchToPatchInterpolation
        (
            patch0, // *this,
            nbrPatch0,
            surfPtr(),
            faceAreaIntersect::tmMesh,
            AMIRequireMatch_,
            AMIMethod,
            AMILowWeightCorrection_,
            AMIReverse_
        )
    );

    // Set the updated flag
    updated_ = true;

    if (debug)
    {
        AMIPtr_->checkSymmetricWeights(true);
    }
}


void Foam::cyclicAMIPolyPatch::setAMIFaces()
{
    DebugInFunction << endl;

    if (!owner())
    {
        return;
    }

    // Create new mesh faces so that there is a 1-to-1 correspondence
    // between faces on each side of the AMI
    const cyclicAMIPolyPatch& nbr = neighbPatch();

    vectorField::subField srcFaceAreas = faceAreas();
    vectorField::subField tgtFaceAreas = nbr.faceAreas();

    // Scale the new face areas and set the centroids
    // Note:
    // - storing local copies so that they can be re-applied after the call to
    //   movePoints that will reset any changes to the areas and centroids
    //
    // - For AMI, src and tgt patches should be the same
    // - For ACMI they are likely to be different!
    faceAreas0_ = srcFaceAreas;
    faceCentres0_.setSize(size());
    nbrFaceAreas0_ = tgtFaceAreas;
    nbrFaceCentres0_.setSize(nbr.size());

    const labelListList& srcToTgtAddr0 = AMIPtr_->srcAddress();
    const labelListList& tgtToSrcAddr0 = AMIPtr_->tgtAddress();
    const pointListList& srcCtr0 = AMIPtr_->srcCentroids();
    const scalarListList& srcToTgtWght0 = AMIPtr_->srcWeights();

    // New addressing on extended mesh
    labelListList srcToTgtAddr1(size(), labelList());
    labelListList tgtToSrcAddr1(nbr.size(), labelList());

    // New maps (mesh has changed since AMI was computed)
    autoPtr<mapDistribute> srcToTgtMap1;
    autoPtr<mapDistribute> tgtToSrcMap1;

    if (AMIPtr_->singlePatchProc() == -1)
    {
        // Parallel running

        // Global index based on old patch sizes (when AMI was computed)
        globalIndex globalSrcFaces0(srcToTgtAddr0.size());
        globalIndex globalTgtFaces0(tgtToSrcAddr0.size());

        // Global index based on new patch sizes
        globalIndex globalSrcFaces1(size());
        globalIndex globalTgtFaces1(nbr.size());


        // Gather source side info
        // =======================

        // Note: using new global index for addressing, and distributed using
        // the old AMI map
        labelListList newTgtGlobalFaces(tgtFaceIDs_);
        forAll(newTgtGlobalFaces, tgtFacei)
        {
            globalTgtFaces1.inplaceToGlobal(newTgtGlobalFaces[tgtFacei]);
        }
        AMIPtr_->tgtMap().distribute(newTgtGlobalFaces);

        // Now have new tgt face indices for each src face

        labelList globalSrcFaceIDs(identity(srcToTgtAddr0.size()));
        globalSrcFaces0.inplaceToGlobal(globalSrcFaceIDs);
        AMIPtr_->srcMap().distribute(globalSrcFaceIDs);
        // globalSrcFaceIDs now has remote data for each srcFacei0 known to the
        // tgt patch

        List<List<point>> globalSrcCtrs0(srcCtr0);
        AMIPtr_->srcMap().distribute(globalSrcCtrs0);

        labelList globalTgtFaceIDs(identity(tgtToSrcAddr0.size()));
        globalTgtFaces0.inplaceToGlobal(globalTgtFaceIDs);
        AMIPtr_->tgtMap().distribute(globalTgtFaceIDs);
        // globalTgtFaceIDs now has remote data for each tgtFacei0 known to the
        // src patch

        // For debug - send tgt face centres and compare against mapped src
        // face centres
        //List<List<point>> globalTgtCtrs0(tgtCtr0);
        //AMIPtr_->tgtMap().distribute(globalTgtCtrs0);

        labelListList globalTgtToSrcAddr(tgtToSrcAddr0);
        forAll(tgtToSrcAddr0, tgtFacei0)
        {
            forAll(tgtToSrcAddr0[tgtFacei0], addri)
            {
                const label globalSrcFacei =
                    globalSrcFaceIDs[tgtToSrcAddr0[tgtFacei0][addri]];
                globalTgtToSrcAddr[tgtFacei0][addri] = globalSrcFacei;
            }
        }
        AMIPtr_->tgtMap().distribute(globalTgtToSrcAddr);

        labelListList globalSrcToTgtAddr(srcToTgtAddr0);
        forAll(srcToTgtAddr0, srcFacei0)
        {
            forAll(srcToTgtAddr0[srcFacei0], addri)
            {
                const label globalTgtFacei =
                    globalTgtFaceIDs[srcToTgtAddr0[srcFacei0][addri]];
                globalSrcToTgtAddr[srcFacei0][addri] = globalTgtFacei;
            }
        }
        AMIPtr_->srcMap().distribute(globalSrcToTgtAddr);

        label nError = 0;
        forAll(srcToTgtAddr0, srcFacei0)
        {
            const labelList& newSrcFaces = srcFaceIDs_[srcFacei0];
            forAll(newSrcFaces, i)
            {
                label srcFacei1 = newSrcFaces[i];

                // What index did srcFacei0 appear in tgtToSrc0 list?
                // - if first index, all ok
                // - else tgt face has been moved to according to tgtFaceIDs_
                label tgtFacei0 = srcToTgtAddr0[srcFacei0][i];
                label addri =
                    globalTgtToSrcAddr[tgtFacei0].find
                    (
                        globalSrcFaceIDs[srcFacei0]
                    );

                if (addri == -1)
                {
                    ++nError;
                    continue;

                    if (debug)
                    {
                        Pout<< "Unable to find global source face "
                            << globalSrcFaceIDs[srcFacei0]
                            << " in globalTgtToSrcAddr[" << tgtFacei0 << "]: "
                            << globalTgtToSrcAddr[tgtFacei0]
                            << endl;
                    }
                }

                label tgtFacei1 = newTgtGlobalFaces[tgtFacei0][addri];

                // Sanity check to see that we've picked the correct face
                // point tgtCtr0(globalTgtCtrs0[tgtFacei0][addri]);
                // Pout<< "srcCtr:" << srcCtr0[srcFacei0][i]
                //     << " tgtCtr:" << tgtCtr0 << endl;

                srcToTgtAddr1[srcFacei1] = labelList(1, tgtFacei1);
                faceAreas0_[srcFacei1] *= srcToTgtWght0[srcFacei0][i];
                faceCentres0_[srcFacei1] = srcCtr0[srcFacei0][i];
            }
        }

        if (nError)
        {
            FatalErrorInFunction
                << "Unable to find " << nError << " global source faces"
                << abort(FatalError);
        }


        // Gather Target side info
        // =======================

        labelListList newSrcGlobalFaces(srcFaceIDs_);
        forAll(newSrcGlobalFaces, srcFacei)
        {
            globalSrcFaces1.inplaceToGlobal(newSrcGlobalFaces[srcFacei]);
        }

        AMIPtr_->srcMap().distribute(newSrcGlobalFaces);

        // Now have new src face indices for each tgt face
        forAll(tgtToSrcAddr0, tgtFacei0)
        {
            const labelList& newTgtFaces = tgtFaceIDs_[tgtFacei0];
            forAll(newTgtFaces, i)
            {
                label srcFacei0 = tgtToSrcAddr0[tgtFacei0][i];

                label addri =
                    globalSrcToTgtAddr[srcFacei0].find
                    (
                        globalTgtFaceIDs[tgtFacei0]
                    );

                if (addri == -1)
                {
                    ++nError;
                    continue;

                    if (debug)
                    {
                        Pout<< "Unable to find global target face "
                            << globalTgtFaceIDs[tgtFacei0]
                            << " in globalSrcToTgtAddr[" << srcFacei0 << "]: "
                            << globalSrcToTgtAddr[srcFacei0]
                            << endl;
                    }
                }

                label srcFacei1 = newSrcGlobalFaces[srcFacei0][addri];

                // Sanity check to see that we've picked the correct face
                point srcCtr0(globalSrcCtrs0[srcFacei0][addri]);
                reverseTransformPosition(srcCtr0, srcFacei0);

                label tgtFacei1 = newTgtFaces[i];
                tgtToSrcAddr1[tgtFacei1] = labelList(1, srcFacei1);
                nbrFaceCentres0_[tgtFacei1] = srcCtr0;
            }
        }

        if (nError)
        {
            FatalErrorInFunction
                << "Unable to find " << nError << " global target faces"
                << abort(FatalError);
        }

        // Update the maps
        {
            List<Map<label>> cMap;
            srcToTgtMap1.reset
            (
                new mapDistribute(globalSrcFaces1, tgtToSrcAddr1, cMap)
            );
        }
        {
            List<Map<label>> cMap;
            tgtToSrcMap1.reset
            (
                new mapDistribute(globalTgtFaces1, srcToTgtAddr1, cMap)
            );
        }

        // Reset tgt patch areas using the new map
        vectorList newSrcGlobalFaceAreas(faceAreas0_);
        srcToTgtMap1->distribute(newSrcGlobalFaceAreas);
        forAll(tgtFaceAreas, tgtFacei)
        {
            label srcFacei = tgtToSrcAddr1[tgtFacei][0];
            nbrFaceAreas0_[tgtFacei] = -newSrcGlobalFaceAreas[srcFacei];
        }
    }
    else
    {
        label nError = 0;
        forAll(srcToTgtAddr0, srcFacei0)
        {
            const labelList& srcFaceTgtAddr = srcToTgtAddr0[srcFacei0];
            const scalarList& srcFaceTgtWght = srcToTgtWght0[srcFacei0];
            const pointList& srcFaceTgtCtr = srcCtr0[srcFacei0];
            forAll(srcFaceTgtAddr, addri)
            {
                label newSrcFacei = srcFaceIDs_[srcFacei0][addri];

                // Find which slot srcFacei0 appears in tgt->src addressing
                label oldTgtFacei = srcFaceTgtAddr[addri];
                label tgtAddri = tgtToSrcAddr0[oldTgtFacei].find(srcFacei0);

                if (tgtAddri == -1)
                {
                    ++nError;
                    continue;

                    if (debug)
                    {
                        Pout<< "Unable to find source face " << srcFacei0
                            << " in tgtToSrcAddr0[" << oldTgtFacei << "]: "
                            << tgtToSrcAddr0[oldTgtFacei]
                            << endl;
                    }
                }

                label newTgtFacei = tgtFaceIDs_[oldTgtFacei][tgtAddri];

                faceAreas0_[newSrcFacei] *= srcFaceTgtWght[addri];
                nbrFaceAreas0_[newTgtFacei] = -faceAreas0_[newSrcFacei];

                point pt(srcFaceTgtCtr[addri]);
                faceCentres0_[newSrcFacei] = pt;
                reverseTransformPosition(pt, srcFacei0);
                nbrFaceCentres0_[newTgtFacei] = pt;

                // SANITY CHECK
                // Info<< "srcPt:" << srcFaceCentres[newSrcFacei]
                //     << " tgtPt:" << tgtFaceCentres[newTgtFacei] << endl;

                srcToTgtAddr1[newSrcFacei] = labelList(1, newTgtFacei);
                tgtToSrcAddr1[newTgtFacei] = labelList(1, newSrcFacei);
            }
        }

        if (nError)
        {
            FatalErrorInFunction
                << "Unable to find " << nError
                << " source faces in tgtToSrcAddr0"
                << abort(FatalError);
        }
    }

    // Update the AMI addressing and weights to reflect the new 1-to-1
    // correspondence
    AMIPtr_->update
    (
        std::move(srcToTgtMap1),
        std::move(tgtToSrcMap1),
        std::move(srcToTgtAddr1),
        scalarListList(srcToTgtAddr1.size(), scalarList(1, scalar(1))),
        std::move(tgtToSrcAddr1),
        scalarListList(tgtToSrcAddr1.size(), scalarList(1, scalar(1)))
    );

    AMIPtr_->setAreas(mag(faceAreas0_), mag(nbrFaceAreas0_));

    if (debug)
    {
        Pout<< "cyclicAMIPolyPatch : " << name()
            << " constructed AMI with " << nl
            << "    " << "srcAddress:" << AMIPtr_().srcAddress().size()
            << nl
            << "    " << "tgAddress :" << AMIPtr_().tgtAddress().size()
            << nl << endl;
    }
}


void Foam::cyclicAMIPolyPatch::calcTransforms()
{
    DebugInFunction << endl;

    const cyclicAMIPolyPatch& half0 = *this;
    vectorField half0Areas(half0.size());
    forAll(half0, facei)
    {
        half0Areas[facei] = half0[facei].areaNormal(half0.points());
    }

    const cyclicAMIPolyPatch& half1 = neighbPatch();
    vectorField half1Areas(half1.size());
    forAll(half1, facei)
    {
        half1Areas[facei] = half1[facei].areaNormal(half1.points());
    }

    calcTransforms
    (
        half0,
        half0.faceCentres(),
        half0Areas,
        half1.faceCentres(),
        half1Areas
    );

    if (debug)
    {
        Pout<< "calcTransforms() : patch: " << name() << nl
            << "    forwardT = " << forwardT() << nl
            << "    reverseT = " << reverseT() << nl
            << "    separation = " << separation() << nl
            << "    collocated = " << collocated() << nl << endl;
    }
}


void Foam::cyclicAMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    if (updatingAMI_)
    {
        resetAMI(AMIMethod_);
    }

    polyPatch::initGeometry(pBufs);


    // Early calculation of transforms so e.g. cyclicACMI can use them.
    // Note: also triggers primitiveMesh face centre. Note that cell
    // centres should -not- be calculated
    // since e.g. cyclicACMI override face areas
    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;
}


bool Foam::cyclicAMIPolyPatch::changeTopology() const
{
    DebugInFunction << endl;

    updatingAMI_ = true;

    createAMIFaces_ = true;

    return true;
}


bool Foam::cyclicAMIPolyPatch::setTopology(polyTopoChange& topoChange)
{
    DebugInFunction << endl;

    if (createAMIFaces_ && updatingAMI_ && owner())
    {
        resetAMI(topoChange.points(), AMIMethod_);

        removeAMIFaces(topoChange);

        addAMIFaces(topoChange);

        return true;
    }

    return false;
}


void Foam::cyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    DebugInFunction << endl;

    // See below. Clear out any local geometry
    primitivePatch::movePoints(p);

    // Note: processorPolyPatch::initMovePoints calls
    // processorPolyPatch::initGeometry which will trigger calculation of
    // patch faceCentres() and cell volumes...


    if (owner())
    {
        if (createAMIFaces_)
        {
            // faceAreas() and faceCentres() have been reset and will be
            // recalculated on-demand using the mesh points and no longer
            // correspond to the scaled areas!
            restoreScaledGeometry();
        }
        else
        {
            resetAMI(p, AMIMethod_);
        }

        // deltas need to be recalculated to use new face centres!
    }

    // Early calculation of transforms. See above.
    calcTransforms();
}


void Foam::cyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    DebugInFunction << endl;

    polyPatch::movePoints(pBufs, p);
/*
    polyPatch::movePoints -> primitivePatch::movePoints -> primitivePatch::clearGeom:
    deleteDemandDrivenData(localPointsPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
    deleteDemandDrivenData(magFaceAreasPtr_);
    deleteDemandDrivenData(faceNormalsPtr_);
    deleteDemandDrivenData(pointNormalsPtr_);
*/
    updatingAMI_ = false;
}


void Foam::cyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    polyPatch::initUpdateMesh(pBufs);

    if (boundaryMesh().mesh().topoChanging() && createAMIFaces_ && owner())
    {
        setAMIFaces();
    }

    updatingAMI_ = false;
}


void Foam::cyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    DebugInFunction << endl;

    // Note: this clears out cellCentres(), faceCentres() and faceAreas()
    polyPatch::updateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::clearGeom()
{
    DebugInFunction << endl;

    if (!updatingAMI_)
    {
//DEBUG("*** CLEARING AMI ***");
        AMIPtr_.clear();
    }
    polyPatch::clearGeom();
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    updated_(false),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    separationVector_(Zero),
    AMIPtr_(nullptr),
    AMIMethod_(AMIPatchToPatchInterpolation::imFaceAreaWeight),
    AMIReverse_(false),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(-1.0),
    surfPtr_(nullptr),
    surfDict_(fileName("surface")),
    createAMIFaces_(false),
    updatingAMI_(true),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_(),
    nbrFaceAreas0_(),
    nbrFaceCentres0_()
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    updated_(false),
    nbrPatchName_(dict.lookupOrDefault<word>("neighbourPatch", "")),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    rotationAngleDefined_(false),
    rotationAngle_(0.0),
    separationVector_(Zero),
    AMIPtr_(nullptr),
    AMIMethod_
    (
        AMIPatchToPatchInterpolation::interpolationMethodNames_
        [
            dict.lookupOrDefault
            (
                "method",
                AMIPatchToPatchInterpolation::interpolationMethodNames_
                [
                    AMIPatchToPatchInterpolation::imFaceAreaWeight
                ]
            )
        ]
    ),
    AMIReverse_(dict.lookupOrDefault("flipNormals", false)),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(dict.lookupOrDefault("lowWeightCorrection", -1.0)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface")),
    createAMIFaces_(false),
    updatingAMI_(true),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_(),
    nbrFaceAreas0_(),
    nbrFaceCentres0_()
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction(dict)
            << "No \"neighbourPatch\" or \"coupleGroup\" provided."
            << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    switch (transform())
    {
        case ROTATIONAL:
        {
            dict.readEntry("rotationAxis", rotationAxis_);
            dict.readEntry("rotationCentre", rotationCentre_);
            if (dict.readIfPresent("rotationAngle", rotationAngle_))
            {
                rotationAngleDefined_ = true;
                rotationAngle_ = degToRad(rotationAngle_);

                if (debug)
                {
                    Info<< "rotationAngle: " << rotationAngle_ << " [rad]"
                        <<  endl;
                }
            }

            scalar magRot = mag(rotationAxis_);
            if (magRot < SMALL)
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal rotationAxis " << rotationAxis_ << endl
                    << "Please supply a non-zero vector."
                    << exit(FatalIOError);
            }
            rotationAxis_ /= magRot;

            break;
        }
        case TRANSLATIONAL:
        {
            dict.readEntry("separationVector", separationVector_);
            break;
        }
        default:
        {
            // No additional info required
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible

    // If topology change, recover the sizes of the original patches
    label srcSize0 = 0;
    if (dict.readIfPresent("srcSize", srcSize0))
    {
        srcFaceIDs_.setSize(srcSize0);
    }
    label tgtSize0 = 0;
    if (dict.readIfPresent("tgtSize", tgtSize0))
    {
        tgtFaceIDs_.setSize(tgtSize0);
    }
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    updated_(false),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_),
    createAMIFaces_(false),
    updatingAMI_(true),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_(),
    nbrFaceAreas0_(),
    nbrFaceCentres0_()
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    updated_(false),
    nbrPatchName_(nbrPatchName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_),
    createAMIFaces_(false),
    updatingAMI_(true),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_(),
    nbrFaceAreas0_(),
    nbrFaceCentres0_()
{
    if (nbrPatchName_ == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    updated_(false),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    rotationAngleDefined_(pp.rotationAngleDefined_),
    rotationAngle_(pp.rotationAngle_),
    separationVector_(pp.separationVector_),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_),
    createAMIFaces_(false),
    updatingAMI_(true),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_(),
    nbrFaceAreas0_(),
    nbrFaceCentres0_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicAMIPolyPatch::neighbPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(neighbPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << neighbPatchName()
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic AMI patch
        const cyclicAMIPolyPatch& nbrPatch =
            refCast<const cyclicAMIPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.neighbPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << neighbPatchName()
                << nl << " but that in return specifies "
                << nbrPatch.neighbPatchName() << endl;
        }
    }

    return nbrPatchID_;
}


bool Foam::cyclicAMIPolyPatch::owner() const
{
    return index() < neighbPatchID();
}


const Foam::cyclicAMIPolyPatch& Foam::cyclicAMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];
    return refCast<const cyclicAMIPolyPatch>(pp);
}


const Foam::autoPtr<Foam::searchableSurface>&
Foam::cyclicAMIPolyPatch::surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && owner() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", name()));

        const polyMesh& mesh = boundaryMesh().mesh();

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


const Foam::AMIPatchToPatchInterpolation& Foam::cyclicAMIPolyPatch::AMI() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_.valid())
    {
        resetAMI(AMIMethod_);
    }

    return *AMIPtr_;
}


bool Foam::cyclicAMIPolyPatch::applyLowWeightCorrection() const
{
    if (owner())
    {
        return AMI().applyLowWeightCorrection();
    }
    else
    {
        return neighbPatch().AMI().applyLowWeightCorrection();
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition(pointField& l) const
{
    if (!parallel())
    {
        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(forwardT(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(forwardT(), l);
        }
    }
    else if (separated())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract

        const vectorField& s = separation();
        if (s.size() == 1)
        {
            forAll(l, i)
            {
                l[i] -= s[0];
            }
        }
        else
        {
            l -= s;
        }
    }
}


void Foam::cyclicAMIPolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            forwardT().size() == 1
          ? forwardT()[0]
          : forwardT()[facei]
        );

        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[facei]
        );

        l -= s;
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformPosition
(
    point& l,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        if (transform() == ROTATIONAL)
        {
            l = Foam::transform(T, l - rotationCentre_) + rotationCentre_;
        }
        else
        {
            l = Foam::transform(T, l);
        }
    }
    else if (separated())
    {
        const vector& s =
        (
            separation().size() == 1
          ? separation()[0]
          : separation()[facei]
        );

        l += s;
    }
}


void Foam::cyclicAMIPolyPatch::reverseTransformDirection
(
    vector& d,
    const label facei
) const
{
    if (!parallel())
    {
        const tensor& T =
        (
            reverseT().size() == 1
          ? reverseT()[0]
          : reverseT()[facei]
        );

        d = Foam::transform(T, d);
    }
}


void Foam::cyclicAMIPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{}


void Foam::cyclicAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::cyclicAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


Foam::label Foam::cyclicAMIPolyPatch::pointFace
(
    const label facei,
    const vector& n,
    point& p
) const
{
    point prt(p);
    reverseTransformPosition(prt, facei);

    vector nrt(n);
    reverseTransformDirection(nrt, facei);

    label nbrFacei = -1;

    if (owner())
    {
        nbrFacei = AMI().tgtPointFace
        (
            *this,
            neighbPatch(),
            nrt,
            facei,
            prt
        );
    }
    else
    {
        nbrFacei = neighbPatch().AMI().srcPointFace
        (
            neighbPatch(),
            *this,
            nrt,
            facei,
            prt
        );
    }

    if (nbrFacei >= 0)
    {
        p = prt;
    }

    return nbrFacei;
}


void Foam::cyclicAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);
    if (!nbrPatchName_.empty())
    {
        os.writeEntry("neighbourPatch", nbrPatchName_);
    }
    coupleGroup_.write(os);

    switch (transform())
    {
        case ROTATIONAL:
        {
            os.writeEntry("rotationAxis", rotationAxis_);
            os.writeEntry("rotationCentre", rotationCentre_);

            if (rotationAngleDefined_)
            {
                os.writeEntry("rotationAngle", radToDeg(rotationAngle_));
            }

            break;
        }
        case TRANSLATIONAL:
        {
            os.writeEntry("separationVector", separationVector_);
            break;
        }
        case NOORDERING:
        {
            break;
        }
        default:
        {
            // No additional info to write
        }
    }

    if (AMIMethod_ != AMIPatchToPatchInterpolation::imFaceAreaWeight)
    {
        os.writeEntry
        (
            "method",
            AMIPatchToPatchInterpolation::interpolationMethodNames_
            [
                AMIMethod_
            ]
        );
    }

    if (AMIReverse_)
    {
        os.writeEntry("flipNormals", AMIReverse_);
    }

    if (AMILowWeightCorrection_ > 0)
    {
        os.writeEntry("lowWeightCorrection", AMILowWeightCorrection_);
    }

    if (!surfDict_.empty())
    {
        surfDict_.writeEntry(surfDict_.dictName(), os);
    }

    if (createAMIFaces_)
    {
        os.writeEntry("srcSize", srcFaceIDs_.size());
        os.writeEntry("tgtSize", tgtFaceIDs_.size());
    }
}


// ************************************************************************* //

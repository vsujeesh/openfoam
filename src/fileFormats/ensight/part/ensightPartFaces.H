/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    Foam::ensightPartFaces

Description
    An implementation of ensightPart to hold mesh faces.

SourceFiles
    ensightPartFaces.C
    ensightPartFacesAddr.C
    ensightPartFacesIO.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartFaces_H
#define ensightPartFaces_H

#include "ensightPart.H"
#include "ensightFaces.H"
#include "polyMesh.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightPartFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightPartFaces
:
    public ensightFaces
{
    // Private Data

        //- The referenced pointField
        const pointField& points_;

        //- The referenced faces
        const faceList& faces_;

        //- Can skip local point renumbering when points are contiguous
        bool contiguousPoints_;


    // Demand-driven data

        //- The number of points used
        mutable label nPointsUsed_;

        //- Map global to local indices, with -1 for ununsed points
        mutable std::unique_ptr<labelList> localPointMap_;


    // Private Member Functions

        //- Calculate the points used
        void calcLocalPoints() const;

        //- Return the number of unique points (demand-driven)
        label getPointsUsed() const;

        //- Return the global to local indices (demand-driven)
        labelList& getPointMap() const;


        //- Helper: write connectivity
        static void writeConnectivity
        (
            ensightGeoFile& os,
            const ensightFaces::elemType etype,
            const UIndirectList<face>& faces,
            const labelUList& pointMap
        );

        //- Write processor-local (serial) geometry
        void writeSerial(ensightGeoFile& os) const;


        //- No copy construct
        ensightPartFaces(const ensightPartFaces&) = delete;

        //- No copy assignment
        void operator=(const ensightPartFaces&) = delete;


public:

    // Constructors

        //- Construct from points and faces.
        //  Can skip local point renumbering when points are contiguous
        ensightPartFaces
        (
            const string& description,
            const pointField& points,
            const faceList& faces,
            const bool contiguousPoints = false
        );

        //- Construct from polyMesh and polyPatch
        //- Part receives the name of the patch unless otherwise specified.
        ensightPartFaces
        (
            const polyMesh& mesh,
            const polyPatch& patch,
            const string& partName = string::null
        );

        //- Construct from polyPatch
        //- Part receives the name of the patch unless otherwise specified.
        explicit ensightPartFaces
        (
            const polyPatch& patch,
            const string& partName = string::null
        );


    //- Destructor
    virtual ~ensightPartFaces() = default;


    // Member Functions

        //- Clear any demand-driven data
        void clearOut();

        //- Write geometry (serial only)
        virtual void write(ensightGeoFile& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

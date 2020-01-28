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
    Foam::ensightPartCells

Description
    An implementation of ensightPart to hold volume mesh cells.

SourceFiles
    ensightPartCells.C
    ensightPartCellsAddr.C
    ensightPartCellsIO.C

\*---------------------------------------------------------------------------*/

#ifndef ensightPartCells_H
#define ensightPartCells_H

#include "ensightCells.H"
#include "polyMesh.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightPartCells Declaration
\*---------------------------------------------------------------------------*/

class ensightPartCells
:
    public ensightCells
{
    // Private Data

        //- The referenced mesh
        const polyMesh& mesh_;

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

        //- Element connectivity for polyhedrals
        void writePolysConnectivity
        (
            ensightGeoFile& os,
            const labelUList& addr,
            const labelUList& pointMap
        ) const;

        //- Element connectivity
        void writeConnectivity
        (
            ensightGeoFile& os,
            const ensightCells::elemType etype,
            const labelUList& idList,
            const labelUList& pointMap
        ) const;

        //- Write processor-local (serial) geometry
        void writeSerial(ensightGeoFile& os) const;


        //- No copy construct
        ensightPartCells(const ensightPartCells&) = delete;

        //- No copy assignment
        void operator=(const ensightPartCells&) = delete;


public:

    // Constructors

        //- Construct from entire polyMesh without zones.
        //- Part receives the specified name (default: "cells").
        explicit ensightPartCells
        (
            const polyMesh& mesh,
            const string& partName = "cells"
        );

        //- Construct a part from polyMesh and list of cells.
        //- Part receives the specified name (default: "cells").
        ensightPartCells
        (
            const polyMesh& mesh,
            const labelUList& cellIds,
            const string& partName = "cells"
        );

        //- Construct a part from polyMesh and selection of cells.
        //- Part receives the specified name (default: "cells").
        ensightPartCells
        (
            const polyMesh& mesh,
            const bitSet& selection,
            const string& partName = "cells"
        );

        //- Construct from cellZone.
        //- Part receives the name of the zone unless otherwise specified.
        ensightPartCells
        (
            const polyMesh& mesh,
            const cellZone& zn,
            const string& partName = string::null
        );


    //- Destructor
    virtual ~ensightPartCells() = default;


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

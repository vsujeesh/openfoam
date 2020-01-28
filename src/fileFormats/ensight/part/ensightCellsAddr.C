/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "ensightCells.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Map<Foam::label>
Foam::ensightCells::meshPointMap(const polyMesh& mesh) const
{
    const label nEstimate = 8*this->size();

    Map<label> pointMap(nEstimate);

    // Pass 1: markup used points from cells

    for (const label celli : this->cellIds())
    {
        for (const label facei : mesh.cells()[celli])
        {
            for (const label pointi : mesh.faces()[facei])
            {
                pointMap.insert(pointi, 0);
            }
        }
    }

    // Compact point numbering, preserves the original order
    label nPoints = 0;
    for (const label pointi : pointMap.sortedToc())
    {
        pointMap(pointi) = nPoints++;
    }

    return pointMap;
}


// ************************************************************************* //

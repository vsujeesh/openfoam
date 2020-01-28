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

#include "ensightPartCells.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPartCells::ensightPartCells
(
    const polyMesh& mesh,
    const string& partName
)
:
    ensightCells(),
    mesh_(mesh),
    contiguousPoints_(true),
    nPointsUsed_(0),
    localPointMap_(nullptr)
{
    if (!partName.empty())
    {
        rename(partName);
    }

    classify(mesh);
}


Foam::ensightPartCells::ensightPartCells
(
    const polyMesh& mesh,
    const labelUList& cellIds,
    const string& partName
)
:
    ensightCells(),
    mesh_(mesh),
    contiguousPoints_(false),
    nPointsUsed_(0),
    localPointMap_(nullptr)
{
    if (!partName.empty())
    {
        rename(partName);
    }

    classify(mesh, cellIds);
}


Foam::ensightPartCells::ensightPartCells
(
    const polyMesh& mesh,
    const bitSet& selection,
    const string& partName
)
:
    ensightCells(),
    mesh_(mesh),
    contiguousPoints_(false),
    nPointsUsed_(0),
    localPointMap_(nullptr)
{
    if (!partName.empty())
    {
        rename(partName);
    }

    classify(mesh, selection);
}


Foam::ensightPartCells::ensightPartCells
(
    const polyMesh& mesh,
    const cellZone& zn,
    const string& partName
)
:
    ensightPartCells
    (
        mesh,
        static_cast<const labelList&>(zn),
        zn.name()
    )
{
    identifier() = zn.index();

    if (!partName.empty())
    {
        rename(partName);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPartCells::clearOut()
{
    ensightCells::clearOut();

    nPointsUsed_ = 0;
    localPointMap_.reset(nullptr);
}


// ************************************************************************* //

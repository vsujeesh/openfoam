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

#include "ensightParts.H"
#include "bitSet.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightParts::ensightParts(const polyMesh& mesh)
{
    recalculate(mesh);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightParts::recalculate(const polyMesh& mesh)
{
    this->clear();

    // Track which cells are in a zone or not
    bitSet selection(mesh.nCells());

    // Do all cell zones
    for (const cellZone& zn : mesh.cellZones())
    {
        if (returnReduce(!zn.empty(), orOp<bool>()))
        {
            selection.set(zn);
            this->append(new ensightPartCells(mesh, zn));
        }
    }

    if (this->empty())
    {
        // No zones? - do entire mesh. Name as per ensightMesh
        this->insert(new ensightPartCells(mesh, "internalMesh"));
    }
    else
    {
        // Flip from zoned to unzoned
        selection.flip();

        if (returnReduce(selection.any(), orOp<bool>()))
        {
            // Place as first in the list
            this->insert
            (
                new ensightPartCells(mesh, selection, "__internal__")
            );
        }
    }


    for (const polyPatch& p : mesh.boundaryMesh())
    {
        // Only do real (non-processor) boundaries.
        if (isA<processorPolyPatch>(p))
        {
            break;
        }

        // Skip empty patch types and zero-sized patches
        // Would ideally like to check for emptyFvPatch,
        // but that is not available here
        if
        (
            !isA<emptyPolyPatch>(p)
         && returnReduce(!p.empty(), orOp<bool>())
        )
        {
            this->append(new ensightPartFaces(p));
        }
    }

    renumber();
}


void Foam::ensightParts::renumber(label start)
{
    for (ensightPart& part : *this)
    {
        part.index() = start;
        ++start;
    }
}


void Foam::ensightParts::write(ensightGeoFile& os) const
{
    // Some feedback
    Info<< "Write geometry part (" << flush;

    for (const ensightPart& part : *this)
    {
        Info<< ' ' << part.index() << flush;
        part.write(os);
    }

    Info<< " )" << endl;
}


void Foam::ensightParts::writeDict(Ostream& os, const bool full) const
{
    for (const ensightPart& part : *this)
    {
        part.writeDict(os, full);
    }
}


// ************************************************************************* //

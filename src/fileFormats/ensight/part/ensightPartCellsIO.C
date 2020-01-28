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

#include "ensightPartCells.H"
#include "ensightOutput.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightPartCells::writePolysConnectivity
(
    ensightGeoFile& os,
    const labelUList& addr,
    const labelUList& pointMap
) const
{
    // The number of faces per element
    {
        labelList send
        (
            ensightOutput::Detail::getPolysNFaces(mesh_, addr)
        );

        os.writeLabels(send);
    }

    // The number of points per element face
    {
        labelList send
        (
            ensightOutput::Detail::getPolysNPointsPerFace(mesh_, addr)
        );

        os.writeLabels(send);
    }

    ensightOutput::writePolysPoints
    (
        os,
        mesh_,
        addr,
        pointMap  // Remap point-ids
    );
}


void Foam::ensightPartCells::writeConnectivity
(
    ensightGeoFile& os,
    const ensightCells::elemType etype,
    const labelUList& addr,
    const labelUList& pointMap
) const
{
    if (addr.empty()) return;

    os.writeKeyword(ensightCells::key(etype));
    os.write(addr.size());
    os.newline();

    if (etype == ensightCells::NFACED)
    {
        writePolysConnectivity(os, addr, pointMap);
        return;
    }


    // Primitive shape - get subset and renumber
    cellShapeList shapes(mesh_.cellShapes(), addr);

    ListListOps::inplaceRenumber(pointMap, shapes);

    ensightOutput::writeCellShapes(os, shapes);
}


void Foam::ensightPartCells::writeSerial
(
    ensightGeoFile& os
) const
{
    const pointField& points = mesh_.points();

    if (size())
    {
        const labelUList& pointMap = getPointMap();
        const label nPoints = getPointsUsed();

        os.beginPart(index(), name());
        os.beginCoordinates(nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            forAll(pointMap, pti)
            {
                if (pointMap[pti] > -1)
                {
                    os.write(points[pti].component(cmpt));
                    os.newline();
                }
            }
        }

        // Write each element type
        for (int typei=0; typei < ensightCells::nTypes; ++typei)
        {
            const auto etype = ensightCells::elemType(typei);

            writeConnectivity
            (
                os,
                etype,
                cellIds(etype),
                pointMap
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPartCells::write(ensightGeoFile& os) const
{
    writeSerial(os);
}


// ************************************************************************* //

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

#include "ensightPartFaces.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightPartFaces::writeConnectivity
(
    ensightGeoFile& os,
    const ensightFaces::elemType etype,
    const UIndirectList<face>& faces,
    const labelUList& pointMap
)
{
    if (faces.empty()) return;

    os.writeKeyword(ensightFaces::key(etype));
    os.write(faces.size());
    os.newline();

    // Write (polygon) face sizes
    if (etype == ensightFaces::NSIDED)
    {
        // Write the number of points per face
        for (const face& f : faces)
        {
            os.write(f.size());
            os.newline();
        }
    }

    // Write the points describing the face
    // Convert global -> local index
    // (note: Ensight indices start with 1)
    for (const face& f : faces)
    {
        for (const label pointi : f)
        {
            os.write(pointMap[pointi] + 1);
        }
        os.newline();
    }
}


void Foam::ensightPartFaces::writeSerial(ensightGeoFile& os) const
{
    if (size())
    {
        const labelUList& pointMap = getPointMap();
        const label nPoints = getPointsUsed();

        os.beginPart(index(), name());
        os.beginCoordinates(nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            forAll(pointMap, pointi)
            {
                if (pointMap[pointi] > -1)
                {
                    os.write(points_[pointi].component(cmpt));
                    os.newline();
                }
            }
        }

        // Write part
        for (int typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const auto etype = ensightFaces::elemType(typei);

            writeConnectivity
            (
                os,
                etype,
                UIndirectList<face>(faces_, faceIds(etype)),
                pointMap
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPartFaces::write(ensightGeoFile& os) const
{
    writeSerial(os);
}



// ************************************************************************* //

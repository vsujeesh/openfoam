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

#include "ensightPartFaces.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightPartFaces::calcLocalPoints() const
{
    if (localPointMap_)
    {
        FatalErrorInFunction
            << "localPointMap_ already allocated"
            << abort(FatalError);
    }

    localPointMap_.reset(new labelList());

    labelList& usedPoints = *localPointMap_;

    usedPoints.resize(points_.size());

    if (contiguousPoints_)
    {
        nPointsUsed_ = usedPoints.size();

        ListOps::identity(usedPoints);

        return;
    }


    // Mark up with -1 for unused entries
    usedPoints = -1;
    label nPoints = 0;

    // Add all points from faces
    for (const label facei : this->faceIds())
    {
        for (const label pointi : faces_[facei])
        {
            if (usedPoints[pointi] == -1)
            {
                usedPoints[pointi] = nPoints++;
            }
        }
    }

    nPointsUsed_ = nPoints;

    // Compact point numbering, preserving the original order
    nPoints = 0;
    for (label& pointi : usedPoints)
    {
        if (pointi > -1)
        {
            pointi = nPoints++;
        }
    }
}


Foam::label Foam::ensightPartFaces::getPointsUsed() const
{
    if (!localPointMap_)
    {
        calcLocalPoints();
    }

    return nPointsUsed_;
}


Foam::labelList& Foam::ensightPartFaces::getPointMap() const
{
    if (!localPointMap_)
    {
        calcLocalPoints();
    }

    return *localPointMap_;
}


// ************************************************************************* //

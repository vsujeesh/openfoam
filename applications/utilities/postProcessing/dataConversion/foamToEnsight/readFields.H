/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Helper routines for reading a field or fields,
    for foamToEnsight

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef readFields_H
#define readFields_H

#include "instantList.H"
#include "IOobjectList.H"
#include "fvMeshSubsetProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Get the field or return nullptr
template<class GeoField>
tmp<GeoField> getField
(
    const IOobject* io,
    const typename GeoField::Mesh& mesh
)
{
    if (io)
    {
        return tmp<GeoField>::New(*io, mesh);
    }

    return nullptr;
}


//- Get the named field from the objects, or return nullptr.
template<class GeoField>
tmp<GeoField> getField
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const word& fieldName
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>(objects.findObject(fieldName), mesh);
}


//- Get internal field and make it a zero-gradient volume field
template<class GeoField>
tmp<GeoField>
getZeroGradField(const IOobject* io, const fvMesh& mesh)
{
    if (io)
    {
        auto tdimfield = tmp<typename GeoField::Internal>::New(*io, mesh);

        auto tfield = fvMeshSubsetProxy::zeroGradientField(tdimfield());
        tdimfield.clear();

        return tfield;
    }

    return nullptr;
}


//- Check if fields are good to use (available at all times)
//  ignore special fields (_0 fields),
//  ignore fields that are not available for all time-steps
label checkData
(
    const fvMesh& mesh,
    const instantList& timeDirs,
    wordList& objectNames
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

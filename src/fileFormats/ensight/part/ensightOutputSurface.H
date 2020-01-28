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

Class
    Foam::ensightOutputSurface

Description
    A variant of ensightFaces that holds references to continuous
    points/faces for writing out. The surface is assumed to have been merged
    prior so the output is serial-only.

SourceFiles
    ensightOutputSurface.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutputSurface_H
#define ensightOutputSurface_H

#include "ensightFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ensightOutputSurface Declaration
\*---------------------------------------------------------------------------*/

class ensightOutputSurface
:
    public ensightFaces
{
    // Private Data

        //- The referenced pointField
        const pointField& points_;

        //- The referenced faces
        const faceList& faces_;


    // Private Member Functions

        //- No copy construct
        ensightOutputSurface(const ensightOutputSurface&) = delete;

        //- No copy assignment
        void operator=(const ensightOutputSurface&) = delete;


public:

    // Constructors

        //- Construct from points and faces.
        ensightOutputSurface
        (
            const pointField& points,
            const faceList& faces
        );

        //- Construct named from points and faces.
        ensightOutputSurface
        (
            const string& description,
            const pointField& points,
            const faceList& faces
        );


    //- Destructor
    virtual ~ensightOutputSurface() = default;


    // Member Functions

        //- Write processor-local (serial) geometry
        virtual void write(ensightGeoFile& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

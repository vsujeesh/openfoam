/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::functionObjects::scale

Group
    grpFieldFunctionObjects

Description
    Multiplies an input field elementwise by an input scaling factor.

    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | PENDING        | PENDING
      output field  | PENDING        | PENDING
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    scale1
    {
        // Mandatory entries (unmodifiable)
        type            scale;
        libs            (fieldFunctionObjects);
        scale           1.0;

        // Mandatory (inherited) entry (runtime modifiable)
        field           <inpField>;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: log                     | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      field        | Name of the operand field          | word |  yes  | -
      scale        | Scaling factor                     | scalar | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::scale

SourceFiles
    scale.C
    scaleTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scale_H
#define functionObjects_scale_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class scale Declaration
\*---------------------------------------------------------------------------*/

class scale
:
    public fieldExpression
{
    // Private Data

        //- Scale factor
        scalar scale_;


    // Private Member Functions

        //- Calculate the scale of the field and register the result
        template<class Type>
        bool calcScale();

        //- Calculate the scale of the field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("scale");


    // Constructors

        //- Construct from Time and dictionary
        scale
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        scale(const scale&) = delete;

        //- No copy assignment
        void operator=(const scale&) = delete;


    //- Destructor
    virtual ~scale() = default;


    // Member Functions

        //- Read the randomise data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "scaleTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

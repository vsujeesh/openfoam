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
    Foam::functionObjects::age

Group
    grpFieldFunctionObjects

Description
    Computes the age of fluid in the domain in time units, i.e. the time taken
    for a fluid particle to travel from an inlet to an input location.

    Solves the following equation for incompressible flows:
    \f[
        \div (\phi t) = 1
    \f]

    where:
    \vartable
        t    | Age             [s]
        \phi | Volumetric flux [m^3/s]
    \endvartable

    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | PENDING        | PEND
      output field      | -              | -
    \endtable

Note
    Boundary conditions are generated automatically as zeroGradient on all
    walls and inletOutlet everywhere else.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    age1
    {
        // Mandatory entries (unmodifiable)
        type            age;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        schemesField    k;
        nCorr           5;
        tolerance       1e-5;
        diffusion       false;
        phi             phi;
        rho             rho;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: age                     | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      schemesField | Name of the field from which schemes are taken <!--
                                                    --> | word  | no   | age
      nCorr        | The maximum number of correctors that determines how <!--
      --> many times the solution is repeated to iterate away any <!--
      --> non-linearities due the choice of scheme      | label | no   | 5
      tolerance    | Solver residual control            | scalar | no  | 1e-5
      diffusion    | Flag to turn on/off the diffusion term | bool | no | false
      phi          | Name of the flux field             | word  | no   | phi
      rho          | Name of the density field          | word  | no   | rho
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Minimal example by using the \c postProcess utility:
    \verbatim
        postProcess -func age -fields "(phi)"
    \endverbatim

Note
    - \b Requirement: Requires a solver set in \c fvSolution and
    a \c div(phi\,age) in \c fvSchemes to be specified.
    - \b Warning: If the divergence scheme is set to \c upwind,
    no corrections will be necessary.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::age

SourceFiles
    age.C

\*---------------------------------------------------------------------------*/

#ifndef age_H
#define age_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                               Class age Declaration
\*---------------------------------------------------------------------------*/

class age
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Flag to turn on/off the diffusion term
        Switch diffusion_;

        //- Number of corrections
        label nCorr_;

        //- Name of the flux field
        word phiName_;

        //- Name of the density field
        word rhoName_;

        //- Name of field from which schemes are taken
        word schemesField_;

        //- Convergence tolerance
        scalar tolerance_;


    // Private Member Functions

        //- The list of patch types for the age field
        wordList patchTypes() const;

        //- Check convergence
        bool converged(const label nCorr, const scalar initialResidual) const;


public:

    //- Runtime type information
    TypeName("age");


    // Constructors

        //- Construct from Time and dictionary
        age
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        age(const age&) = delete;

        //- No copy assignment
        void operator=(const age&) = delete;


    //- Destructor
    virtual ~age() = default;


    // Member Functions

        //- Read the data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

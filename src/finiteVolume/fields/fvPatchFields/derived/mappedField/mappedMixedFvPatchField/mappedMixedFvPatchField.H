/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::mappedMixedFvPatchField

Group
    grpGenericBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition provides a self-contained version of e.g. mapped
    mixed boundary conditions.
    It does not use information on the patch; instead it holds the data locally.

Usage
    \table
        Property     | Description                      | Required | Default
        field        | Name of field to be mapped       | no  | this field name
        weight       | Name of field to use as weight   | yes |
        setAverage   | Use average value                | no  | false |
        average      | Average value to use if \c setAverage = yes | partly |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedMixed;
        field           T;              // optional field name
        setAverage      no;             // apply an average value
        average         0;              // average to apply if setAverage
        value           uniform 0;      // place holder
    }
    \endverbatim

Note
    Since this condition can be applied on a per-field and per-patch basis,
    it is possible to duplicate the mapping information.  If possible, employ
    the \c mapped condition in preference to avoid this situation, and only
    employ this condition if it is not possible to change the underlying
    geometric (poly) patch type to \c mapped.

See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::mappedFieldFvPatchField
    Foam::mixedFvPatchField

SourceFiles
    mappedMixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedMixedFvPatchField_H
#define mappedMixedFvPatchField_H

#include "mappedPatchBase.H"
#include "mappedPatchFieldBase.H"
#include "mixedFvPatchFields.H"
#include "interpolation.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mappedMixedFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedMixedFvPatchField
:
    public mixedFvPatchField<Type>,
    public mappedPatchBase,
    public mappedPatchFieldBase<Type>
{

        //- Name of weight field to sample
        word weightFieldName_;


public:

    //- Runtime type information
    TypeName("mappedMixed");


    // Constructors

        //- Construct from patch and internal field
        mappedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

//        //- Construct from patch, internal field and distance for normal type
//        //  sampling
//        mappedMixedFvPatchField
//        (
//            const fvPatch&,
//            const DimensionedField<Type, volMesh>&,
//
//            // mappedPatchBase
//            const word& sampleRegion,
//            const sampleMode sampleMode,
//            const word& samplePatch,
//            const scalar distance,
//
//            // My settings
//            const word& fieldName,
//            const bool setAverage,
//            const Type average,
//            const word& interpolationScheme
//        );

        //- Construct by mapping given
        //  mappedMixedFvPatchField
        //  onto a new patch
        mappedMixedFvPatchField
        (
            const mappedMixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedMixedFvPatchField
        (
            const mappedMixedFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedMixedFvPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        mappedMixedFvPatchField
        (
            const mappedMixedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedMixedFvPatchField<Type>
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

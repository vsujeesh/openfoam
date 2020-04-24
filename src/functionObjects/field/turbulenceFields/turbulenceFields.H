/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::functionObjects::turbulenceFields

Group
    grpFieldFunctionObjects

Description
    Computes and stores various turbulence-related quantities
    on the mesh database for further manipulation.

    Fields are stored as copies of the original, with the prefix
    "tubulenceModel:", e.g.:

    \verbatim
      turbulenceModel:R
    \endverbatim

    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | -              | -
      output field      | vol<Type>Field | $FOAM_CASE/<time>/<outField>
    \endtable

Usage
    Minimal example by using \c system/controlDict:functions comprise:
    \verbatim
    turbulenceFields1
    {
        // Mandatory entries (unmodifiable)
        type        turbulenceFields;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        // Either field or fields entries
        fields      (R devRhoReff);
        field       R;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: turbulenceFields           | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      fields       | Names of fields to store (see below)  | wordList | yes | -
      field        | Name of a field to store (see below)  | word | yes | -
    \endtable

    where \c fields can include:
    \plaintable
      k           | turbulent kinetic energy
      epsilon     | turbulent kinetic energy dissipation rate
      omega       | specific dissipation rate
      nuTilda     | modified turbulent viscosity
      nut         | turbulent viscosity (incompressible)
      nuEff       | effective turbulent viscosity (incompressible)
      mut         | turbulent viscosity (compressible)
      muEff       | effective turbulent viscosity (compressible)
      alphat      | turbulence thermal diffusivity (compressible)
      alphaEff    | effective turbulence thermal diffusivity (compressible)
      R           | Reynolds stress tensor
      devReff     | deviatoric part of the effective Reynolds stress
      devRhoReff  | divergence of the Reynolds stress
      L           | integral-length scale / mixing-length scale
      I           | turbulence intensity
    \endplaintable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
    The complete documentation is given in the Extended Code Guide.


    Minimal example by using the \c postProcess utility is as follows:
    \verbatim
        <solver> -postProcess -func turbulenceFields
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl

SourceFiles
    turbulenceFields.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_turbulenceFields_H
#define functionObjects_turbulenceFields_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"
#include "Enum.H"
#include "volFieldsFwd.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class turbulenceFields Declaration
\*---------------------------------------------------------------------------*/

class turbulenceFields
:
    public fvMeshFunctionObject
{
public:

    // Public Enumerations

        //- Options for the turbulence fields (compressible)
        enum compressibleField
        {
            cfK,            //!< "Turbulent kinetic energy"
            cfEpsilon,      //!< "Turbulent kinetic energy dissipation rate"
            cfOmega,        //!< "Specific dissipation rate"
            cfNuTilda,      //!< "Modified turbulent viscosity"
            cfMut,          //!< "Turbulent viscosity"
            cfMuEff,        //!< "Effective turbulent viscosity"
            cfAlphat,       //!< "Turbulence thermal diffusivity"
            cfAlphaEff,     //!< "Effective turbulence thermal diffusivity"
            cfR,            //!< "Reynolds stress tensor"
            cfDevRhoReff,   //!< "Divergence of the Reynolds stress"
            cfI             //!< "Turbulence intensity"
            cfL,            //!< "Integral-length/Mixing-length scale"
            cfT,            //!< "Integral time scale"
            cfU,            //!< "Integral velocity scale"
            cfLt,           //!< "Taylor length scale"
            cfTt,           //!< "Taylor time scale"
            cfUt,           //!< "Taylor velocity scale"
            cfLk,           //!< "Kolmogorov length scale"
            cfTk,           //!< "Kolmogorov time scale"
            cfUk            //!< "Kolmogorov velocity scale"
        };

        //- Names for compressibleField turbulence fields
        static const Enum<compressibleField> compressibleFieldNames_;

        //- Options for the turbulence fields (incompressible)
        enum incompressibleField
        {
            ifK,            //!< "Turbulent kinetic energy"
            ifEpsilon,      //!< "Turbulent kinetic energy dissipation rate"
            ifOmega,        //!< "Specific dissipation rate"
            ifNuTilda,      //!< "Modified turbulent viscosity"
            ifNut,          //!< "Turbulent viscosity"
            ifNuEff,        //!< "Effective turbulent viscosity"
            ifR,            //!< "Reynolds stress tensor"
            ifDevReff,   //!< "Deviatoric part of the effective Reynolds stress"
            ifI             //!< "Turbulence intensity"
            ifL,            //!< "Integral-length/Mixing-length scale"
            ifT,            //!< "Integral time scale"
            ifU,            //!< "Integral velocity scale"
            ifLt,           //!< "Taylor length scale"
            ifTt,           //!< "Taylor time scale"
            ifUt,           //!< "Taylor velocity scale"
            ifLk,           //!< "Kolmogorov length scale"
            ifTk,           //!< "Kolmogorov time scale"
            ifUk            //!< "Kolmogorov velocity scale"
        };

        //- Names for incompressibleField turbulence fields
        static const Enum<incompressibleField> incompressibleFieldNames_;

        //- Turbulence closure model name
        static const word modelName_;


private:

    // Private Data

        //- Internal flags to determine the base parameter of turbulence model
        enum class modelBase
        {
            EPSILON,
            OMEGA,
            OTHER
        };

        //- Internal flag to determine the base parameter of turbulence model
        modelBase modelBase_;


    // Private Member Functions

        //- No copy construct
        turbulenceFields(const turbulenceFields&) = delete;

        //- No copy assignment
        void operator=(const turbulenceFields&) = delete;


protected:

    // Protected Data

        //- Fields to load
        wordHashSet fieldSet_;


    // Protected Member Functions

        //- Return true if compressible turbulence model is identified
        bool compressible();

        //- Process the turbulence field
        template<class Type>
        void processField
        (
            const word& fieldName,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
        );

        //- Return nuTilda calculated from k and omega
        template<class Model>
        tmp<volScalarField> nuTilda(const Model& model) const;

        //- Return turbulence intensity, I, calculated from k and U
        template<class Model>
        tmp<volScalarField> I(const Model& model) const;

        //- Return integral length scale, L
        template<class Model>
        tmp<volScalarField> L(const Model& model) const;

        //- Return integral time scale, T
        template<class Model>
        tmp<volScalarField> T(const Model& model) const;

        //- Return integral velocity scale, U
        template<class Model>
        tmp<volScalarField> U(const Model& model) const;

        //- Return Taylor length scale, Lt
        template<class Model>
        tmp<volScalarField> Lt(const Model& model) const;

        //- Return Taylor time scale, Tt
        template<class Model>
        tmp<volScalarField> Tt(const Model& model) const;

        //- Return Taylor velocity scale, Ut
        template<class Model>
        tmp<volScalarField> Ut(const Model& model) const;

        //- Return Kolmogorov length scale, Lk
        template<class Model>
        tmp<volScalarField> Lk(const Model& model) const;

        //- Return Kolmogorov time scale, Tk
        template<class Model>
        tmp<volScalarField> Tk(const Model& model) const;

        //- Return Kolmogorov velocity scale, Uk
        template<class Model>
        tmp<volScalarField> Uk(const Model& model) const;


public:

    //- Runtime type information
    TypeName("turbulenceFields");


    // Constructors

        //- Construct from Time and dictionary
        turbulenceFields
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~turbulenceFields() = default;


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Calculate turbulence fields
        virtual bool execute();

        //- Do nothing.
        //  The turbulence fields are registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulenceFieldsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

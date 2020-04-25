/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::functionObjects::comfort

Group
    grpFieldFunctionObjects

Description
    Computes the thermal comfort indicators, i.e. the predicted mean vote (PMV)
    and predicted percentage of dissatisfaction (PPD) based on
    DIN ISO EN 7730:2005.

    The values of the \c PMV mean:
    \table
      Predicted Mean Vote (PMV)   | Evaluation
      + 3                         | hot
      + 2                         | warm
      + 1                         | slightly warm
      + 0                         | neutral
      - 1                         | slightly cool
      - 2                         | cool
      - 3                         | cold
    \endtable

    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | PENDING        | PENDING
      output field      | PENDING        | PENDING
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    comfort1
    {
        // Mandatory entries (unmodifiable)
        type            comfort;
        libs            (fieldFunctionObjects);

        // Optional entries (runtime modifiable)
        clothing        0.0;
        metabolicRate   0.8;
        extWork         0.0;
        Trad            0.0;
        relHumidity     0.5;
        pSat            0.0;
        tolerance       1e-4;
        maxClothIter    100;
        meanVelocity    false;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                        | Type | Req'd | Dflt
      type          | Type name: comfort                 | word |  yes  | -
      libs          | Library name: fieldFunctionObjects | word |  yes  | -
      clothing      | The insulation value of the cloth [-] | scalar | no  | 0.0
      metabolicRate | The metabolic rate           [kg/s^3] | scalar | no  | 0.8
      extWork       | The external work            [kg/s^3] | scalar | no  | 0.0
      Trad          | Mean radiation temperature   [K]      | scalar | no  | 0.0
      relHumidity   | Relative humidity of the air [%]      | scalar | no  | 0.5
      pSat          | Saturation pressure of water [Pa]     | scalar | no  | 0.0
      tolerance     | Tolerance criteria for iterative process to find Tcl <!--
                                                        --> | scalar | no | 1e-4
      maxClothIter  | Maximum number of correctors for cloth temperature   <!--
                                                        --> | label  | no | 100
      meanVelocity  | Flag to use volume weighted velocity field or a <!--
                  --> constant mean velocity for calculation in the whole <!--
                  --> domain                                | bool | no | false
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::comfort

SourceFiles
    comfort.C

\*---------------------------------------------------------------------------*/

#ifndef comfort_H
#define comfort_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class comfort Declaration
\*---------------------------------------------------------------------------*/

class comfort
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Flag to set to true if the radiation temperature (Trad) is provided
        Switch TradSet_;

        //- Flag to use volume weighted velocity field for calculation
        Switch meanVelocity_;

        //- Maximum number of correctors for cloth temperature
        label maxClothIter_;

        //- Tolerance criteria for iterative process to find Tcl
        scalar tolerance_;

        //- The insulation value of the cloth
        dimensionedScalar clothing_;

        //- Metabolic rate
        dimensionedScalar metabolicRate_;

        //- External work
        dimensionedScalar extWork_;

        //- Mean radiation temperature
        dimensionedScalar Trad_;

        //- Relative humidity
        dimensionedScalar relHumidity_;

        //- Saturation pressure of water
        dimensionedScalar pSat_;

        //- Thermal insulation of clothing [W/m^2/K]
        dimensionedScalar Icl_;

        //- Prefactor of cloth area [-]
        dimensionedScalar fcl_;


    // Private Member Functions

        //- Calculate the magnitude of the velocity [m/s]
        tmp<volScalarField> magU() const;

        //- Calculate the radiation temperature in the domain using a simple
        //- approach [K]
        dimensionedScalar Trad() const;

        //- Calculate the saturation pressure based on 7730:2006
        //  Possible options: adding different calculation methods such as
        //  the formulation based on Magnus or others [Pa]
        tmp<volScalarField> pSat() const;

        //- Calculate and return the surface temperature of the cloth [K]
        //- and the heat transfer coefficient hc [W/m^2/K]
        tmp<volScalarField> Tcloth
        (
            volScalarField& hc,
            const dimensionedScalar& metabolicRateSI,
            const dimensionedScalar& extWorkSI,
            const volScalarField& TdegC,
            const dimensionedScalar& Trad
        );

        //- Return true if the cloth temperature iteration has converged
        bool converged(const volScalarField&) const;


public:

    //- Runtime type information
    TypeName("comfort");


    // Constructors

        //- Construct from Time and dictionary
        comfort
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        comfort(const comfort&) = delete;

        //- No copy assignment
        void operator=(const comfort&) = delete;


    //- Destructor
    virtual ~comfort() = default;


    // Member Functions

        //- Read the data needed for the comfort calculation
        virtual bool read(const dictionary&);

        //- Calculate the predicted mean vote (PMV)
        //- and predicted percentage dissatisfaction (PPD) fields
        virtual bool execute();

        //- Write the PPD and PMV fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

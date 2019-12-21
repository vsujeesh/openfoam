/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::RASModels::realizableKE

Group
    grpRASTurbulence

Description
    Realizable k-epsilon turbulence model for incompressible and compressible
    flows.

    \heading Input fields
    \plaintable
        k        | Turbulent kinetic energy [m2/s2]
        epsilon  | Turbulent kinetic energy dissipation rate [m2/s3]
    \endplaintable

    Reference:
    \verbatim
        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1995).
        A new k-epsilon eddy viscosity model for high Reynolds number
        turbulent flows.
        Computers and Fluids, 24(3), 227-238.
        DOI:10.1016/0045-7930(94)00032-T
    \endverbatim

    The default model coefficients are (SLSYZ:Table 2):
    \verbatim
        realizableKECoeffs
        {
            A0          4.0;
            C2          1.9;
            sigmak      1.0;
            sigmaEps    1.2;
        }
    \endverbatim

SourceFiles
    realizableKE.C

\*---------------------------------------------------------------------------*/

#ifndef realizableKE_H
#define realizableKE_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class realizableKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class realizableKE
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        realizableKE(const realizableKE&) = delete;

        //- No copy assignment
        void operator=(const realizableKE&) = delete;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar A0_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;

            //-
            volScalarField Cmu_;


   // Protected Member Functions

        volScalarField Cmu
        (
            const volSymmTensorField& S,
            const volScalarField& S2,
            const volScalarField& magS,
            const volScalarField& magOmega2
        );

        //- Update nut with the latest available k, and epsilon
        virtual void correctNut();

        virtual tmp<fvScalarMatrix> kSource() const;

        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("realizableKE");

    // Constructors

        //- Construct from components
        realizableKE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~realizableKE() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/sigmak_ + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    (this->nut_/sigmaEps_ + this->nu())
                )
            );
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the (estimated) specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return tmp<volScalarField>::New
            (
                IOobject
                (
                    IOobject::groupName("omega", this->alphaRhoPhi_.group()),
                    this->runTime_.timeName(),
                    this->mesh_
                ),
                epsilon_/(0.09*k_)
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "realizableKE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

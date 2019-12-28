/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::RASModels::kOmega

Group
    grpRASTurbulence

Description
    The standard k-omega turbulence closure model for high Reynolds-number
    incompressible and compressible flows.

    The model is a two-transport-equation linear-eddy-viscosity turbulence
    closure model:
      - Turbulent kinetic energy, \c k,
      - Specific dissipation rate, \c omega.

    Reference:
    \verbatim
        Standard model (Tag:W):
            Wilcox, D. C. (2008).
            Formulation of the k-w turbulence model revisited.
            AIAA journal, 46(11), 2823-2838.
            DOI:10.2514/1.36541

        Low Reynolds-number modification (Tag:R):
            Rumsey, C. (2016).
            The Wilcox k-omega Turbulence Model.
            [online] Turbulence Modeling Resource.
            Available at: https://turbmodels.larc.nasa.gov/wilcox.html
            [Accessed 28 Dec. 2019].

            Wilcox, D. C. (2006).
            Turbulence Modeling for CFD (3rd ed.).
            La Canada, California: DCW Industries.
            ISBN-13:978-1-928729-08-2
    \endverbatim

    The default model coefficients are (W:Eqs. 7, 10-13):
    \verbatim
        kOmegaCoeffs
        {
            betaStar    0.09;    // equivalent to 'Cmu'
            gamma       0.52;    // equivalent to 'alpha'
            beta0       0.0708;
            CLim        0.875;
            sigmaD      0.125;
            sigmaK      0.6;
            sigmaOmega  0.5;
        }
    \endverbatim

Note:
    Please note that prior to OpenFOAM-v2006, the implementation of k-omega
    model was based on the previous version of the model:
    \verbatim
        Wilcox, D. C. (1998).
        Turbulence modeling for CFD
        (Vol. 2, pp. 103-217). La Canada, CA: DCW industries.
    \endverbatim

SourceFiles
    kOmega.C

\*---------------------------------------------------------------------------*/

#ifndef kOmega_H
#define kOmega_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmega
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kOmega(const kOmega&) = delete;

        //- No copy assignment
        void operator=(const kOmega&) = delete;

        //- Round-jet function (W:Eqs. 12-13)
        volScalarField::Internal fBetaFunc
        (
            const volTensorField::Internal& gradU
        ) const;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar betaStar_;
            dimensionedScalar gamma_;
            dimensionedScalar beta0_;
            dimensionedScalar CLim_;
            dimensionedScalar sigmaD_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaOmega_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Specific dissipation rate [1/s]
            volScalarField omega_;


    // Protected Member Functions

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmega");


    // Constructors

        //- Construct from components
        kOmega
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
    virtual ~kOmega() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k (W:Eq. 8)
        //  Variable change: rho*k/omega ~ nut
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nu() + sigmaK_*k_/omega_
                )
            );
        }

        //- Return the effective diffusivity for omega (W:Eq. 9)

        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    this->nu() + sigmaOmega_*k_/omega_
                )
            );
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the specific dissipation rate field
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    betaStar_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

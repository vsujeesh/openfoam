/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::incompressible::RASModels::kkLOmega

Group
    grpIcoRASTurbulence

Description
    Low Reynolds-number k-kl-omega turbulence model for
    incompressible transitional flows.

    The model is a three-transport-equation linear-eddy-viscosity turbulence
    closure model:
      - Laminar kinetic energy, \c kl,
      - Turbulent kinetic energy, \c kt,
      - Specific dissipation rate, \c omega.

    Reference:
    \verbatim
        Standard model (tag:WC):
            Walters, D. K., & Cokljat, D. (2008).
            A three-equation eddy-viscosity model for Reynolds-averaged
            Navier–Stokes simulations of transitional flow.
            Journal of Fluids Engineering, 130(12), 121401.
            DOI:10.1115/1.2979230

        Typo corrections (tag:F):
            F{\"u}rst, J. (2013).
            Numerical simulation of transitional flows
            with laminar kinetic energy.
            Engineering MECHANICS, 20(5), 379-388.
    \endverbatim

    The default model coefficients are (F:p. 383):
    \verbatim
        kkLOmegaCoeffs
        {
            A0             4.04;
            As             2.12;
            Av             6.75;
            Abp            0.6;
            Anat           200;
            Ats            200;
            CbpCrit        1.2;
            Cnc            0.1;
            CnatCrit       1250;
            Cint           0.75;
            CtsCrit        1000;
            CrNat          0.02;
            C11            3.4e-6;
            C12            1.0e-10;
            CR             0.12;
            CalphaTheta    0.035;
            Css            1.5;
            CtauL          4360;
            Cw1            0.44;
            Cw2            0.92;
            Cw3            0.3;
            CwR            1.5;
            Clambda        2.495;
            CmuStd         0.09;
            Prtheta        0.85;
            sigmak         1;
            sigmaOmega     1.17;
        }
    \endverbatim

SourceFiles
    kkLOmega.C

\*---------------------------------------------------------------------------*/

#ifndef kkLOmega_H
#define kkLOmega_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kkLOmega Declaration
\*---------------------------------------------------------------------------*/

class kkLOmega
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private Member Functions

        //- No copy construct
        kkLOmega(const kkLOmega&) = delete;

        //- No copy assignment
        void operator=(const kkLOmega&) = delete;

        tmp<volScalarField::Internal> fv
        (
            const volScalarField::Internal& Ret
        ) const;

        tmp<volScalarField::Internal> fINT() const;

        tmp<volScalarField::Internal> fSS
        (
            const volScalarField::Internal& Omega
        ) const;

        tmp<volScalarField> Cmu(const volScalarField& S) const;

        tmp<volScalarField> BetaTS(const volScalarField& Rew) const;

        tmp<volScalarField::Internal> fTaul
        (
            const volScalarField::Internal& lambdaEff,
            const volScalarField::Internal& ktL,
            const volScalarField& omega
        ) const;

        tmp<volScalarField> alphaT
        (
            const volScalarField& lambdaEff,
            const volScalarField& fv,
            const volScalarField& ktS
        ) const;

        tmp<volScalarField> fOmega
        (
            const volScalarField& lambdaEff,
            const volScalarField& lambdaT
        ) const;

        tmp<volScalarField> phiBP(const volScalarField& omega) const;

        tmp<volScalarField> phiNAT
        (
            const volScalarField& ReOmega,
            const volScalarField& fNatCrit
        ) const;

        tmp<volScalarField> D(const volScalarField& k) const;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar A0_;
            dimensionedScalar As_;
            dimensionedScalar Av_;
            dimensionedScalar Abp_;
            dimensionedScalar Anat_;
            dimensionedScalar Ats_;
            dimensionedScalar CbpCrit_;
            dimensionedScalar Cnc_;
            dimensionedScalar CnatCrit_;
            dimensionedScalar Cint_;
            dimensionedScalar CtsCrit_;
            dimensionedScalar CrNat_;
            dimensionedScalar C11_;
            dimensionedScalar C12_;
            dimensionedScalar CR_;
            dimensionedScalar CalphaTheta_;
            dimensionedScalar Css_;
            dimensionedScalar CtauL_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar CwR_;
            dimensionedScalar Clambda_;
            dimensionedScalar CmuStd_;
            dimensionedScalar Prtheta_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaOmega_;


        // Fields

            //- Laminar kinetic energy [m2/s2]
            volScalarField kt_;

            //- Turbulent kinetic energy [m2/s2]
            volScalarField kl_;

            //- Specific dissipation rate [1/s]
            volScalarField omega_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("kkLOmega");


    // Constructors

        //- Construct from components
        kkLOmega
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kkLOmega() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphaT/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DomegaEff", alphaT/sigmaOmega_ + nu())
            );
        }

        //- Return the laminar kinetic energy field
        virtual tmp<volScalarField> kl() const
        {
            return kl_;
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> kt() const
        {
            return kt_;
        }

        //- Return the specific dissipation rate field
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the total kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    kt_ + kl_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the total kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    kt_*omega_ + D(kl_) + D(kt_),
                    omega_.boundaryField().types()
                )
            );
        }

        //- Validate the turbulence fields after construction
        //  Update turbulent viscosity and other derived fields as requires
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

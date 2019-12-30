/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::RASModels::kEpsilonLRN

Group
    grpRASTurbulence

Description
    A suite of low Reynolds-number (LRN) k-epsilon turbulence closure models
    for incompressible and compressible flows.

    The suite consists of the following model variants:
    \vartable
        LaunderSharma    | Launder, B. E., & Sharma, B. I. (1974).
        LamBremhorst     | Lam, C. K. G., & Bremhorst, K. (1981).
        LienLeschziner   | Lien, F. S., & Leschziner, M. A. (1993).
        Chien            | Chien (1982).
        MyongKasagi      | Myong and Kasagi (1990).
        AbeKondohNagano  | Abe, Kondoh and Nagano (1994).
    \endvartable

    \heading Input fields
    \plaintable
        k        | Turbulent kinetic energy [m2/s2]
        epsilon  | Turbulent kinetic energy dissipation rate [m2/s3]
    \endplaintable

    References:
    \verbatim
        Summary of model coefficients and
        damping functions of variants (tags:HBCS; WM):
            Hrenya, C.M., Bolio, E.J., Chakrabarti, D., & Sinclair, J.L. (1995).
            Comparison of low Reynolds number k−ε turbulence models in
            predicting fully developed pipe flow.
            Chemical Engineering Science, 50(12), 1923-1941.
            DOI:10.1016/0009-2509(95)00035-4

            Wang, S. J., & Mujumdar, A. S. (2005).
            A comparative study of five low Reynolds number
            k–ε models for impingement heat transfer.
            Applied thermal engineering, 25(1), 31-44.
            DOI: 10.1016/j.applthermaleng.2004.06.001

        LaunderSharma model (tag:LS):
            Launder, B. E., & Sharma, B. I. (1974).
            Application of the energy-dissipation model of turbulence to the
            calculation of flow near a spinning disc.
            Letters in heat and mass transfer, 1(2), 131-137.

        LamBremhorst model (tag:LB):
            Lam, C. K. G., & Bremhorst, K. (1981).
            A modified form of the k-ε model for predicting wall turbulence.
            Journal of Fluids Engineering, 103(3), 456-460.
            DOI:10.1115/1.3240815

        LienLeschziner model (tag:LL):
            Lien, F. S., & Leschziner, M. A. (1993).
            A pressure-velocity solution strategy for compressible flow
            and its application to shock/boundary-layer interaction
            using second-moment turbulence closure.
            Journal of fluids engineering, 115(4), 717-725.
            DOI:10.1115/1.2910204

        Chien model (tag:C):
            Chien, K. Y. (1982).
            Predictions of channel and boundary-layer flows with
            a low-Reynolds-number turbulence model.
            AIAA journal, 20(1), 33-38.
            DOI:10.2514/3.51043

        MyongKasagi model (tag:MK):
            Myong, H. K., & Kasagi, N. (1990).
            A new approach to the improvement of k-ε turbulence model
            for wall-bounded shear flows.
            JSME Int. J. (Series II) 33, 63-72.
            DOI:10.1299/jsmeb1988.33.1_63

        AbeKondohNagano model (tag:AKN):
            Abe, K., Kondoh, T., & Nagano, Y. (1994).
            A new turbulence model for predicting fluid flow and heat transfer
            in separating and reattaching flows—I. Flow field calculations.
            International journal of heat and mass transfer, 37(1), 139-151.
            DOI:10.1016/0017-9310(94)90168-6
    \endverbatim

    The default model coefficients are:
    \verbatim
        kEpsilonLRNCoeffs
        {
            // Model variant
            variant     LaunderSharma;    // please refer above for options

            // Common model coefficients
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          0.0;
            sigmak      1.0;
            sigmaEps    1.0;
        }
    \endverbatim

SourceFiles
    kEpsilonLRN.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonLRN_H
#define kEpsilonLRN_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kEpsilonLRN Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonLRN
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kEpsilonLRN(const kEpsilonLRN&) = delete;

        //- No copy assignment
        void operator=(const kEpsilonLRN&) = delete;


    // Private Enumerations

        //- Options for the low Reynolds-number kEpsilon model variants
        enum variantType
        {
            LAUNDER_SHARMA,    //!< "Launder-Sharma (1974)"
            LAM_BREMHORST,     //!< "Lam-Bremhorst (1981)"
            LIEN_LESCHZINER,   //!< "Lien-Leschziner (1993)"
            CHIEN,             //!< "Chien (1982)"
            MYONG_KASAGI,      //!< "Myong-Kasagi (1990)"
            ABE_KONDOH_NAGANO, //!< "Abe-Kondoh-Nagano (1994)""
        };

        //- Names for variantType
        static const Enum<variantType> variantTypeNames;


    // Private Data

        //- Model variant (default=LAUNDER_SHARMA)
        enum variantType variant_;


protected:

    // Protected Data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;


        // Fields (input)

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;


        // Fields (internal)

            //- Model functions to account for low Reynolds-number effects
            volScalarField fMu_;
            volScalarField f1_;
            volScalarField f2_;

            //- Diffusion-balancing for near-wall region 'k_' (HBCS:p. 1928)
            volScalarField D_;

            //- Dissipation-rate-balancing
            //- for near-wall region 'k_' (HBCS: p. 1928)
            volScalarField E_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Private Member Functions

        //- Get input values, or use default values for model coefficients
        void variantCoeffs();

        //- Compute model functions, i.e. fMu_, f1_ and f2_
        void evaluateModelFuncs();

        //- Compute diffusion/dissipation rate balancers, i.e. D_ and E_
        void evaluateBalancers();

        //- Update nut with the latest available k, epsilon, and model funcs
        virtual void correctNut();

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonLRN");


    // Constructors

        //- Construct from components
        kEpsilonLRN
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
    virtual ~kEpsilonLRN() = default;


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

        //- Solve the turbulence equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonLRN.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

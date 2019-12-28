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

\*---------------------------------------------------------------------------*/

#include "kOmega.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
volScalarField::Internal kOmega<BasicTurbulenceModel>::fBetaFunc
(
    const volTensorField::Internal& gradU
) const
{
    // Rotation tensor (W:p. 2823)
    volTensorField::Internal Omega(skew(gradU));

    // Galilean-invariant Favre-averaged strain-rate tensor (W:p. 2823)
    volSymmTensorField::Internal SHat(symm(gradU) - 0.5*tr(gradU)*I);

    // Absolute value of Pope’s nondimensional measure of
    // vortex stretching parameter (W:p. 2823; Eq. 13)
    volScalarField::Internal chiOmega
    (
        mag(((Omega & Omega) && SHat)/pow3(betaStar_*omega_()))
    );

    // Round-jet function (W:Eq. 12)
    return (1.0 + 85.0*chiOmega)/(1.0 + 100.0*chiOmega);
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kOmega<BasicTurbulenceModel>::correctNut()
{
    // (W:Eq. 6)
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmega<BasicTurbulenceModel>::kOmega
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    gamma_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52    // = 13/25
        )
    ),
    beta0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta0",
            this->coeffDict_,
            0.0708
        )
    ),
    CLim_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "CLim",
            this->coeffDict_,
            0.875    // = 7/8
        )
    ),
    sigmaD_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaD",
            this->coeffDict_,
            0.125    // = 1/8
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            0.6
        )
    ),
    sigmaOmega_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaOmega",
            this->coeffDict_,
            0.5
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmega<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        betaStar_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        beta0_.readIfPresent(this->coeffDict());
        CLim_.readIfPresent(this->coeffDict());
        sigmaD_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaOmega_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kOmega<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    // Zero-trace Favre-averaged strain-rate tensor (W:Eq. 5)
    const volSymmTensorField SBar(dev(symm(gradU)));

    // Round-jet function
    const volScalarField::Internal fBeta(fBetaFunc(gradU()));

    const volScalarField::Internal GbyNu((2.0*SBar()) && gradU());

    const volScalarField::Internal G(this->GName(), nut()*GbyNu);

    // Cross-diffusion term (W:Eq. 11; p. 2825)
    const volScalarField::Internal crossDiffusion
    (
        max
        (
            sigmaD_/omega_()*(fvc::grad(k_) & fvc::grad(omega_))(),
            dimensionedScalar(dimless/sqr(dimTime), Zero)
        )
    );


    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Specific dissipation rate equation (W:Eq. 9)
    //  Variable changes: omega/k ~ 1/nut; P ~ G - 2/3 k divU
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        alpha()*rho()*gamma_*GbyNu
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(alpha()*rho()*beta0_*fBeta*omega_(), omega_)
      + alpha()*rho()*crossDiffusion
      + fvOptions(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation (W:Eq. 8)
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*betaStar_*omega_(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    // Effective specific dissipation rate used to compute nut (W:Eq. 6)
    omega_ = max(omega_, CLim_*sqrt(2.0/betaStar_)*mag(SBar));

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

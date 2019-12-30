/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "kEpsilonLRN.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
const Foam::Enum
<
    typename Foam::RASModels::kEpsilonLRN<BasicTurbulenceModel>::variantType
>
Foam::RASModels::kEpsilonLRN<BasicTurbulenceModel>::variantTypeNames
({
    { variantType::LAUNDER_SHARMA , "LaunderSharma" },
    { variantType::LAM_BREMHORST , "LamBremhorst" },
    { variantType::LIEN_LESCHZINER , "LienLeschziner" },
    { variantType::CHIEN , "Chien" },
    { variantType::MYONG_KASAGI , "MyongKasagi" },
    { variantType::ABE_KONDOH_NAGANO , "AbeKondohNagano" },
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilonLRN<BasicTurbulenceModel>::variantCoeffs()
{
    using ds = dimensioned<scalar>;

    Cmu_ = ds::getOrAddToDict("Cmu", this->coeffDict_, 0.09);
    C1_ = ds::getOrAddToDict("C1", this->coeffDict_, 1.44);
    C2_ = ds::getOrAddToDict("C2", this->coeffDict_, 1.92);
    C3_ = ds::getOrAddToDict("C3", this->coeffDict_, 0.0);
    sigmak_ = ds::getOrAddToDict("sigmak", this->coeffDict_, 1.0);
    sigmaEps_ = ds::getOrAddToDict("sigmaEps", this->coeffDict_, 1.3);

    switch (variant_)
    {
        case variantType::CHIEN:
        {
            C1_ = ds::getOrAddToDict("C1", this->coeffDict_, 1.35);
            C2_ = ds::getOrAddToDict("C2", this->coeffDict_, 1.8);
            break;
        }

        case variantType::MYONG_KASAGI:
        {
            C1_ = ds::getOrAddToDict("C1", this->coeffDict_, 1.4);
            C2_ = ds::getOrAddToDict("C2", this->coeffDict_, 1.8);
            sigmak_ = ds::getOrAddToDict("sigmak", this->coeffDict_, 1.4);
            break;
        }

        default:
            break;
    }
}


template<class BasicTurbulenceModel>
void kEpsilonLRN<BasicTurbulenceModel>::evaluateModelFuncs()
{
    switch (variant_)
    {
        case variantType::LAUNDER_SHARMA:
        {
            fMu_ = exp(-3.4/sqr(1.0 + sqr(k_)/(this->nu()*epsilon_)/50.0));
            f2_ = 1.0 - 0.3*exp(-min(sqr(sqr(k_)/(this->nu()*epsilon_)), 50.0));
            break;
        }

        case variantType::LAM_BREMHORST:
        {
            const tmp<volScalarField> Rt(sqr(k_)/(this->nu()*epsilon_));
            const tmp<volScalarField> Ry(sqrt(k_)*y_/this->nu());

            fMu_ = sqr(1.0 - exp(-0.0165*Ry))*(1.0 + 20.5/(Rt() + SMALL));
            f1_ = 1.0 + pow3(0.05/(fMu_ + SMALL));
            f2_ = 1.0 - exp(-sqr(Rt));
            break;
        }

        case variantType::LIEN_LESCHZINER:
        {
            const volScalarField yStar(sqrt(k_)*y_/this->nu());
            const tmp<volScalarField> Rt(sqr(k_)/(this->nu()*epsilon_));

            fMu_ = (1.0 - exp(-0.016*yStar))/(1.0 + SMALL - exp(-0.263*yStar));
            f2_ = 1.0 - 0.3*exp(-sqr(Rt));
            break;
        }

        case variantType::CHIEN:
        {
            const tmp<volScalarField> yStar(sqrt(k_)*y_/this->nu());
            const tmp<volScalarField> RT(sqr(k_)/(this->nu()*epsilon_));

            fMu_ = 1.0 - exp(-0.0115*yStar);
            f2_ = 1.0 - 2.0/9.0*exp(-sqr(RT/6.0));
            break;
        }

        case variantType::MYONG_KASAGI:
        {
            const tmp<volScalarField> yStar(sqrt(k_)*y_/this->nu());
            const tmp<volScalarField> RT(sqr(k_)/(this->nu()*epsilon_));

            fMu_ = (1.0 - exp(-yStar/70.0))*(1.0 + 3.45/(sqrt(RT) + SMALL));
            f2_ = (1.0 - 2.0/9.0*exp(-sqr(RT/6.0)))*sqr(1.0 - exp(-yStar/5.0));
            break;
        }

        case variantType::ABE_KONDOH_NAGANO:
        {
            const tmp<volScalarField> ReT(sqr(k_)/(this->nu()*epsilon_));
            const tmp<volScalarField> ReEpsilon
            (
                pow025(epsilon_)*y_/pow(this->nu(), 3.0/4.0)
            );

            fMu_ =
                (1.0 + 5.0/pow(ReT, 3.0/4.0))
               *exp(-sqr(ReT/200.0))
               *sqr(1.0 - exp(-ReEpsilon/14.0));
            
            f2_ = (1.0 - 0.3*exp(-sqr(ReT/6.5)))*sqr(1.0 - ReEpsilon/3.1);
            break;
        }
    }
}


template<class BasicTurbulenceModel>
void kEpsilonLRN<BasicTurbulenceModel>::evaluateBalancers()
{
    switch (variant_)
    {
        case variantType::LAUNDER_SHARMA:
        {
            D_ = 2.0*this->nu()*magSqr(fvc::grad(sqrt(k_)));
            E_ = 2.0*this->nu()*nut*fvc::magSqrGradGrad(U);
            break;
        }

        case variantType::LIEN_LESCHZINER:
        {
            const volScalarField yStar(sqrt(k_)*y_/this->nu());

            const tmp<volScalarField> le
            (
                0.41*y_*(1.0 + SMALL - exp(-0.263*yStar))
            );

            E_ =
                (C2_*pow(Cmu_, 0.75))
               *(f2_*sqrt(k_)*epsilon_/le)
               *exp(-0.00222*sqr(yStar));
            break;
        }

        case variantType::CHIEN:
        {
            const tmp<volScalarField> yStar(sqrt(k_)*y_/this->nu());

            D_ = 2.0*this->nu()*k_/(sqr(y_) + SMALL);

            E_ = (-2.0*this->nu()*epsilon_/(sqr(y_) + SMALL))*exp(-0.5*yStar);
            break;
        }

        default:
            break;
    }

}

template<class BasicTurbulenceModel>
void kEpsilonLRN<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*fMu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonLRN<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> kEpsilonLRN<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonLRN<BasicTurbulenceModel>::kEpsilonLRN
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
    variant_
    (
        variantTypeNames.getOrDefault
        (
            "variant",
            this->coeffDict_,
            variantType::LAUNDER_SHARMA
        )
    ),

    Cmu_(Zero),
    C1_(Zero),
    C2_(Zero),
    C3_(Zero),
    sigmak_(Zero),
    sigmaEps_(Zero),

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
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    fMu_
    (
        IOobject
        (
            "fMu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar("fMu", dimless, 1.0)
    ),
    f1_
    (
        IOobject
        (
            "f1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar("f1", dimless, 1.0)
    ),
    f2_
    (
        IOobject
        (
            "f2",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar("f2", dimless, 1.0)
    ),
    D_
    (
        IOobject
        (
            "D",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimLength)/pow4(dimTime), Zero)
    ),
    y_(wallDist::New(this->mesh_).y())
{
    variantCoeffs();

    if (mag(sigmak_.value()) < VSMALL || mag(sigmaEps_.value()) < VSMALL)
    {
        FatalErrorInFunction
            << "Non-zero values are required for the model constants:" << nl
            << "sigmak = " << sigmak_ << nl
            << "sigmaEps = " << sigmaEps_ << nl
            << exit(FatalError);
    }

    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    evaluateModelFuncs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonLRN<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kEpsilonLRN<BasicTurbulenceModel>::correct()
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

    const volScalarField divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volScalarField G
    (
        this->GName(), nut*(dev(twoSymm(tgradU())) && tgradU())
    );
    tgradU.clear();

    evaluateBalancers();


    // Turbulent kinetic energy dissipation rate equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha*rho*f1_*G*epsilon_/k_
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha*rho*divU, epsilon_)
      - fvm::Sp(C2_*f2_*alpha*rho*epsilon_/k_, epsilon_)
      + alpha*rho*E_
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G - fvm::SuSp(2.0/3.0*alpha*rho*divU, k_)
      - fvm::Sp(alpha*rho*(epsilon_ + D_)/k_, k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    evaluateModelFuncs();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

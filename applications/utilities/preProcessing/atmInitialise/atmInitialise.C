/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 CENER
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

Application
    atmInitialise

Group
    grpPreProcessingUtilities

Description
    \c atmInitialise is a preprocessing utility to initialise internal and/or
    patch fields of the following field types (each of which is optional) for
    an atmospheric boundary layer simulation.

    \table
      Field         | Description                                  | Units
      U             | Velocity                                     | [m/s]
      T             | Temperature                                  | [K]
      k             | Turbulent kinetic energy                     | [m^2/s^2]
      epsilon       | Turbulent kinetic energy dissipation rate    | [m^2/s^3]
      omega         | Specific dissipation rate                    | [1/s]
      p_\text{rgh}  | Perturbation pressure                        | [m^2/s^-2]
    \endtable

    \c atmInitialise utility is executed according to the settings in
    \f$ system/atmInitialiseDict \f$.

    where the common entries for each field mean:
    \table
      Property     | Description                         | Type | Req'd | Dflt
      initMethod   | Initialisation method               | word |  yes  | -
      interpMethod | Interpolation method for table input | word |  no  | linear
      initInternal | Flag to initialise internal field   | bool |  no   | true
      initPatch    | Flag to initialise patch field      | bool |  no   | true
    \endtable

    where the entries specific to \c U mean:
    \table
      Property   | Description                        | Type | Req'd | Dflt
      Uh         | Reference flow speed at height h   | scalar | yes | -
      h          | Reference height                   | scalar | yes | -
      z0         | Characteristic roughness length    | scalar | yes | -
      flowDir    | Streamwise flow direction          | vector | yes | -
      d          | Displacement height                | scalar | no  | 0
      kappa      | von Kármán constant                | scalar | no  | 0.40
    \endtable

    where the entries specific to \c T mean:
    \table
      Property   | Description                         | Type   | Req'd | Dflt
      Tground    | Potential temperature at the ground | scalar |  yes  | -
      Tground    | Potential temperature at the top of the capping inversion \\
                 | scalar | yes | -
      zInversion | Initial height of the capping inversion centre \\
                 | scalar | yes | -
      wInversion | Width of the capping inversion      | scalar |  yes  | -
      dTdz       | Rate of change of potential temperature with respect to \\
                   height above the capping inversion  | scalar |  yes  | -
    \endtable

    Options for the \c initMethod entry for \c U
    \table
      Property      | Description
      logNeutral    | -
      geostrophic   | -
      ekman         | -
      table         | -
    \endtable

    References:
    \verbatim
        Homogeneous two-dimensional ABL expressions (tag:RH):
            Richards, P. J., & Hoxey, R. P. (1993).
            Appropriate boundary conditions for computational wind
            engineering models using the k-ε turbulence model.
            In Computational Wind Engineering 1 (pp. 145-153).
            DOI:10.1016/B978-0-444-81688-7.50018-8

        Logarithmic profiles for unstable and stable stratifications (tag:E):
            Emeis, S. (2013).
            Wind Energy Meteorology: Atmospheric
            Physics for Wind Power Generation.
            Springer-Verlag Berlin Heidelberg
            DOI:10.1007/978-3-642-30523-8
    \endverbatim

Usage
    Minimal example by using \c atmInitialiseDict is as follows:
    \verbatim
        foamGetDict atmInitialiseDict
        atmInitialise
    \endverbatim

See also
    $FOAM_ETC/caseDicts/annotated/atmInitialiseDict

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "interpolateXY.H"
#include "interpolateSplineXY.H"
#include "Enum.H"
#include "wallDist.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum initMethod
{
    LOG_NEUTRAL,
//    LOG_UNSTABLE,
//    LOG_STABLE,
    GEOSTROPHIC,
//    EKMAN,
    TABLE
};

const Enum<initMethod> initMethodNames
{
    { initMethod::LOG_NEUTRAL, "logNeutral" },
    //{ initMethod::LOG_UNSTABLE, "logUnstable" },
    //{ initMethod::LOG_STABLE, "logStable" },
    { initMethod::GEOSTROPHIC, "geostrophic" },
    //{ initMethod::EKMAN, "ekman" },
    { initMethod::TABLE, "table" },
};

enum interpMethod
{
    LINEAR,
    CUBIC
};

const Enum<interpMethod> interpMethodNames
{
    { interpMethod::LINEAR, "linear" },
    { interpMethod::CUBIC, "cubic" },
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeField>
void calcTable(TypeField& fld, const enum interpMethod interpMode)
{
    typedef TypeField::cmptType cmptType;

    const List<Tuple2<scalar, cmptType>> inpTable
    (
        dict.lookup("inpTable")
    );
    scalarField zProfile(inpTable.size(), Zero);
    TypeField fldProfile(inpTable.size(), Zero);
    forAll(zProfile, i)
    {
        zProfile[i] = inpTable[i].first();
        fldProfile[i] = inpTable[i].second();
    }

    if (!inpTable.empty())
    {
        if (interpMode == interpMethod::LINEAR)
        {
            fld = interpolateXY(z, zProfile, fldProfile);
        }
        else if (interpMode == interpMethod::CUBIC)
        {
            fld = interpolateSplineXY(z, zProfile, fldProfile);
        }
    }
}


void initU(const dictionary& dict, const fvMesh& mesh)
{
    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const enum initMethod initMode(initMethodNames.get("initMethod", dict));
    const enum interpMethod interpMode
    (
        interpMethodNames.getOrDefault
        (
            "interpMethod",
            dict,
            interpMethod::LINEAR
        )
    );
    const scalar Uh(dict.getCheck<scalar>("Uh", scalarMinMax::ge(SMALL)));
    const scalar h(dict.getCheck<scalar>("h", scalarMinMax::ge(SMALL)));
    const scalar z0(dict.getCheck<scalar>("z0", scalarMinMax::ge(SMALL)));
    const vector flowDir
    (
        dict.getCheck<vector>
        (
            "flowDir",
            [&](const vector x){ return mag(x) > SMALL; }
        ).normalise()
    );
    const scalar kappa
    (
        dict.getCheckOrDefault<scalar>("kappa", 0.40, scalarMinMax::ge(SMALL))
    );
    const scalar d
    (
        dict.getCheckOrDefault<scalar>("d", 0.0, scalarMinMax::ge(0))
    );
    const bool initInternal(dict.getOrDefault<bool>("initInternal", true));
    const bool initPatch(dict.getOrDefault<bool>("initPatch", true));

    if (initInternal)
    {
        const scalarField& z = wallDist::New(mesh).y().primitiveField();
        vectorField& Ui = U.primitiveFieldRef();

        switch (initMode)
        {
            case LOG_NEUTRAL:
            {
                // (RH:Eqs. 23-24)
                const scalar Ustar = (kappa*Uh)/(Foam::log((h + z0)/z0));

                // (RH:Eq. 6)
                Ui = (Ustar/kappa*flowDir)*Foam::log((z + z0 - d)/z0);
                break;
            }

            case GEOSTROPHIC:
                Ui = Uh*flowDir;
                break;

            case TABLE:
                calcTable(U, interpMode);
                break;
            }

            default:
                break;
        }
    }

    if (initPatch)
    {
        U.correctBoundaryConditions();
    }

    if (initInternal || initPatch)
    {
        U.write();
    }
}


void initT(const dictionary& dict, const fvMesh& mesh)
{
    Info << "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const enum initMethod initMode(initMethodNames.get("initMethod", dict));
    const enum interpMethod interpMode
    (
        interpMethodNames.getOrDefault
        (
            "interpMethod",
            dict,
            interpMethod::CUBIC
        )
    );
    const bool initInternal(dict.getOrDefault<bool>("initInternal", true));
    const bool initPatch(dict.getOrDefault<bool>("initPatch", true));

    if (initInternal)
    {
        const scalarField& z = wallDist::New(mesh).y().primitiveField();
        scalarField& Ti = T.primitiveFieldRef();

        switch (initMode)
        {
            case GEOSTROPHIC:
            {
                const scalar Tground
                (
                    dict.getCheck<scalar>("Tground", scalarMinMax::ge(SMALL))
                );
                const scalar Ttop
                (
                    dict.getCheck<scalar>("Ttop", scalarMinMax::ge(SMALL))
                );
                const scalar zInversion
                (
                    dict.getCheck<scalar>("zInversion", scalarMinMax::ge(SMALL))
                );
                const scalar wInversion
                (
                    dict.getCheck<scalar>("wInversion", scalarMinMax::ge(SMALL))
                );
                const scalar dTdz
                (
                    dict.getCheck<scalar>("dTdz", scalarMinMax::ge(SMALL))
                );

                Ti = Tground;

                forAll(z, i)
                {
                    if
                    (
                        (z[i] >= zInversion - 0.5*wInversion) &&
                        (z[i] <= zInversion + 0.5*wInversion)
                    )
                    {
                        Ti[i] = Tground
                        +((Ttop - Tground)/wInversion)
                        *(z[i] - (zInversion - 0.5*wInversion));
                    }
                    else if (z[i] > zInversion + 0.5*wInversion)
                    {
                        Ti[i] =
                            Ttop + dTdz*(z[i] - (zInversion + 0.5*wInversion));
                    }
                }
                break;
            }

            case TABLE:
                calcTable(T, interpMode);
                break;

            default:
                break;
        }
    }

    if (initPatch)
    {
        T.correctBoundaryConditions();
    }

    if (initInternal || initPatch)
    {
        T.write();
    }
}


void initK(const dictionary& dict, const fvMesh& mesh)
{
    Info << "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const enum initMethod initMode(initMethodNames.get("initMethod", dict));
    const enum interpMethod interpMode
    (
        interpMethodNames.getOrDefault
        (
            "interpMethod",
            dict,
            interpMethod::CUBIC
        )
    );

    const scalar kAmb(dict.getCheck<scalar>("kAmb", scalarMinMax::ge(SMALL)));
    const bool initInternal(dict.getOrDefault<bool>("initInternal", true));
    const bool initPatch(dict.getOrDefault<bool>("initPatch", true));

    if (initInternal)
    {
        const scalarField& z = wallDist::New(mesh).y().primitiveField();
        scalarField& ki = k.primitiveFieldRef();

        switch (initMethod)
        {
            case LOG_NEUTRAL:
            {
                const scalar Uh(dict.getCheck<scalar>("Uh", scalarMinMax::ge(SMALL)));
                const scalar h(dict.getCheck<scalar>("h", scalarMinMax::ge(SMALL)));
                const scalar z0(dict.getCheck<scalar>("z0", scalarMinMax::ge(SMALL)));
                const scalar kappa
                (
                    dict.getCheck<scalar>("kappa", scalarMinMax::ge(SMALL))
                );
                const scalar Cmu
                (
                    dict.getCheck<scalar>("Cmu", scalarMinMax::ge(SMALL))
                );

                // (RH:Eqs. 23-24)
                const scalar Ustar = (kappa*Uh)/(Foam::log((h + z0)/z0));

                ki = max(sqr(Ustar)/Foam::sqrt(Cmu), kAmb);
                break;
            }

            case GEOSTROPHIC:
            {
                const scalar zInversion
                (
                    dict.getCheck<scalar>("zInversion", scalarMinMax::ge(SMALL))
                );

                forAll(z, i)
                {
                    if (z[i] <= zInversion)
                    {
                        // Taken from GABLS1 initial condition
                        ki[i] = max(0.4*pow3(1.0 - z[i]/zInversion), kAmb);
                    }
                    else
                    {
                        ki[i] = kAmb; 
                    }
                }
                break;
            }

            case TABLE:
                calcTable(k, interpMode);
                break;

            default:
                break;
        }
    }

    if (initPatch)
    {
        k.correctBoundaryConditions();
    }

    if (initInternal || initPatch)
    {
        k.write();
    }
}


void initEpsilon(const dictionary& dict, const fvMesh& mesh)
{
    Info << "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField k
    (
        IOobject
        (
            "k",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    const enum initMethod initMode(initMethodNames.get("initMethod", dict));
    const enum interpMethod interpMode
    (
        interpMethodNames.getOrDefault
        (
            "interpMethod",
            dict,
            interpMethod::CUBIC
        )
    );
    const scalar z0(dict.getCheck<scalar>("z0", scalarMinMax::ge(SMALL)));
    const scalar epsilonAmb
    (
        dict.getCheck<scalar>("epsilonAmb", scalarMinMax::ge(SMALL))
    );
    const bool initInternal(dict.getOrDefault<bool>("initInternal", true));
    const bool initPatch(dict.getOrDefault<bool>("initPatch", true));


    if (initInternal)
    {
        const scalarField& z = wallDist::New(mesh).y().primitiveField();
        scalarField& epsiloni = epsilon.primitiveFieldRef();
        scalarField& ki = k.primitiveFieldRef();

        switch (initMethod)
        {
            case LOG_NEUTRAL:
            {
                const scalar Uh(dict.getCheck<scalar>("Uh", scalarMinMax::ge(SMALL)));
                const scalar h(dict.getCheck<scalar>("h", scalarMinMax::ge(SMALL)));
                const scalar kappa(dict.getCheck<scalar>("kappa", scalarMinMax::ge(SMALL)));

                // (RH:Eqs. 23-24)
                const scalar Ustar = (kappa*Uh)/(Foam::log((h + z0)/z0));

                epsiloni = max(pow3(Ustar)/(kappa*(z + z0)), epsilonAmb);
                break;
            }

            case GEOSTROPHIC:
            {
                const scalar phim
                (
                    dict.getCheck<scalar>("phim", scalarMinMax::ge(SMALL))
                );
                const scalar Lmax
                (
                    dict.getCheck<scalar>("Lmax", scalarMinMax::ge(SMALL))
                );
                const scalar zInversion
                (
                    dict.getCheck<scalar>("zInversion", scalarMinMax::ge(SMALL))
                );
                const scalar Cmu
                (
                    dict.getCheck<scalar>("Cmu", scalarMinMax::ge(SMALL))
                );

                forAll(z, i)
                {
                    if (z[i] <= zTKEInversion)
                    {
                        // Taken from Blackadars parameterization
                        const scalar Lm = kappa*(z[i] + z0)/(phim + kappa*(z[i] + z0)/Lmax);
                        epsiloni[i] = max(pow(Cmu, 0.75)*pow(ki[i], 1.5)/Lm , epsilonAmb);
                    }
                    else
                    {
                        epsiloni[i] = epsilonAmb;
                    }
                }
            }

            case TABLE:
                calcTable(epsilon, interpMode);
                break;

            default:
                break;
    }

    if (initPatch)
    {
        epsilon.correctBoundaryConditions();
    }

    if (initInternal || initPatch)
    {
        epsilon.write();
    }
}


void initOmega(const dictionary& dict, const fvMesh& mesh)
{
    Info << "Reading field omega\n" << endl;
    volScalarField omega
    (
        IOobject
        (
            "omega",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    const enum initMethod initMode(initMethodNames.get("initMethod", dict));
    const enum interpMethod interpMode
    (
        interpMethodNames.getOrDefault
        (
            "interpMethod",
            dict,
            interpMethod::CUBIC
        )
    );

    const scalar omegaAmb
    (
        dict.getCheck<scalar>("omegaAmb", scalarMinMax::ge(SMALL))
    );
    const bool initInternal(dict.getOrDefault<bool>("initInternal", true));
    const bool initPatch(dict.getOrDefault<bool>("initPatch", true));


    if (initInternal)
    {
        const scalarField& z = wallDist::New(mesh).y().primitiveField();
        scalarField& omegai = omega.primitiveFieldRef();
        scalarField& ki = k.primitiveFieldRef();

        switch (initMethod)
        {
            case LOG_NEUTRAL:
            {
                const scalar Uh(dict.getCheck<scalar>("Uh", scalarMinMax::ge(SMALL)));
                const scalar h(dict.getCheck<scalar>("h", scalarMinMax::ge(SMALL)));
                const scalar z0(dict.getCheck<scalar>("z0", scalarMinMax::ge(SMALL)));
                const scalar kappa
                (
                    dict.getCheck<scalar>("kappa", scalarMinMax::ge(SMALL))
                );
                // (RH:Eqs. 23-24)
                const scalar Ustar = (kappa*Uh)/(Foam::log((h + z0)/z0));

                const scalarField epsilon(pow3(Ustar)/(kappa*(z + z0)));

                omegai = max(epsilon/ki, omegaAmb);
                break;
            }

            case GEOSTROPHIC:
            {
                const scalar phim
                (
                    dict.getCheck<scalar>("phim", scalarMinMax::ge(SMALL))
                );
                const scalar Lmax
                (
                    dict.getCheck<scalar>("Lmax", scalarMinMax::ge(SMALL))
                );
                const scalar zInversion
                (
                    dict.getCheck<scalar>("zInversion", scalarMinMax::ge(SMALL))
                );
                const scalar Cmu
                (
                    dict.getCheck<scalar>("Cmu", scalarMinMax::ge(SMALL))
                );

                forAll(z, i)
                {
                    if (z[i] <= zInversion)
                    {
                        // Taken from Blackadars parameterization
                        const scalar Lm = kappa*(z[i] + z0)/(phim + kappa*(z[i] + z0)/Lmax);
                        const scalar epsilon = pow(Cmu, 0.75)*pow(ki[i], 1.5)/Lm;
                        omegai[i] = max(epsilon/ki[i], omegaAmb);
                    }
                    else
                    {
                        omegai[i] = omegaAmb;
                    }
                }
            }

            case TABLE:
                calcTable(omega, interpMode);
                break;

            default:
                break;
    }

    if (initPatch)
    {
        omega.correctBoundaryConditions();
    }

    if (initInternal || initPatch)
    {
        omega.write();
    }
}


void initPrgh(const dictionary& dict)
{
    Info << "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    if (initInternal)
    {
        const scalar pRef(dict.get<scalar>("pRef"));
        p_rgh = pRef;
    }

    if (initPatch)
    {
        p_rgh.correctBoundaryConditions();
    }

    if (initInternal || initPatch)
    {
        p_rgh.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Initialise internal fields for an atmospheric boundary layer flow"
        " simulation by using system/atmInitialiseDict dictionary"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "atmInitialiseDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (dict.found("U"))
    {
        const dictionary subDict(dict.subDict("U"));
        initU(subDict, mesh);
    }

    if (dict.found("T"))
    {
        const dictionary subDict(dict.subDict("T"));
        initT(subDict, mesh);
    }

    if (dict.found("k"))
    {
        const dictionary subDict(dict.subDict("k"));
        initK(subDict);
    }

    if (dict.found("epsilon"))
    {
        const dictionary subDict(dict.subDict("epsilon"));
        initEpsilon(subDict);
    }

    if (dict.found("omega"))
    {
        const dictionary subDict(dict.subDict("omega"));
        initOmega(subDict);
    }

    if (dict.found("p_rgh"))
    {
        const dictionary subDict(dict.subDict("p_rgh"));
        initPrgh(subDict);
    }


    runTime.printExecutionTime(Info);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

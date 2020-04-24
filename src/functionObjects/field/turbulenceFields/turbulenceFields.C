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

\*---------------------------------------------------------------------------*/

#include "turbulenceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(turbulenceFields, 0);
    addToRunTimeSelectionTable(functionObject, turbulenceFields, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::turbulenceFields::compressibleField
>
Foam::functionObjects::turbulenceFields::compressibleFieldNames_
({
    { compressibleField::cfK, "k" },
    { compressibleField::cfEpsilon, "epsilon" },
    { compressibleField::cfOmega, "omega" },
    { compressibleField::cfNuTilda, "nuTilda" },
    { compressibleField::cfMut, "mut" },
    { compressibleField::cfMuEff, "muEff" },
    { compressibleField::cfAlphat, "alphat" },
    { compressibleField::cfAlphaEff, "alphaEff" },
    { compressibleField::cfR, "R" },
    { compressibleField::cfDevRhoReff, "devRhoReff" },
    { compressibleField::cfL, "L" },
    { compressibleField::cfI, "I" },
});


const Foam::Enum
<
    Foam::functionObjects::turbulenceFields::incompressibleField
>
Foam::functionObjects::turbulenceFields::incompressibleFieldNames_
({
    { incompressibleField::ifK, "k" },
    { incompressibleField::ifEpsilon, "epsilon" },
    { incompressibleField::ifOmega, "omega" },
    { incompressibleField::ifNuTilda, "nuTilda" },
    { incompressibleField::ifNut, "nut" },
    { incompressibleField::ifNuEff, "nuEff" },
    { incompressibleField::ifR, "R" },
    { incompressibleField::ifDevReff, "devReff" },
    { incompressibleField::ifL, "L" },
    { incompressibleField::ifI, "I" },
});


const Foam::word Foam::functionObjects::turbulenceFields::modelName_
(
    Foam::turbulenceModel::propertiesName
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFields::compressible()
{
    if (obr_.foundObject<compressible::turbulenceModel>(modelName_))
    {
        return true;
    }
    else if (obr_.foundObject<incompressible::turbulenceModel>(modelName_))
    {
        return false;
    }

    FatalErrorInFunction
        << "Turbulence model not found in database, deactivating"
        << exit(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::turbulenceFields::turbulenceFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_()
{
    const auto* turbPtr = obr_.findObject<Foam::turbulenceModel>(modelName_);

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (!tepsilon.isTmp())
    {
        modelBase_ = modelBase::EPSILON;
    }
    else if (!tomega.isTmp())
    {
        modelBase_ = modelBase::OMEGA;
    }
    else
    {
        modelBase_ = modelBase::OTHER;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::turbulenceFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    if (dict.found("field"))
    {
        fieldSet_.insert(dict.get<word>("field"));
    }
    else
    {
        fieldSet_.insert(dict.get<wordList>("fields"));
    }

    Info<< type() << " " << name() << ": ";
    if (fieldSet_.size())
    {
        Info<< "storing fields:" << nl;
        for (const word& f : fieldSet_)
        {
            Info<< "    " << modelName_ << ':' << f << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no fields requested to be stored" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::turbulenceFields::execute()
{
    bool comp = compressible();

    if (comp)
    {
        const compressible::turbulenceModel& model =
            obr_.lookupObject<compressible::turbulenceModel>(modelName_);

        for (const word& f : fieldSet_)
        {
            switch (compressibleFieldNames_[f])
            {
                case cfK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case cfEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case cfOmega:
                {
                    processField<scalar>(f, model.omega());
                    break;
                }
                case cfNuTilda:
                {
                    processField<scalar>(f, nuTilda(model));
                    break;
                }
                case cfMut:
                {
                    processField<scalar>(f, model.mut());
                    break;
                }
                case cfMuEff:
                {
                    processField<scalar>(f, model.muEff());
                    break;
                }
                case cfAlphat:
                {
                    processField<scalar>(f, model.alphat());
                    break;
                }
                case cfAlphaEff:
                {
                    processField<scalar>(f, model.alphaEff());
                    break;
                }
                case cfR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case cfDevRhoReff:
                {
                    processField<symmTensor>(f, model.devRhoReff());
                    break;
                }
                case cfI:
                {
                    processField<scalar>(f, I(model));
                    break;
                }
                case cfL:
                {
                    processField<scalar>(f, L(model));
                    break;
                }
                case cfT:
                {
                    processField<scalar>(f, T(model));
                    break;
                }
                case cfU:
                {
                    processField<scalar>(f, U(model));
                    break;
                }
                case cfLt:
                {
                    processField<scalar>(f, Lt(model));
                    break;
                }
                case cfTt:
                {
                    processField<scalar>(f, Tt(model));
                    break;
                }
                case cfUt:
                {
                    processField<scalar>(f, Ut(model));
                    break;
                }
                case cfLk:
                {
                    processField<scalar>(f, Lk(model));
                    break;
                }
                case cfTk:
                {
                    processField<scalar>(f, Tk(model));
                    break;
                }
                case cfUk:
                {
                    processField<scalar>(f, Uk(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }
    else
    {
        const incompressible::turbulenceModel& model =
            obr_.lookupObject<incompressible::turbulenceModel>(modelName_);

        for (const word& f : fieldSet_)
        {
            switch (incompressibleFieldNames_[f])
            {
                case ifK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case ifEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case ifOmega:
                {
                    processField<scalar>(f, model.omega());
                    break;
                }
                case ifNuTilda:
                {
                    processField<scalar>(f, nuTilda(model));
                    break;
                }
                case ifNut:
                {
                    processField<scalar>(f, model.nut());
                    break;
                }
                case ifNuEff:
                {
                    processField<scalar>(f, model.nuEff());
                    break;
                }
                case ifR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case ifDevReff:
                {
                    processField<symmTensor>(f, model.devReff());
                    break;
                }
                case ifI:
                {
                    processField<scalar>(f, I(model));
                    break;
                }
                case ifL:
                {
                    processField<scalar>(f, L(model));
                    break;
                }
                case ifT:
                {
                    processField<scalar>(f, T(model));
                    break;
                }
                case ifU:
                {
                    processField<scalar>(f, U(model));
                    break;
                }
                case ifLt:
                {
                    processField<scalar>(f, Lt(model));
                    break;
                }
                case ifTt:
                {
                    processField<scalar>(f, Tt(model));
                    break;
                }
                case ifUt:
                {
                    processField<scalar>(f, Ut(model));
                    break;
                }
                case ifLk:
                {
                    processField<scalar>(f, Lk(model));
                    break;
                }
                case ifTk:
                {
                    processField<scalar>(f, Tk(model));
                    break;
                }
                case ifUk:
                {
                    processField<scalar>(f, Uk(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::turbulenceFields::write()
{
    for (const word& f : fieldSet_)
    {
        const word fieldName = modelName_ + ':' + f;
        writeObject(fieldName);
    }

    return true;
}


// ************************************************************************* //

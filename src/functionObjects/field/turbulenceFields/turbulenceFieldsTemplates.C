/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::turbulenceFields::processField
(
    const word& fieldName,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvalue
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    const word scopedName = modelName + ':' + fieldName;

    FieldType* fldPtr = obr_.getObjectPtr<FieldType>(scopedName);

    if (fldPtr)
    {
        (*fldPtr) == tvalue();
    }
    else if (obr_.found(scopedName))
    {
        WarningInFunction
            << "Cannot store turbulence field " << scopedName
            << " since an object with that name already exists"
            << nl << endl;
    }
    else
    {
        obr_.store
        (
            new FieldType
            (
                IOobject
                (
                    scopedName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                tvalue
            )
        );
    }
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::nuTilda
(
    const Model& model
) const
{
    return tmp<volScalarField>::New
    (
        "nuTilda.tmp",
        model.k()/omega(model)
    );
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::I
(
    const Model& model
) const
{
    // Assume k is available
    const volScalarField uPrime(sqrt((2.0/3.0)*model.k()));
    const dimensionedScalar U0("U0", dimVelocity, SMALL);

    return tmp<volScalarField>::New
    (
        "I.tmp",
        uPrime/max(max(uPrime, mag(model.U())), U0)
    );
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::L
(
    const Model& model
) const
{
    const dictionary& modelDict = model.coeffDict();

    if (modelBase_ == modelBase::EPSILON)
    {
        tmp<volScalarField> tepsilon = model.epsilon();
        tmp<volScalarField> tk = model.k();
        const dimensionedScalar Cmu(read("Cmu", modelDict));
        const dimensionedScalar eps0("e0", sqr(dimLength)/pow3(dimTime), SMALL);

        return tmp<volScalarField>::New
        (
            "L.tmp",
            pow(Cmu, 0.75)*pow(tk, 1.5)/(tepsilon + eps0)
        );
    }
    else if (modelBase_ == modelBase::OMEGA)
    {
        tmp<volScalarField> tomega = model.omega();
        tmp<volScalarField> tk = model.k();
        const dimensionedScalar Cmu(read("Cmu", modelDict));
        const dimensionedScalar omg0("o0", dimless/dimTime, SMALL);

        return tmp<volScalarField>::New
        (
            "L.tmp",
            Foam::sqrt(tk)/(pow025(Cmu)*(tomega + omg0))
        );
    }

    // pending: what would happen if model=SpalartAllmaras/WrayAgarwal
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::T
(
    const Model& model
) const
{
    const dictionary& modelDict = model.coeffDict();

    if (modelBase_ == modelBase::EPSILON)
    {
        tmp<volScalarField> tepsilon = model.epsilon();
        tmp<volScalarField> tk = model.k();
        const dimensionedScalar Cmu(read("Cmu", modelDict));
        const dimensionedScalar eps0("e0", sqr(dimLength)/pow3(dimTime), SMALL);

        return tmp<volScalarField>::New
        (
            "T.tmp",
            Cmu*tk/(tepsilon + eps0)
        );
    }
    else if (modelBase_ == modelBase::OMEGA)
    {
        tmp<volScalarField> tomega = model.omega();
        const dimensionedScalar Cmu(read("Cmu", modelDict));
        const dimensionedScalar omg0("o0", dimless/dimTime, SMALL);

        return tmp<volScalarField>::New
        (
            "T.tmp",
            Cmu/(tomega + omg0)
        );
    }
}


template<class Model>
Foam::tmp<Foam::volScalarField>
Foam::functionObjects::turbulenceFields::U
(
    const Model& model
) const
{
    return tmp<volScalarField>::New
    (
        "U.tmp",
        Foam::sqrt(model.k())
    );
}


// ************************************************************************* //

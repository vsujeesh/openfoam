/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "mappedMixedFvPatchField.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    mappedPatchBase(p.patch()),
    mappedPatchFieldBase<Type>(*this, *this),
    weightFieldName_(word::null)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict),
    mappedPatchBase(p.patch(), dict),
    mappedPatchFieldBase<Type>(*this, *this, dict),
    weightFieldName_(dict.get<word>("weightField"))
{
    mixedFvPatchField<Type>::operator=
    (
        Field<Type>("value", dict, p.size())
    );

    if (dict.found("refValue"))
    {
        // Full restart
        this->refValue() = Field<Type>("refValue", dict, p.size());
        this->refGrad() = Field<Type>("refGradient", dict, p.size());
        this->valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        this->refValue() = *this;
        this->refGrad() = Zero;
        this->valueFraction() = 1.0;
    }
}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const mappedMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    mappedPatchBase(p.patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


//template<class Type>
//Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
//(
//    const fvPatch& p,
//    const DimensionedField<Type, volMesh>& iF,
//
//    // mappedPatchBase
//    const word& sampleRegion,
//    const sampleMode sampleMode,
//    const word& samplePatch,
//    const scalar distance,
//
//    // My settings
//    const word& fieldName,
//    const bool setAverage,
//    const Type average,
//    const word& interpolationScheme
//)
//:
//    mixedFvPatchField<Type>(p, iF),
//    mappedPatchBase
//    (
//        p.patch(),
//        sampleRegion,
//        sampleMode,
//        samplePatch,
//        distance
//    ),
//    mappedPatchFieldBase<Type>
//    (
//        *this,
//        *this,
//        fieldName,
//        setAverage,
//        average,
//        interpolationScheme
//    )
//    weightFieldName_(ptf.weightFieldName_)
//{}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const mappedMixedFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


template<class Type>
Foam::mappedMixedFvPatchField<Type>::mappedMixedFvPatchField
(
    const mappedMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    mappedPatchBase(ptf.patch().patch(), ptf),
    mappedPatchFieldBase<Type>(*this, *this, ptf),
    weightFieldName_(ptf.weightFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedMixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<Type>::autoMap(m);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedMixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<Type>::rmap(ptf, addr);
    mappedPatchBase::clearOut();
}


template<class Type>
void Foam::mappedMixedFvPatchField<Type>::updateCoeffs()
{
    //typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (this->updated())
    {
        return;
    }
    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = *this;

    //const bool isSampleWorld(UPstream::myWorld() == mapper_.sampleWorld());
    const bool isSampleWorld = true;

    // Swap to obtain full local values of neighbour internal field
    Field<Type> nbrIntFld;
    scalarField nbrKDelta;

    if (isSampleWorld)
    {
        const auto& nbrMesh = refCast<const fvMesh>(this->mapper_.sampleMesh());
        const label nbrPatchID = mpp.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];
        const auto& nbrField = this->sampleField();

        nbrIntFld = nbrField.boundaryField()[nbrPatchID].patchInternalField();

        //Pout<< "Sampled patch:" << nbrPatch.name() << nl
        //    << "    patchValue:" << *this << nl
        //    << "    patchIntValue:" << nbrIntFld
        //    << endl;
        // Weightfield is volScalarField
        const auto& nbrWeightField =
            //this->sampleField<scalar>(weightFieldName_);
            nbrMesh.template lookupObject<volScalarField>(weightFieldName_);

        nbrKDelta =
            nbrWeightField.boundaryField()[nbrPatchID].patchInternalField()
           *nbrPatch.deltaCoeffs();
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    mpp.distribute(nbrIntFld);
    mpp.distribute(nbrKDelta);

DebugVar(nbrIntFld);
DebugVar(nbrKDelta);


    // Restore tag
    UPstream::msgType() = oldTag;

    const auto& myWeightField =
    this->db().objectRegistry::lookupObject<volScalarField>
    (
        weightFieldName_
    );

    const label patchi = this->patch().index();

    tmp<scalarField> myKDelta =
        myWeightField.boundaryField()[patchi].patchInternalField()
       *this->patch().deltaCoeffs();


    // Both sides agree on
    // - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
    // - gradient    : (temperature-fld)*delta
    // We've got a degree of freedom in how to implement this in a mixed bc.
    // (what gradient, what fixedValue and mixing coefficient)
    // Two reasonable choices:
    // 1. specify above temperature on one side (preferentially the high side)
    //    and above gradient on the other. So this will switch between pure
    //    fixedvalue and pure fixedgradient
    // 2. specify gradient and temperature such that the equations are the
    //    same on both sides. This leads to the choice of
    //    - refGradient = zero gradient
    //    - refValue = neighbour value
    //    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    this->refValue() = nbrIntFld;
    this->refGrad() = Zero;
    this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta());

    mixedFvPatchScalarField::updateCoeffs();

    //if (debug)
    //{
    //    scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());
    //
    //    Info<< patch().boundaryMesh().mesh().name() << ':'
    //        << patch().name() << ':'
    //        << this->internalField().name() << " <- "
    //        << nbrMesh.name() << ':'
    //        << nbrPatch.name() << ':'
    //        << this->internalField().name() << " :"
    //        << " heat transfer rate:" << Q
    //        << " walltemperature "
    //        << " min:" << gMin(*this)
    //        << " max:" << gMax(*this)
    //        << " avg:" << gAverage(*this)
    //        << endl;
    //}

    mixedFvPatchField<Type>::updateCoeffs();

    Pout<< "** now this value:" << static_cast<const Field<Type>&>(*this)
        << endl;
}


template<class Type>
void Foam::mappedMixedFvPatchField<Type>::write(Ostream& os) const
{
    mappedPatchBase::write(os);
    mappedPatchFieldBase<Type>::write(os);
    os.writeEntry("weightField", weightFieldName_);
    mixedFvPatchField<Type>::write(os);
}


// ************************************************************************* //

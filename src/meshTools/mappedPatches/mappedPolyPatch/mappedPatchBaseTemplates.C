/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

template<class Type>
void Foam::mappedPatchBase::distribute(List<Type>& lst) const
{
    const label oldComm(Pstream::warnComm);
    Pstream::warnComm = map().comm();
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm = Pstream::worldComm;
            Pstream::worldComm = comm_;
            lst = AMI().interpolateToSource(Field<Type>(std::move(lst)));
            Pstream::worldComm = oldWorldComm;
            break;
        }
        default:
        {
            map().distribute(lst);
        }
    }
    Pstream::warnComm = oldComm;
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::distribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    const label oldComm(Pstream::warnComm);
    Pstream::warnComm = comm_;
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm = Pstream::worldComm;
            Pstream::worldComm = comm_;
            lst = AMI().interpolateToSource(Field<Type>(std::move(lst)), cop);
            Pstream::worldComm = oldWorldComm;
            break;
        }
        default:
        {
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                map().schedule(),
                map().constructSize(),
                map().subMap(),
                false,
                map().constructMap(),
                false,
                lst,
                Type(Zero),
                cop,
                flipOp(),
                UPstream::msgType(),
                comm_
            );
        }
    }
    Pstream::warnComm = oldComm;
}


template<class Type>
void Foam::mappedPatchBase::reverseDistribute(List<Type>& lst) const
{
    const label oldComm(Pstream::warnComm);
    Pstream::warnComm = map().comm();
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm = Pstream::worldComm;
            Pstream::worldComm = comm_;
            lst = AMI().interpolateToTarget(Field<Type>(std::move(lst)));
            Pstream::worldComm = oldWorldComm;
            break;
        }
        default:
        {
            map().reverseDistribute(sampleSize(), lst);
            break;
        }
    }
    Pstream::warnComm = oldComm;
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::reverseDistribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    const label oldComm(Pstream::warnComm);
    Pstream::warnComm = map().comm();
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            const label oldWorldComm = Pstream::worldComm;
            Pstream::worldComm = comm_;
            lst = AMI().interpolateToTarget(Field<Type>(std::move(lst)), cop);
            Pstream::worldComm = oldWorldComm;
            break;
        }
        default:
        {
            label cSize = sampleSize();
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                map().schedule(),
                cSize,
                map().constructMap(),
                false,
                map().subMap(),
                false,
                lst,
                Type(Zero),
                cop,
                flipOp(),
                UPstream::msgType(),
                comm_
            );
            break;
        }
    }
    Pstream::warnComm = oldComm;
}


// ************************************************************************* //

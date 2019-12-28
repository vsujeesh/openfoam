/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "BinSum.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IndexType, class List, class CombineOp>
Foam::BinSum<IndexType, List, CombineOp>::BinSum
(
    const IndexType min,
    const IndexType max,
    const IndexType delta
)
:
    List
    (
        (min < max && delta <= ceil(max - min) && SMALL < delta)
      ? List(ceil((max - min)/delta), Zero)
      : throw std::domain_error("Inconsistent min/max/delta combination.")
    ),
    min_(min),
    max_(max),
    delta_(delta),
    lowSum_(Zero),
    highSum_(Zero)
{}


template<class IndexType, class List, class CombineOp>
Foam::BinSum<IndexType, List, CombineOp>::BinSum
(
    const IndexType min,
    const IndexType max,
    const IndexType delta,
    const List& mainVals,
    const CombineOp& cop
)
:
    BinSum<IndexType, List>(min, max, delta, mainVals, mainVals, cop)
{}


template<class IndexType, class List, class CombineOp>
Foam::BinSum<IndexType, List, CombineOp>::BinSum
(
    const IndexType min,
    const IndexType max,
    const IndexType delta,
    const UList<IndexType>& auxVals,
    const List& mainVals,
    const CombineOp& cop
)
:
    List
    (
        (min < max && delta <= ceil(max - min) && SMALL < delta)
      ? List(ceil((max - min)/delta), Zero)
      : throw std::domain_error("Inconsistent min/max/delta combination.")
    ),
    min_(min),
    max_(max),
    delta_(delta),
    lowSum_(Zero),
    highSum_(Zero)
{
    sumIntoBins(auxVals, mainVals, cop);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class IndexType, class List, class CombineOp>
void Foam::BinSum<IndexType, List, CombineOp>::sumIntoBins
(
    const IndexType& auxVal,
    const typename List::const_reference mainVal,
    const CombineOp& cop
)
{
    if (auxVal < min_)
    {
        cop(lowSum_, mainVal);
    }
    else if (auxVal >= max_)
    {
        cop(highSum_, mainVal);
    }
    else
    {
        label index = (auxVal - min_)/delta_;
        cop(this->operator[](index), mainVal);
    }
}


template<class IndexType, class List, class CombineOp>
void Foam::BinSum<IndexType, List, CombineOp>::sumIntoBins
(
    const UList<IndexType>& auxVals,
    const List& mainVals,
    const CombineOp& cop
)
{
    forAll(auxVals, i)
    {
        sumIntoBins(auxVals[i], mainVals[i], cop);
    }
}

template<class IndexType, class List, class CombineOp>
void Foam::BinSum<IndexType, List, CombineOp>::sumIntoBins
(
    const List& mainVals,
    const CombineOp& cop
)
{
    sumIntoBins(mainVals, mainVals, cop);
}

// ************************************************************************* //

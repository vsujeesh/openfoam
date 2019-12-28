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

Class
    Foam::BinSum

Description
    List container that partitions a given dataset into equal-width
    bins of a specified data window, and returns the sum of dataset
    elements in each bin.

    The container can also accept an auxillary dataset and a main dataset,
    wherein the specified data window can be applied onto the auxillary
    dataset, and the corresponding-in-order elements of the main dataset
    can be summed into the bins.

SourceFiles
    BinSum.C

SeeAlso
    Test-BinSum.C

\*---------------------------------------------------------------------------*/

#ifndef BinSum_H
#define BinSum_H

#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class BinSum Declaration
\*---------------------------------------------------------------------------*/

template
<
    class IndexType,
    class List,
    class CombineOp = plusEqOp<typename List::value_type>
>
class BinSum
:
    public List
{
    static_assert
    (
        std::is_same<IndexType, label>::value
     || std::is_same<IndexType, scalar>::value,
        "BinSum operates only in label/scalar arithmetic."
    );

    static_assert
    (
        std::is_same<CombineOp, plusEqOp<typename List::value_type>>::value,
        "BinSum operates only for sum operations."
    );


    // Private Data

        //- The lowest-edge value of the specified data window
        const IndexType min_;

        //- The highest-edge value of the specified data window
        const IndexType max_;

        //- Width of a bin
        //  All bins have the same bin-width
        //  Number of bins = window width / delta = (max - min) / delta
        const IndexType delta_;

        //- Sum of all dataset elements smaller than
        //- 'min_' outside of the data window
        typename List::value_type lowSum_;

        //- Sum of all dataset elements larger or equal
        //- than 'max_' outside of the data window
        typename List::value_type highSum_;


public:

    // Constructors

        //- Construct given min, max, and delta
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta
        );

        //- Construct given min, max, delta, and dataset
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta,
            const List& mainVals,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );

        //- Construct given min, max, delta, auxillary dataset, and main dataset
        BinSum
        (
            const IndexType min,
            const IndexType max,
            const IndexType delta,
            const UList<IndexType>& auxVals,
            const List& mainVals,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );


    //- Destructor
    ~BinSum() = default;


    // Member Functions

        //- Sum a dataset element into the corresponding bin
        //- with an auxillary dataset element
        void sumIntoBins
        (
            const IndexType& auxVal,
            const typename List::const_reference mainVal,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );

        //- Sum dataset elements into corresponding bins
        //- with an auxillary dataset
        void sumIntoBins
        (
            const UList<IndexType>& auxVals,
            const List& mainVals,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );

        //- Sum dataset elements into corresponding bins
        void sumIntoBins
        (
            const List& mainVals,
            const CombineOp& cop = plusEqOp<typename List::value_type>()
        );


        // Access

            //- Return delta_
            inline IndexType delta() const
            {
                return delta_;
            }

            //- Return lowSum_
            inline const IndexType& lowSum() const
            {
                return lowSum_;
            }

            //- Return highSum_
            inline const IndexType& highSum() const
            {
                return highSum_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BinSum.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

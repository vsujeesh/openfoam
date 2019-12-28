/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2013 OpenFOAM Foundation
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

Application
    Test-BinSum

Description
    Test application for 'BinSum' container

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "BinSum.H"
#include "IOstreams.H"
#include "Random.H"
#include "scalarField.H"
#include "complex.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class IndexType, class List>
void test_constructors
(
    const IndexType min,
    const IndexType max,
    const IndexType delta,
    const List& dataset
)
{
    Info<< "## Test constructors ##" << endl;
    {
        Info<< "# Construct given min, max, and delta" << endl;
        BinSum<IndexType, Field<IndexType>> container(min, max, delta);
        Info<< "container:" << tab << container << endl;
    }

    {
        Info<< "# Construct given min, max, delta and dataset" << endl;
        BinSum<IndexType, Field<IndexType>> container(min, max, delta, dataset);
        Info<< "container:" << tab << container << endl;
    }

    {
        Info<< "# Construct given min, max, delta, auxillary dataset, "
            << "and main dataset" << endl;
        const List auxDataset(dataset);
        BinSum<IndexType, Field<IndexType>> container
        (
            min, max, delta, auxDataset, dataset
        );
        Info<< "container:" << tab << container << endl;
    }
}


template<class IndexType, class List>
void test_member_funcs
(
    BinSum<IndexType, List>& window,
    const List& dataset
)
{
    Info<< "## Test member functions ##" << endl;

    BinSum<IndexType, List> binsum(window);
    binsum.sumIntoBins(dataset);

    const List auxDataset(dataset);
    BinSum<IndexType, List> binsumWithAuxDataset(window);
    binsumWithAuxDataset.sumIntoBins(auxDataset, dataset);

    BinSum<IndexType, List> bincounts(window);
    BinSum<IndexType, List> binsumElementwise(window);
    for (const auto& sample : dataset)
    {
        bincounts.sumIntoBins(sample, 1);
        binsumElementwise.sumIntoBins(sample, sample);
    }


    Info<< "binsum                     : " << binsum << endl;
    Info<< "binsum with aux dataset    : " << binsumWithAuxDataset << endl;
    Info<< "binsum elementwise         : " << binsumElementwise << endl;
    Info<< "bin counts                 : " << bincounts << endl;
    Info<< "bin arithmetic average     : " << binsum/bincounts << endl;
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    const scalar min = 0;
    const scalar max = 1;
    const scalar delta = 0.1;

    Random rndGen(0);
    scalarField dataset(10000000);
    for (auto& sample: dataset)
    {
        sample = rndGen.sample01<scalar>();
    }

    BinSum<scalar, scalarField> window(min, max, delta);


    test_constructors(min, max, delta, dataset);

    test_member_funcs(window, dataset);


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

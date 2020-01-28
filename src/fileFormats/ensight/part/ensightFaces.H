/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::ensightFaces

Description
    Sorting/classification of faces (2D) into corresponding ensight types.

    Some caution may be required when handling face addressing into fields.
    If it is a boundaryField, the local face ids should be used.

\*---------------------------------------------------------------------------*/

#ifndef ensightFaces_H
#define ensightFaces_H

#include "ensightPart.H"
#include "boolList.H"
#include "faceList.H"
#include "FixedList.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightFaces
:
    public ensightPart
{
public:

    // Public Data

        //- Supported ensight 'Face' element types.
        //  Must be zero-based since they are also for internal bookkeeping.
        enum elemType
        {
            TRIA3 = 0,  //!< "tria3"
            QUAD4,      //!< "quad4"
            NSIDED      //!< "nsided"
        };

        //- Number of 'Face' element types (3)
        static constexpr int nTypes = 3;

        //- The ensight 'Face' element type names
        static const char* elemNames[nTypes];


    // Static Functions

        //- The ensight element name for the specified 'Face' type
        static inline const char* key(const elemType etype);


private:

    // Private Data

        //- Linear list of face-flips
        boolList flipMap_;

        //- Begin/end offsets for address/flips of each element type
        FixedList<label, nTypes+1> offsets_;

        //- List of global sizes for each element type.
        //  Used temporarily for local sizes when building the element lists.
        FixedList<label, nTypes> sizes_;


    // Private Member Functions

        //- Low-level internal addition routine
        inline void add(const elemType etype, label id, bool flip=false);

        //- Use temporarily stored sizes to redimension the element lists
        void resizeAll();


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("ensightFaces");


    // Generated Methods

        //- Copy construct
        ensightFaces(const ensightFaces&) = default;

        //- No copy assignment
        void operator=(const ensightFaces&) = delete;

        //- Move construct
        ensightFaces(ensightFaces&&) = default;

        //- Move assignment
        ensightFaces& operator=(ensightFaces&&) = default;

        //- Destructor
        virtual ~ensightFaces() = default;


    // Constructors

        //- Default construct, with part index 0
        ensightFaces();

        //- Default construct, with description/partName
        explicit ensightFaces(const string& description);

        //- Default construct, with specified part index
        explicit ensightFaces(const label partIndex);


    // Member Functions

    // Access

        //- Processor-local size of all elements.
        using ensightPart::size;

        //- Processor-local size of the specified element type.
        inline label size(const elemType etype) const;

        //- Processor-local offset/size of element type.
        inline labelRange range(const elemType etype) const;

        //- The global size of all element types.
        //  This value is only meaningful after a reduce operation.
        label total() const;

        //- The global size of the specified element type.
        //  This value is only meaningful after a reduce operation.
        inline label total(const elemType etype) const;

        //- The global sizes for each element type.
        //  This value is only meaningful after a reduce operation.
        inline const FixedList<label, nTypes>& totals() const;

        //- Processor-local sizes per element type.
        FixedList<label, nTypes> sizes() const;

        //- Processor-local face ids of all elements
        inline const labelList& faceIds() const;

        //- Processor-local face ids of the specified element type
        inline const labelUList faceIds(const elemType etype) const;

        //- Processor-local flip-map of all elements
        inline const boolList& flipMap() const;

        //- Processor-local (face) ids relative to start()
        inline labelList localIds() const;

        //- Processor-local (face) ids relative to start()
        inline labelList localIds(const elemType etype) const;


    // Edit

        //- Classify the face types and set the element lists.
        void classify(const UList<face>& faces);

        //- Classify face types (for a sublist) and set element lists.
        void classify(const UList<face>& faces, const labelRange& range);

        //- Classify the face types and set the element lists.
        //  The indirect addressing can be used when classifying groups of
        //  face (eg, from a faceZone etc) with an optional flipMap.
        //  The optional exclude marker can be used to skip faces on particular
        //  boundary types or regions.
        void classify
        (
            const UList<face>& faces,
            const labelUList& addr,
            const boolList& flipMap = boolList(),
            const bitSet& exclude = bitSet()
        );


        //- Clear any demand-driven data
        void clearOut();

        //- Set addressable sizes to zero, free up addressing memory.
        void clear();

        //- Sum element counts across all processes.
        void reduce();

        //- Sort element lists numerically.
        void sort();


    // Output

        //- Write information about the object as a dictionary,
        //- optionally write all element addresses
        virtual void writeDict(Ostream& os, const bool full=false) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightFacesI.H"

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::ensightMesh

Description
    Encapsulation of volume meshes for writing in ensight format.

SourceFiles
    ensightMesh.C
    ensightMeshIO.C
    ensightMeshOptions.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "ensightCells.H"
#include "ensightFaces.H"
#include "ensightGeoFile.H"
#include "Map.H"
#include "scalarField.H"
#include "wordRes.H"
#include "globalIndex.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class ensightMesh;

/*---------------------------------------------------------------------------*\
                         Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
public:

    // Forward Declarations
    class options;

    //- The zone-id for internal mesh or unzoned cells.
    static const label internalZone;


private:

    // Private Data

        //- Writer options
        const std::unique_ptr<options> options_;

        //- Reference to the OpenFOAM mesh
        const fvMesh& mesh_;

        //- Volume elements per cellZone, lookup by zone index.
        //  The zone -1 is reserved for internal mesh (unzoned cells)
        Map<ensightCells> cellZoneParts_;

        //- Face elements per faceZone, lookup by zone index.
        Map<ensightFaces> faceZoneParts_;

        //- Face elements per selected patch, lookup by patch index
        Map<ensightFaces> boundaryParts_;

        //- Track if it needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Clear all storage
        void clear();

        //- Enforce consistent index/part numbering
        void renumber();

        //- Write the poly connectivity
        void writePolysConnectivity
        (
            ensightGeoFile& os,
            const labelUList& polys,
            const labelList& pointToGlobal
        ) const;

        //- Write cell connectivity for specified type
        void writeCellConnectivity
        (
            ensightGeoFile& os,
            const ensightCells::elemType etype,
            const ensightCells& part,
            const labelList& pointToGlobal
        ) const;

        //- Write the regular cell connectivity for all types
        void writeCellConnectivity
        (
            ensightGeoFile& os,
            const ensightCells& part,
            const labelList& pointToGlobal
        ) const;


        //- No copy construct
        ensightMesh(const ensightMesh&) = delete;

        //- No copy assignment
        void operator=(const ensightMesh&) = delete;


public:

    // Constructors

        //- Construct from components
        ensightMesh(const fvMesh& mesh, const options& opts);

        //- Construct from fvMesh with all default options, binary output
        explicit ensightMesh(const fvMesh& mesh);

        //- Construct from fvMesh with all default options and specified format
        ensightMesh(const fvMesh& mesh, const IOstream::streamFormat format);


    // Member Functions

    // Access

        //- Reference to the underlying fvMesh
        inline const fvMesh& mesh() const;

        //- Reference to the writer/mesh options
        inline const ensightMesh::options& option() const;

        //- Ascii/Binary file output
        inline IOstream::streamFormat format() const;

        //- Face elements per selected patch, lookup by patch index
        //  Process in sorted order.
        //  May require special treatment for zone -1 (internal).
        inline const Map<ensightCells>& cellZoneParts() const;

        //- Face elements per faceZone, lookup by zone index.
        //  Process in sorted order.
        inline const Map<ensightFaces>& faceZoneParts() const;

        //- Face elements per selected patch, lookup by patch index
        //  Process in sorted order.
        inline const Map<ensightFaces>& boundaryParts() const;


    // Other

        //- Does the content need an update?
        bool needsUpdate() const;

        //- Mark as needing an update.
        //  May also free up unneeded data.
        //  Return false if already marked as expired.
        bool expire();

        //- Update for new mesh
        void correct();


    // Output

        //- Write to file
        inline void write(autoPtr<ensightGeoFile>& os) const;

        //- Write to file
        void write(ensightGeoFile& os) const;
};


//- Configuration options for the ensightMesh
class ensightMesh::options
{
    // Private Data

        //- Ascii/Binary file output
        IOstream::streamFormat format_;

        //- Create in 'expired' mode
        bool lazy_;

        //- Use the internal mesh
        bool internal_;

        //- Use the boundary mesh
        bool boundary_;

        //- Handle cellZones (if internal_ is true)
        bool cellZones_;

        //- Selected patches only
        wordRes patchInclude_;

        //- Deselected patches
        wordRes patchExclude_;

        //- Selected cellZones
        wordRes cellZoneInclude_;

        //- Selected faceZones
        wordRes faceZoneInclude_;


public:

    // Constructors

        //- Default construct, use binary output
        options();

        //- Construct for specified format
        explicit options(IOstream::streamFormat format);


    // Member Functions

    // Access

        //- File output format (ascii | binary)
        IOstream::streamFormat format() const;

        //- Lazy creation? (ie, ensightMesh starts as needsUpdate)
        bool lazy() const;

        //- Using internal?
        bool useInternalMesh() const;

        //- Using boundary?
        bool useBoundaryMesh() const;

        //- Using faceZones?
        bool useFaceZones() const;

        //- Using cellZones?
        bool useCellZones() const;

        //- Selection of patches. Empty if unspecified.
        const wordRes& patchSelection() const;

        //- Selection of black listed patches. Empty if unspecified.
        const wordRes& patchExclude() const;

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& faceZoneSelection() const;

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& cellZoneSelection() const;


    // Edit

        //- Reset to defaults
        void reset();

        //- Lazy creation - ensightMesh starts as needsUpdate
        void lazy(bool beLazy);

        //- Alter the useBoundaryMesh state
        void useInternalMesh(bool on);

        //- Alter the useBoundaryMesh state
        void useBoundaryMesh(bool on);

        //- Alter the useCellZones state
        void useCellZones(bool on);

        //- Define patch selection matcher
        void patchSelection(const UList<wordRe>& patterns);

        //- Define patch selection matcher
        void patchSelection(List<wordRe>&& patterns);

        //- Define patch selection blacklist
        void patchExclude(const UList<wordRe>& patterns);

        //- Define patch selection blacklist
        void patchExclude(List<wordRe>&& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(const UList<wordRe>& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(List<wordRe>&& patterns);

        //- Define cellZone selection matcher
        void cellZoneSelection(const UList<wordRe>& patterns);

        //- Define cellZone selection matcher
        void cellZoneSelection(List<wordRe>&& patterns);


    // Edit

        //- Report values
        void print(Ostream& os) const;


    // Housekeeping

        //- Older name for useBoundaryMesh()
        //  \deprecated(2018-10)
        bool usePatches() const { return useBoundaryMesh(); }

        //- Older name for useBoundaryMesh()
        //  \deprecated(2018-10)
        void noPatches(bool off) { useBoundaryMesh(!off); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightMeshI.H"

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "IOmanip.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit each component
    template<class Type>
    static inline void writeData(Ostream& os, const Type& val)
    {
        for (direction cmpt=0; cmpt <  pTraits<Type>::nComponents; ++cmpt)
        {
            os  << ' ' << component(val, cmpt);
        }
        os  << nl;
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::surfaceWriters::abaqusWriter::writeFaceValue
(
    Ostream& os,
    const Type& value,
    const label elemId
) const
{
    os << elemId;
    writeData(os, value);

    return os;
}


template<class Type>
Foam::fileName Foam::surfaceWriters::abaqusWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Field:  rootdir/<TIME>/field/surfaceName.abq

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }
    outputFile /= fieldName / outputPath_.name();
    outputFile.ext("raw");  // stop-gap


    // Output scaling for the variable, but not for integer types.
    // could also solve with clever templating

    const scalar varScale =
    (
        std::is_integral<Type>::value
      ? scalar(1)
      : fieldScale_.getOrDefault<scalar>(fieldName, 1)
    );

    if (verbose_)
    {
        Info<< "Writing field " << fieldName;
        if (!equal(varScale, 1))
        {
            Info<< " (scaling " << varScale << ')';
        }
        Info<< " to " << outputFile << endl;
    }


    // geometry merge() implicit
    tmp<Field<Type>> tfield = mergeField(localValues);

    const meshedSurf& surf = surface();

    if (Pstream::master() || !parallel_)
    {
        const auto& values = tfield();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        // const scalar timeValue(0);


        // Additional bookkeeping for decomposing non tri/quad
        labelList decompOffsets;
        DynamicList<face> decompFaces;

        // Separate geometry

        if (noGeometry_)
        {
            fileFormats::ABAQUSCore::faceDecomposition
            (
                surf.points(),
                surf.faces(),
                decompOffsets,
                decompFaces
            );
        }
        else
        {
            OFstream osGeom(outputFile.lessExt().ext("abq"));
            writeGeometry(osGeom, surf, decompOffsets, decompFaces);
        }

        // Raw output values. To be improved

        OFstream os(outputFile);

        if (verbose_)
        {
            Info<< "Writing raw (abaqus) file to " << os.name() << endl;
            // Info<< "Writing raw abaqus file to " << os.name() << endl;
        }


        // Regular (undecomposed) faces
        const faceList&    faces = surf.faces();
        const labelList& elemIds = surf.faceIds();


        // Possible to use faceIds?
        // Not possible with on-the-fly face decomposition
        const bool useOrigFaceIds =
        (
            elemIds.size() == faces.size()
         && decompFaces.empty()
        );


        label elemId = 0;

        if (this->isPointData())
        {
            forAll(faces, facei)
            {
                if (useOrigFaceIds)
                {
                    // When available and not decomposed
                    elemId = elemIds[facei];
                }

                const label beginElemId = elemId;

                // Any face decomposition
                for
                (
                    label decompi = decompOffsets[facei];
                    decompi < decompOffsets[facei+1];
                    ++decompi
                )
                {
                    const face& f = decompFaces[decompi];

                    Type v = Zero;
                    for (const label verti : f)
                    {
                        v += values[verti];
                    }
                    v *= (varScale / f.size());

                    writeFaceValue(os, v, ++elemId);
                }


                // Face not decomposed
                if (beginElemId == elemId)
                {
                    const face& f = faces[facei];

                    Type v = Zero;
                    for (const label verti : f)
                    {
                        v += values[verti];
                    }
                    v *= (varScale / f.size());

                    writeFaceValue(os, v, ++elemId);
                }
            }
        }
        else
        {
            auto valIter = values.cbegin();

            forAll(faces, facei)
            {
                if (useOrigFaceIds)
                {
                    // When available and not decomposed
                    elemId = elemIds[facei];
                }

                const Type v(varScale * *valIter);
                ++valIter;

                label nValues =
                    max
                    (
                        label(1),
                        (decompOffsets[facei+1] - decompOffsets[facei])
                    );

                while (nValues--)
                {
                    writeFaceValue(os, v, ++elemId);
                }
            }
        }

        writeFooter(os, surf);
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

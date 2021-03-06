/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::lumpedPointState

Description
    The \a state of lumped points corresponds to positions and rotations.

    This class encapsulates the response from the external application and
    serves as the entry point for applying relaxation, sub-stepping etc.

    \heading Dictionary input format
    \table
        Property    | Description                           | Required | Default
        points      | List of points                        | yes |
        angles      | List of Euler rotation angles         | yes |
        order       | The Euler-angle rotation order        | no  | zxz
        degrees     | Rotation angles in degrees            | no  | false
    \endtable

    \heading Plain input format.
    Blank and comment lines starting with a '#' character are ignored.
    The angles are always in radians.
    \verbatim
    NumPoints
    x0 y0 z0 eulerz0 eulerx'0 eulerz''0
    x1 y1 z1 eulerz1 eulerx'1 eulerz''1
    ...
    \endverbatim

SeeAlso
    Foam::coordinateRotations::euler, Foam::quaternion

SourceFiles
    lumpedPointState.C
    lumpedPointStateI.H

\*---------------------------------------------------------------------------*/

#ifndef lumpedPointState_H
#define lumpedPointState_H

#include "dictionary.H"
#include "scalarList.H"
#include "pointField.H"
#include "scalarField.H"
#include "vectorField.H"
#include "tensorField.H"
#include "quaternion.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                      Class lumpedPointState Declaration
\*---------------------------------------------------------------------------*/

class lumpedPointState
{
public:

    //- Input format types
    enum class inputFormatType
    {
        PLAIN,      //!< "plain" is a simple ASCII format
        DICTIONARY  //!< "dictionary" is the OpenFOAM dictionary format
    };

    // Static data

    //- Names for the input format types
    static const Enum<inputFormatType> formatNames;


private:

    // Private Data

        //- Positions of lumped points
        pointField  points_;

        //- Orientation of lumped points (as Euler angles)
        vectorField angles_;

        //- The Euler-angle rotation order (default: zxz)
        quaternion::eulerOrder order_;

        //- Euler angles measured in degrees
        bool degrees_;

        //- Tensor rotation of lumped points
        mutable tensorField* rotationPtr_;


    // Private Member Functions

        void calcRotations() const;

        void readDict(const dictionary& dict);

public:

    // Constructors

        //- Construct null
        lumpedPointState();

        //- Copy constructor
        lumpedPointState(const lumpedPointState& rhs);

        //- Construct from points with zero-rotation
        lumpedPointState(const pointField& pts);

        //- Construct from points with zero-rotation
        lumpedPointState(tmp<pointField>& pts);

        //- Construct from dictionary
        lumpedPointState(const dictionary& dict);


    //- Destructor
    virtual ~lumpedPointState();


    // Member Functions

        //- Has positions and consistent number of rotations?
        inline bool valid() const;

        //- If no points were specified
        inline bool empty() const;

        //- The number of points
        inline label size() const;

        //- The points corresponding to mass centres
        inline const pointField& points() const;

        //- The orientation of the points (mass centres)
        inline const vectorField& angles() const;

        //- The local-to-global transformation for each point
        inline const tensorField& rotations() const;

        //- Scale points by given factor.
        //  Zero and negative values are ignored.
        void scalePoints(const scalar scaleFactor);

        //- Relax the state
        //  alpha = 1 : no relaxation
        //  alpha < 1 : relaxation
        //  alpha = 0 : do nothing
        void relax(const scalar alpha, const lumpedPointState& prev);

        //- Read input as dictionary content
        bool readData(Istream& is);

        //- Output as dictionary content
        bool writeData(Ostream& os) const;

        //- Output as dictionary content
        void writeDict(Ostream& os) const;

        //- Read input as plain content
        bool readPlain(Istream& is);

        //- Output as plain content
        void writePlain(Ostream& os) const;

        //- Read input from file (master process only) using specified format
        bool readData(const inputFormatType& fmt, const fileName& file);

        //- Output as VTK file for debugging/visualization
        //  The points are joined as lines, the rotation is visualized
        //  by planes, write as VTK PolyData format.
        void writeVTP(const fileName& outputFile, const vector& axis) const;


    // Member Operators

        //- Assignment operator
        void operator=(const lumpedPointState& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lumpedPointStateI.H"

#endif

// ************************************************************************* //

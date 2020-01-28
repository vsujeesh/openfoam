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

\*---------------------------------------------------------------------------*/

#include "ensightMesh.H"
#include "ensightOutput.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightMesh::writePolysConnectivity
(
    ensightGeoFile& os,
    const labelUList& addr,
    const labelList& pointToGlobal
) const
{
    // Number of faces per polyhedral (1/line in ASCII)
    {
        labelList send
        (
            ensightOutput::Detail::getPolysNFaces(mesh_, addr)
        );

        if (Pstream::master())
        {
            // Master
            os.writeLabels(send);

            // Slaves
            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

                labelList recv(fromSlave);
                os.writeLabels(recv);
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }


    // Number of points for each polyhedral face (1/line in ASCII)
    {
        labelList send
        (
            ensightOutput::Detail::getPolysNPointsPerFace(mesh_, addr)
        );

        if (Pstream::master())
        {
            // Master
            os.writeLabels(send);

            // Slaves
            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

                labelList recv(fromSlave);
                os.writeLabels(recv);
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }


    // List of points id for each face of the above list
    if (Pstream::master())
    {
        // Master
        ensightOutput::writePolysPoints
        (
            os,
            mesh_,
            addr,
            pointToGlobal
        );

        // Slaves
        for (int slave=1; slave<Pstream::nProcs(); ++slave)
        {
            IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

            cellList  cells(fromSlave);
            labelList addr(fromSlave);
            faceList  faces(fromSlave);
            labelList owner(fromSlave);

            ensightOutput::writePolysPoints
            (
                os,
                cells,
                addr,
                faces,
                owner
            );
        }
    }
    else
    {
        // Renumber faces to use global point numbers
        faceList faces(mesh_.faces());
        ListListOps::inplaceRenumber(pointToGlobal, faces);

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        toMaster
            << mesh_.cells()
            << addr
            << faces
            << mesh_.faceOwner();
    }
}


void Foam::ensightMesh::writeCellConnectivity
(
    ensightGeoFile& os,
    const ensightCells::elemType etype,
    const ensightCells& part,
    const labelList& pointToGlobal
) const
{
    const label nTotal = part.total(etype);

    if (nTotal)
    {
        const labelUList& addr = part.cellIds(etype);

        if (Pstream::master())
        {
            os.writeKeyword(ensightCells::key(etype));
            os.write(nTotal);
            os.newline();
        }

        if (etype == ensightCells::NFACED)
        {
            writePolysConnectivity(os, addr, pointToGlobal);
            return;
        }


        // Primitive shape - get subset and renumber
        cellShapeList shapes(mesh_.cellShapes(), addr);

        ListListOps::inplaceRenumber(pointToGlobal, shapes);

        if (Pstream::master())
        {
            ensightOutput::writeCellShapes(os, shapes);

            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                cellShapeList recv(fromSlave);

                ensightOutput::writeCellShapes(os, recv);
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << shapes;
        }
    }
}


void Foam::ensightMesh::writeCellConnectivity
(
    ensightGeoFile& os,
    const ensightCells& part,
    const labelList& pointToGlobal
) const
{
    for (label typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const auto etype = ensightCells::elemType(typei);

        writeCellConnectivity(os, etype, part, pointToGlobal);
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "facePointInterpolation.H"
#include "emptyFvPatch.H"
#include "pointConstraints.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::facePointInterpolation::interpolate
(
    const label patchi,
    const UList<Type>& ff
) const
{
    const polyPatch& pp = mesh().boundaryMesh()[patchi];

    const labelList& meshPoints = pp.meshPoints();
    const scalarListList& pointWeights = pointWeights_[patchi];

    tmp<Field<Type>> tpf(new Field<Type>(meshPoints.size(), Zero));
    Field<Type>& pf = tpf.ref();

    // Do points on 'normal' patches from the surrounding patch faces
    forAll(meshPoints, i)
    {
        const labelList& pFaces = pp.pointFaces()[i];
        const scalarList& pWeights = pointWeights[i];

        Type& val = pf[i];

        val = Zero;
        forAll(pFaces, j)
        {
            val += pWeights[j]*ff[pFaces[j]];
        }
    }

    // Reduce coupled points
    if (Pstream::parRun())
    {
        const Map<label>& coupledPoints = coupledPoints_[patchi];
        Map<Type> pfMap;
        forAllConstIter(Map<label>, coupledPoints, iter)
        {
            pfMap.insert
            (
                iter.key(),
                pf[iter()]
            );
        }

        syncTools::syncPointMap
        (
            mesh(),
            pfMap,
            plusEqOp<Type>()
        );

        forAllConstIter(Map<label>, coupledPoints, iter)
        {
            pf[iter()] = pfMap[iter.key()];
        }
    }
    return tpf;
}

// ************************************************************************* //

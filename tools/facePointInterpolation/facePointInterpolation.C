/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
#include "polyMesh.H"
#include "pointConstraints.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(facePointInterpolation, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::facePointInterpolation::makeWeights() const
{
    if (debug)
    {
        Pout<< "newVolPointInterpolation::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    coupledPoints_.setSize(mesh().boundaryMesh().size());
    pointWeights_.setSize(mesh().boundaryMesh().size());

    HashSet<label> coupledPoints;
    if (Pstream::parRun())
    {
        forAll(mesh().boundaryMesh(), patchi)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchi];

            if (pp.coupled())
            {
                coupledPoints.insert(pp.meshPoints());
            }
        }
    }

    forAll(mesh().boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        if (pp.coupled() || isA<emptyPolyPatch>(pp))
        {
            continue;
        }

        const labelList& meshPoints = pp.meshPoints();

        if (Pstream::parRun())
        {
            // Set the coupled points
            coupledPoints_.set(patchi, new Map<label>());
            Map<label>& cPoints = coupledPoints_[patchi];

            // Add coupled points and their relative index
            forAll(meshPoints, i)
            {
                const label pointi = meshPoints[i];
                if (coupledPoints.found(pointi))
                {
                    cPoints.insert(pointi, i);
                }
            }
        }

        pointWeights_.set(patchi, new scalarListList(meshPoints.size()));
        scalarListList& pWeights = pointWeights_[patchi];

        scalarField sumWeights(pWeights.size(), 0.0);


        const pointField& points = mesh().points();
        const vectorField::subField faceCentres = pp.faceCentres();
        forAll(meshPoints, i)
        {
            label pointi = meshPoints[i];

            const labelList& pFaces = pp.pointFaces()[i];

            scalarList& pw = pWeights[i];
            pw.setSize(pFaces.size());

            sumWeights[i] = 0.0;

            forAll(pFaces, j)
            {
                pw[j] = 1.0/mag(points[pointi] - faceCentres[pFaces[j]]);
                sumWeights[i] += pw[j];
            }
        }

        // Reduce coupled points
        if (Pstream::parRun())
        {
            const Map<label>& coupledPoints = coupledPoints_[patchi];
            Map<scalar> sumWeightMap;
            forAllConstIter(Map<label>, coupledPoints, iter)
            {
                sumWeightMap.insert
                (
                    iter.key(),
                    sumWeights[iter()]
                );
            }

            syncTools::syncPointMap
            (
                mesh(),
                sumWeightMap,
                plusEqOp<scalar>()
            );

            forAllConstIter(Map<label>, coupledPoints, iter)
            {
                sumWeights[iter()] = sumWeightMap[iter.key()];
            }
        }

        // Normalise weights
        forAll(pWeights, i)
        {
            scalarList& pw = pWeights[i];
            forAll(pw, j)
            {
                pw[j] /= sumWeights[i];
            }
        }
    }


    if (debug)
    {
        Pout<< "newVolPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::facePointInterpolation::facePointInterpolation(const polyMesh& pm)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, facePointInterpolation>(pm)
{
    makeWeights();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::facePointInterpolation::~facePointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::facePointInterpolation::updateMesh(const mapPolyMesh&)
{
    makeWeights();
}


bool Foam::facePointInterpolation::movePoints()
{
    makeWeights();

    return true;
}

// ************************************************************************* //

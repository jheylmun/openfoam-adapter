/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::facePointInterpolation

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting

SourceFiles
    facePointInterpolation.C
    newVolPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef facePointInterpolation_H
#define facePointInterpolation_H

#include "MeshObject.H"
#include "scalarList.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                    Class facePointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class facePointInterpolation
:
    public MeshObject<polyMesh, UpdateableMeshObject, facePointInterpolation>
{
    // Private Data

        //- Interpolation scheme weighting factor array.
        mutable PtrList<scalarListList> pointWeights_;

        mutable PtrList<Map<label>> coupledPoints_;


    // Private Member Functions

        //- Construct all point weighting factors
        void makeWeights() const;


public:

    // Declare name of the class and its debug switch
    ClassName("facePointInterpolation");


    // Constructors

        //- Constructor given fvMesh and pointMesh.
        explicit facePointInterpolation(const polyMesh&);

        //- Disallow default bitwise copy construction
        facePointInterpolation(const facePointInterpolation&) = delete;


    //- Destructor
    ~facePointInterpolation();


    // Member Functions

        // Edit

            //- Update mesh topology using the morph engine
            void updateMesh(const mapPolyMesh&);

            //- Correct weighting factors for moving mesh.
            bool movePoints();


    // Interpolation functions

        //- Interpolate volField using inverse distance weighting
        //  returning pointField
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const label patchi,
            const UList<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const facePointInterpolation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "facePointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

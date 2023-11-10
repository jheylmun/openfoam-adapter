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
    Foam::timeVaryingFixedValuePointPatchVectorField

Description
    Foam::timeVaryingFixedValuePointPatchVectorField

SourceFiles
    timeVaryingFixedValuePointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingFixedValuePointPatchVectorField_H
#define timeVaryingFixedValuePointPatchVectorField_H

#include "fixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class timeVaryingFixedValuePointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingFixedValuePointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private Data

        //- Reference displacement
        vectorField refDisp_;

        //- Patch velocity
        vectorField velocity_;

        //- Reference time
        scalar tRef_;

        //- delta T
        scalar deltaT_;


public:

    //- Runtime type information
    TypeName("timeVaryingFixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingFixedValuePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingFixedValuePointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        timeVaryingFixedValuePointPatchVectorField
        (
            const timeVaryingFixedValuePointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        timeVaryingFixedValuePointPatchVectorField
        (
            const timeVaryingFixedValuePointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new timeVaryingFixedValuePointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const pointPatchFieldMapper&);

            //- Reverse map the given pointPatchField onto this pointPatchField
            //  Used to reconstruct fields
            virtual void rmap(const pointPatchField<vector>&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Store the current patch values
            virtual void saveRefState(const scalar deltaT);

            //- Store ref displacement, calculate velocity using current
            //  reference state
            virtual void saveRefState
            (
                const scalar deltaT,
                const vectorField& refDisp
            );

            //- Store ref displacement, velocity
            virtual void saveRefState
            (
                const scalar deltaT,
                const vectorField& refDisp,
                const vectorField& velocity
            );
            //- Update the velocity
            virtual void updateVelocity();



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

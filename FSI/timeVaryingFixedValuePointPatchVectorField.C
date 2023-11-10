/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "timeVaryingFixedValuePointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingFixedValuePointPatchVectorField::
timeVaryingFixedValuePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    refDisp_(p.size(), Zero),
    velocity_(p.size(), Zero),
    tRef_(0.0),
    deltaT_(small)
{}


timeVaryingFixedValuePointPatchVectorField::
timeVaryingFixedValuePointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    refDisp_(p.size(), Zero),
    velocity_(p.size(), Zero),
    tRef_
    (
        dict.lookupOrDefault
        (
            "tRef",
            iF.mesh().time().value()
        )
    ),
    deltaT_(iF.mesh().time().deltaTValue())
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("refDisp"))
    {
        refDisp_ = vectorField("refDisp", dict , p.size());
    }
    else
    {
        refDisp_ = *this;
    }
    if (dict.found("velocity"))
    {
        velocity_ = vectorField("velocity", dict , p.size());
    }
}


timeVaryingFixedValuePointPatchVectorField::
timeVaryingFixedValuePointPatchVectorField
(
    const timeVaryingFixedValuePointPatchVectorField& tvfv,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(tvfv, p, iF, mapper),
    refDisp_(mapper(tvfv.refDisp_)),
    velocity_(mapper(tvfv.velocity_)),
    tRef_(tvfv.tRef_),
    deltaT_(tvfv.deltaT_)
{}


timeVaryingFixedValuePointPatchVectorField::
timeVaryingFixedValuePointPatchVectorField
(
    const timeVaryingFixedValuePointPatchVectorField& tvfv,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(tvfv, iF),
    refDisp_(tvfv.refDisp_),
    velocity_(tvfv.velocity_),
    tRef_(tvfv.tRef_),
    deltaT_(tvfv.deltaT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingFixedValuePointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);
    m(refDisp_, refDisp_);
    m(velocity_, velocity_);
}


void timeVaryingFixedValuePointPatchVectorField::rmap
(
    const pointPatchField<vector>& ptf,
    const labelList& addr
)
{
    const timeVaryingFixedValuePointPatchVectorField& tvfv =
        refCast<const timeVaryingFixedValuePointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::rmap(tvfv, addr);
    refDisp_.rmap(tvfv.refDisp_, addr);
    velocity_.rmap(tvfv.velocity_, addr);
}


void timeVaryingFixedValuePointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
    Field<vector>::operator=(refDisp_ + (t.value() - tRef_)*velocity_);

    fixedValuePointPatchField<vector>::updateCoeffs();
}

void timeVaryingFixedValuePointPatchVectorField::saveRefState
(
    const scalar deltaT
)
{
    const scalar t = this->internalField().mesh().time().value();
    if (t >= tRef_)
    {
        deltaT_ = deltaT;
        velocity_ = ((*this) - refDisp_)/deltaT;
        refDisp_ = *this;
        tRef_ = t;
    }
    updateCoeffs();
}


void timeVaryingFixedValuePointPatchVectorField::saveRefState
(
    const scalar deltaT,
    const vectorField& refDisp
)
{
    const scalar t = this->internalField().mesh().time().value();

    deltaT_ = deltaT;
    velocity_ = (refDisp - refDisp_)/deltaT;
    refDisp_ = refDisp;
    tRef_ = t;

    updateCoeffs();
}


void timeVaryingFixedValuePointPatchVectorField::saveRefState
(
    const scalar deltaT,
    const vectorField& refDisp,
    const vectorField& velocity
)
{
    deltaT_ = deltaT;
    tRef_ = this->internalField().mesh().time().value();
    refDisp_ = refDisp;
    velocity_ = velocity;

    updateCoeffs();
}


void timeVaryingFixedValuePointPatchVectorField::updateVelocity()
{
    velocity_ = ((*this) - refDisp_)/deltaT_;
    updateCoeffs();
}

void timeVaryingFixedValuePointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    writeEntry(os, "refDisp", refDisp_);
    writeEntry(os, "velocity", velocity_);
    writeEntry(os, "tRef", tRef_);
    writeEntry(os, "deltaT", deltaT_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    timeVaryingFixedValuePointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

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
    p0_(*this),
    pDelta_(p.size(), Zero),
    t0_(0.0),
    deltaT_(1.0)
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
    pDelta_(p.size(), Zero),
    t0_
    (
        dict.lookupOrDefault
        (
            "t0",
            iF.mesh().time().value()
        )
    ),
    deltaT_
    (
        dict.lookupOrDefault
        (
            "deltaT",
            iF.mesh().time().deltaTValue()
        )
    )
{
    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    if (dict.found("p0"))
    {
        p0_ = vectorField("p0", dict , p.size());
    }
    else
    {
        p0_ = *this;
    }
    if (dict.found("pDelta"))
    {
        pDelta_ = vectorField("pDelta", dict , p.size());
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
    p0_(mapper(tvfv.p0_)),
    pDelta_(mapper(tvfv.pDelta_)),
    t0_(tvfv.t0_),
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
    p0_(tvfv.p0_),
    pDelta_(tvfv.pDelta_),
    t0_(tvfv.t0_),
    deltaT_(tvfv.deltaT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingFixedValuePointPatchVectorField::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<vector>::autoMap(m);
    m(p0_, p0_);
    m(pDelta_, pDelta_);
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
    p0_.rmap(tvfv.p0_, addr);
    pDelta_.rmap(tvfv.pDelta_, addr);
}


void timeVaryingFixedValuePointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
    const scalar f = (t.value() - t0_)/deltaT_;
    if (f > 1.0)
    {
        WarningInFunction
            << "Current time change is greater than the supported time step" << nl
            << t0_ << " "<<t0_ + deltaT_<< nl
            << endl;
    }
    else if (f < 0)
    {
        WarningInFunction
            << "Current time change is less than the initial time"
            << endl;
    }
    Field<vector>::operator=(p0_ + f*pDelta_);

    fixedValuePointPatchField<vector>::updateCoeffs();
}

void timeVaryingFixedValuePointPatchVectorField::save(const scalar deltaT)
{
    const scalar t = this->internalField().mesh().time().value();
    if (t >= t0_)
    {
        deltaT_ = deltaT;
        pDelta_ = (*this) - p0_;
        p0_ = *this;
        t0_ = t;
    }
    updateCoeffs();
}

void timeVaryingFixedValuePointPatchVectorField::update()
{
    pDelta_ = (*this) - p0_;
    updateCoeffs();
}

void timeVaryingFixedValuePointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    writeEntry(os, "p0", p0_);
    writeEntry(os, "pDelta", pDelta_);
    writeEntry(os, "t0", t0_);
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

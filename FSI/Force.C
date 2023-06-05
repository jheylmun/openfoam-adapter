#include "Force.H"

#if defined (BLAST_FOAM)
    #include "tractionBase.H"
    #include "lookupSolidModel.H"
#elif defined (SOLIDS_4_FOAM)
    #include "solidTractionFvPatchVectorField.H"
    #include "lookupSolidModel.H"
#endif

using namespace Foam;

preciceAdapter::FSI::Force::Force(
    const Foam::fvMesh& mesh,
    const std::string solverType,
    const bool usePoint)
: ForceBase(mesh, solverType, "Force", dimensionSet(1, 1, -2, 0, 0, 0, 0), usePoint)
{}

std::string preciceAdapter::FSI::Force::getDataName() const
{
    return "Force";
}

void preciceAdapter::FSI::Force::read(double* buffer, const unsigned int dim)
{
#if defined(BLAST_FOAM)
    // Lookup solid model
    const solidModel& solid = lookupSolidModel(mesh_);

    // Get the solution displacement field
    volVectorField& solD =
        const_cast<volVectorField&>(solid.solutionD());

    // Set boundary traction
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if (this->locationType_ == LocationType::faceCenters)
        {
            // Make a force field
            vectorField force(mesh_.boundaryMesh()[patchID].size(), vector::zero);

            // Copy the forces from the buffer into the force field
            forAll(force, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                {
                    force[i][d] = buffer[i * dim + d];
                }
            }

            // Calculate the traction field
            vectorField traction(force.size(), vector::zero);

            // Check if it is a linear or nonlinear geometry case
            if (solid.nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
            {
                traction = force/mesh_.boundary()[patchID].magSf();
            }
            else if (solid.nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
            {
                // Patch area vectors
                const vectorField& patchSf = mesh_.Sf().boundaryField()[patchID];

                // Lookup the inverse of the deformation gradient
                const tensorField& Finv =
                    mesh_.lookupObject<volTensorField>
                    (
                        "Finv"
                    ).boundaryField()[patchID];

                // Lookup the Jacobian
                const scalarField& J =
                    mesh_.lookupObject<volScalarField>
                    (
                        "J"
                    ).boundaryField()[patchID];

                // Calculate area vectors in the deformed configuration
                const scalarField patchDeformMagSf(mag(J*Finv.T() & patchSf));

                traction = force/patchDeformMagSf;
            }
            else if (solid.nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
            {
                // Patch area vectors
                const vectorField& patchSf = mesh_.Sf().boundaryField()[patchID];

                // Lookup the inverse of the relative deformation gradient
                const tensorField& relFinv =
                    mesh_.lookupObject<volTensorField>
                    (
                        "relFinv"
                    ).boundaryField()[patchID];

                // Lookup the relative Jacobian
                const scalarField& relJ =
                    mesh_.lookupObject<volScalarField>
                    (
                        "relJ"
                    ).boundaryField()[patchID];

                // Calculate area vectors in the deformed configuration
                const scalarField patchDeformMagSf
                (
                    mag(relJ*relFinv.T() & patchSf)
                );

                traction = force/patchDeformMagSf;
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown solidModel nonLinGeom type = "
                    << solid.nonLinGeom() << abort(FatalError);
            }

            // Cast boundary condition to solidTraction
            tractionBase& tractionPatch =
                refCast<tractionBase>
                (
                    solD.boundaryFieldRef()[patchID]
                );

            // Apply traction field
            tractionPatch.traction() = traction;
        }
        else if (this->locationType_ == LocationType::faceNodes)
        {
            NotImplemented;
        }
    }
#elif defined (SOLIDS_4_FOAM)
    // Lookup displacement field name
    const solidModel& solid = lookupSolidModel(mesh_);

    // Lookup the cell-centre displacement field
    volVectorField& solD =
        const_cast<solidModel&>(solid).solutionD();

    // Set boundary traction
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if (this->locationType_ == LocationType::faceCenters)
        {
            // Make a force field
            vectorField force(mesh_.boundaryMesh()[patchID].size(), vector::zero);

            // Copy the forces from the buffer into the force field
            forAll(force, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                {
                    force[i][d] = buffer[i * dim + d];
                }
            }

            // Calculate the traction field
            vectorField traction(force.size(), vector::zero);

            // Check if it is a linear or nonlinear geometry case
            if (solid.nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
            {
                traction = force/mesh_.boundary()[patchID].magSf();
            }
            else if (solid.nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
            {
                // Patch area vectors
                const vectorField& patchSf = mesh_.Sf().boundaryField()[patchID];

                // Lookup the inverse of the deformation gradient
                const tensorField& Finv =
                    mesh_.lookupObject<volTensorField>
                    (
                        "Finv"
                    ).boundaryField()[patchID];

                // Lookup the Jacobian
                const scalarField& J =
                    mesh_.lookupObject<volScalarField>
                    (
                        "J"
                    ).boundaryField()[patchID];

                // Calculate area vectors in the deformed configuration
                const scalarField patchDeformMagSf(mag(J*Finv.T() & patchSf));

                traction = force/patchDeformMagSf;
            }
            else if (solid.nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
            {
                // Patch area vectors
                const vectorField& patchSf = mesh_.Sf().boundaryField()[patchID];

                // Lookup the inverse of the relative deformation gradient
                const tensorField& relFinv =
                    mesh_.lookupObject<volTensorField>
                    (
                        "relFinv"
                    ).boundaryField()[patchID];

                // Lookup the relative Jacobian
                const scalarField& relJ =
                    mesh_.lookupObject<volScalarField>
                    (
                        "relJ"
                    ).boundaryField()[patchID];

                // Calculate area vectors in the deformed configuration
                const scalarField patchDeformMagSf
                (
                    mag(relJ*relFinv.T() & patchSf)
                );

                traction = force/patchDeformMagSf;
            }
            else
            {
                FatalErrorInFunction
                    << "Unknown solidModel nonLinGeom type = "
                    << solid.nonLinGeom() << abort(FatalError);
            }

            // Cast boundary condition to solidTraction
            solidTractionFvPatchVectorField& tractionPatch =
                refCast<solidTractionFvPatchVectorField>
                (
                    solD.boundaryFieldRef()[patchID]
                );

            // Apply traction field
            tractionPatch.traction() = traction;
        }
        else if (this->locationType_ == LocationType::faceNodes)
        {
            NotImplemented;
        }
    }
#else
    NotImplemented;
#endif
}

Foam::tmp<Foam::vectorField> preciceAdapter::FSI::Force::getFaceVectors(const unsigned int patchID) const
{
    // Normal vectors multiplied by face area
    return mesh_.boundary()[patchID].Sf();
}

preciceAdapter::FSI::Force::~Force()
{}

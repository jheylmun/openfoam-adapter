#include "DisplacementDelta.H"
#include "timeVaryingFixedValuePointPatchVectorField.H"

using namespace Foam;

preciceAdapter::FSI::DisplacementDelta::DisplacementDelta(
    const Foam::fvMesh& mesh,
    const std::string namePointDisplacementDelta,
    const std::string nameCellDisplacementDelta)
: pointDisplacementDelta_(
    const_cast<pointVectorField*>(
        &mesh.lookupObject<pointVectorField>(namePointDisplacementDelta))),
  cellDisplacementDelta_(
      const_cast<volVectorField*>(
          &mesh.lookupObject<volVectorField>(nameCellDisplacementDelta))),
  mesh_(mesh)
{
    dataType_ = vector;
}

// We cannot do this step in the constructor by design of the adapter since the information of the CouplingDataUser is
// defined later. Hence, we call this method after the CouplingDaaUser has been configured
void preciceAdapter::FSI::DisplacementDelta::initialize()
{
    // Initialize appropriate objects for each interface patch, namely the volField and the interpolation object
    // this is only necessary for face based FSI
    if (this->locationType_ == LocationType::faceCenters)
    {
        for (unsigned int j = 0; j < patchIDs_.size(); ++j)
        {
            const unsigned int patchID = patchIDs_.at(j);
            interpolationObjects_.emplace_back(new primitivePatchInterpolation(mesh_.boundaryMesh()[patchID]));
        }
    }
}

void preciceAdapter::FSI::DisplacementDelta::setDeltaT
(
    const double deltaT,
    const bool complete
)
{
    if (complete)
    {
        pointVectorField::Boundary& bpointDisplacementDelta =
            pointDisplacementDelta_->boundaryFieldRef();
        for (unsigned int j = 0; j < patchIDs_.size(); j++)
        {
            // Get the ID of the current patch
            const unsigned int patchID = patchIDs_.at(j);

            if
            (
                isA<timeVaryingFixedValuePointPatchVectorField>
                (
                    bpointDisplacementDelta[patchID]
                )
            )
            {
                refCast<timeVaryingFixedValuePointPatchVectorField>
                (
                    bpointDisplacementDelta[patchID]
                ).saveRefState(deltaT);
            }
        }
    }
    else
    {
        pointDisplacementDelta_->correctBoundaryConditions();
    }
}

void preciceAdapter::FSI::DisplacementDelta::update()
{
    pointVectorField::Boundary& bpointDisplacementDelta =
        pointDisplacementDelta_->boundaryFieldRef();
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if
        (
            isA<timeVaryingFixedValuePointPatchVectorField>
            (
                bpointDisplacementDelta[patchID]
            )
        )
        {
            refCast<timeVaryingFixedValuePointPatchVectorField>
            (
                bpointDisplacementDelta[patchID]
            ).updateVelocity();
        }
    }
}

void preciceAdapter::FSI::DisplacementDelta::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    if (this->locationType_ == LocationType::faceCenters)
    {
        // For every boundary patch of the interface
        for (const label patchID : patchIDs_)
        {
            const vectorField& pD =
                cellDisplacementDelta_->boundaryField()[patchID];

            // Write the forces to the preCICE buffer
            // For every face of the patch
            forAll(pD, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                {
                    buffer[i * dim + d] = pD[i][d];
                }
            }
        }
    }
    else if (this->locationType_ == LocationType::faceNodes)
    {
        // For every boundary patch of the interface
        for (const label patchID : patchIDs_)
        {
            // Write the forces to the preCICE buffer
            // For every point of the patch
            const labelList& meshPoints =
                mesh_.boundaryMesh()[patchID].meshPoints();
            forAll(meshPoints, i)
            {
                const label pointID = meshPoints[i];
                for (unsigned int d = 0; d < dim; ++d)
                {
                    buffer[i * dim + d] =
                        (*pointDisplacementDelta_)[pointID][d];
                }
            }
        }
    }
}

// return the displacement to use later in the velocity?
void preciceAdapter::FSI::DisplacementDelta::read(double* buffer, const unsigned int dim)
{
    DynamicList<Foam::vector> mappedDisplacementDelta;
    DynamicList<Foam::scalar> area;
    if (this->locationType_ == LocationType::faceCenters)
    {
        volVectorField::Boundary& bcellDisplacementDelta =
            cellDisplacementDelta_->boundaryFieldRef();
        pointVectorField::Boundary& bpointDisplacementDelta =
            pointDisplacementDelta_->boundaryFieldRef();

        for (unsigned int j = 0; j < patchIDs_.size(); j++)
        {
            // Get the ID of the current patch
            const unsigned int patchID = patchIDs_.at(j);

            // the boundaryCellDisplacement is a vector and ordered according to the iterator j
            // and not according to the patchID
            // First, copy the buffer data into the center based vectorFields on each interface patch
            // For DisplacementDelta, set absolute values here and sum the interpolated values up to the point field
            // since the temporary field in this class is not reloaded in the implicit coupling
            vectorField& pcellDisplacementDelta =
                bcellDisplacementDelta[patchID];

            forAll(bcellDisplacementDelta[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                {
                    pcellDisplacementDelta[i][d] = buffer[i * dim + d];
                }
            }
            if (log) mappedDisplacementDelta.append(pcellDisplacementDelta);
            if (log) area.append(mesh_.boundary()[patchID].magSf());

            // Get a reference to the displacement on the point patch in order to overwrite it
            vectorField& ppointDisplacementDelta
            (
                refCast<vectorField>
                (
                    bpointDisplacementDelta[patchID]
                )
            );

            // Overwrite the node based patch using the interpolation objects and the cell based vector field
            // Afterwards, continue as usual
            ppointDisplacementDelta =
                interpolationObjects_[j]->faceToPointInterpolate
                (
                    pcellDisplacementDelta
                );
        }
    }
    if (this->locationType_ == LocationType::faceNodes)
    {
        pointVectorField::Boundary& bpointDisplacementDelta =
            pointDisplacementDelta_->boundaryFieldRef();

        for (unsigned int j = 0; j < patchIDs_.size(); j++)
        {
            // Get the ID of the current patch
            const unsigned int patchID = patchIDs_.at(j);

            // Get the displacement on the patch
            vectorField& ppointDisplacementDelta =
                refCast<vectorField>
                (
                    bpointDisplacementDelta[patchID]
                );

            // Overwrite the nodes on the interface directly
            forAll(ppointDisplacementDelta, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                {
                    ppointDisplacementDelta[i][d] = buffer[i * dim + d];
                }
            }
            if (log) mappedDisplacementDelta.append(ppointDisplacementDelta);
        }
    }

    if (log)
    {
        if (this->locationType_ == LocationType::faceNodes)
        {
            printListStatistics
            (
                "pointDisplacementDelta",
                mappedDisplacementDelta
            );
        }
        else
        {
            printListStatistics
            (
                "cellDisplacementDelta",
                mappedDisplacementDelta,
                area
            );
        }
    }
}

bool preciceAdapter::FSI::DisplacementDelta::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::faceNodes);
}

std::string preciceAdapter::FSI::DisplacementDelta::getDataName() const
{
    return "DisplacementDelta";
}

#include "Displacement.H"
#include "timeVaryingFixedValuePointPatchVectorField.H"

using namespace Foam;

preciceAdapter::FSI::Displacement::Displacement(
    const Foam::fvMesh& mesh,
    const std::string namePointDisplacement,
    const std::string nameCellDisplacement)
: pointDisplacement_(
    const_cast<pointVectorField*>(
        &mesh.lookupObject<pointVectorField>(namePointDisplacement))),
  cellDisplacement_(
      const_cast<volVectorField*>(
          &mesh.lookupObject<volVectorField>(nameCellDisplacement))),
  mesh_(mesh)
{
    dataType_ = vector;
    if (logToFile)
    {
        if (this->locationType_ == LocationType::faceCenters)
        {
            osPtr_.set(new OFstream("log." + pointDisplacement_->name()));
        }
        else
        {
            osPtr_.set(new OFstream("log." + cellDisplacement_->name()));
        }
        LogListHeader(osPtr_(), false);
    }
}

// We cannot do this step in the constructor by design of the adapter since the information of the CouplingDataUser is
// defined later. Hence, we call this method after the CouplingDaaUser has been configured
void preciceAdapter::FSI::Displacement::initialize()
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


void preciceAdapter::FSI::Displacement::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    if (this->locationType_ == LocationType::faceCenters)
    {
        // For every boundary patch of the interface
        for (const label patchID : patchIDs_)
        {
            const vectorField& pD =
                cellDisplacement_->boundaryField()[patchID];

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
                    buffer[i * dim + d] = (*pointDisplacement_)[pointID][d];
                }
            }
        }
    }
}


// return the displacement to use later in the velocity?
void preciceAdapter::FSI::Displacement::read(double* buffer, const unsigned int dim)
{
    pointVectorField::Boundary& bpointDisplacement =
        pointDisplacement_->boundaryFieldRef();
    DynamicList<Foam::vector> mappedDisplacement;
    DynamicList<Foam::scalar> area;
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        // Get a reference to the displacement on the point patch in order to overwrite it
        vectorField& ppointDisplacement =
            refCast<vectorField>(bpointDisplacement[patchID]);
        if (this->locationType_ == LocationType::faceCenters)
        {
            vectorField& pcellDisplacement =
                cellDisplacement_->boundaryFieldRef()[patchID];

            // the boundaryCellDisplacement is a vector and ordered according to the iterator j
            // and not according to the patchID
            // First, copy the buffer data into the center based vectorFields on each interface patch
            forAll(pcellDisplacement, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    pcellDisplacement[i][d] = buffer[i * dim + d];
            }
            if (log()) mappedDisplacement.append(pcellDisplacement);
            if (log()) area.append(mesh_.boundary()[patchID].magSf());

            // Overwrite the node based patch using the interpolation objects and the cell based vector field
            // Afterwards, continue as usual
            ppointDisplacement =
                interpolationObjects_[j]->faceToPointInterpolate
                (
                    cellDisplacement_->boundaryField()[patchID]
                );
        }
        else if (this->locationType_ == LocationType::faceNodes)
        {
            // Overwrite the nodes on the interface directly
            forAll(ppointDisplacement, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    ppointDisplacement[i][d] = buffer[i * dim + d];
            }
            if (log()) mappedDisplacement.append(ppointDisplacement);
        }
    }
    if (logToTerminal)
    {
        if (this->locationType_ == LocationType::faceNodes)
        {
            printListStatistics("pointDisplacement", mappedDisplacement);
        }
        else
        {
            printListStatistics("cellDisplacement", mappedDisplacement, area);
        }
    }
    if (logToFile)
    {
        if (this->locationType_ == LocationType::faceNodes)
        {
            LogListStatistics(osPtr_(), mesh_.time(), mappedDisplacement);
        }
        else
        {
            LogListStatistics(osPtr_(), mesh_.time(), mappedDisplacement, area);
        }
    }
}

void preciceAdapter::FSI::Displacement::setDeltaT
(
    const double deltaT,
    const bool complete
)
{
    if (complete)
    {
        pointVectorField::Boundary& bpointDisplacement =
            pointDisplacement_->boundaryFieldRef();
        for (unsigned int j = 0; j < patchIDs_.size(); j++)
        {
            // Get the ID of the current patch
            const unsigned int patchID = patchIDs_.at(j);

            if
            (
                isA<timeVaryingFixedValuePointPatchVectorField>
                (
                    bpointDisplacement[patchID]
                )
            )
            {
                refCast<timeVaryingFixedValuePointPatchVectorField>
                (
                    bpointDisplacement[patchID]
                ).saveRefState(deltaT);
            }
        }
    }

    pointDisplacement_->correctBoundaryConditions();
}

void preciceAdapter::FSI::Displacement::update()
{
    pointVectorField::Boundary& bpointDisplacement =
        pointDisplacement_->boundaryFieldRef();
    for (unsigned int j = 0; j < patchIDs_.size(); j++)
    {
        // Get the ID of the current patch
        const unsigned int patchID = patchIDs_.at(j);

        if
        (
            isA<timeVaryingFixedValuePointPatchVectorField>
            (
                bpointDisplacement[patchID]
            )
        )
        {
            refCast<timeVaryingFixedValuePointPatchVectorField>
            (
                bpointDisplacement[patchID]
            ).updateVelocity();
        }
    }
}

bool preciceAdapter::FSI::Displacement::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::faceNodes);
}

std::string preciceAdapter::FSI::Displacement::getDataName() const
{
    return "Displacement";
}

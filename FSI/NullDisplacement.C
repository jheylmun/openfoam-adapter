#include "NullDisplacement.H"

using namespace Foam;

preciceAdapter::FSI::NullDisplacement::NullDisplacement(
    const Foam::fvMesh& mesh,
    const std::string namePointDisplacement,
    const std::string nameCellDisplacement)
{}

// We cannot do this step in the constructor by design of the adapter since the information of the CouplingDataUser is
// defined later. Hence, we call this method after the CouplingDaaUser has been configured
void preciceAdapter::FSI::NullDisplacement::initialize()
{}


void preciceAdapter::FSI::NullDisplacement::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    NotImplemented;
}


// return the displacement to use later in the velocity?
void preciceAdapter::FSI::NullDisplacement::read(double* buffer, const unsigned int dim)
{}

bool preciceAdapter::FSI::NullDisplacement::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters || this->locationType_ == LocationType::faceNodes);
}

std::string preciceAdapter::FSI::NullDisplacement::getDataName() const
{
    return "Displacement";
}

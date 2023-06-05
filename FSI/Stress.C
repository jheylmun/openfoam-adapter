#include "Stress.H"

using namespace Foam;

preciceAdapter::FSI::Stress::Stress(
    const Foam::fvMesh& mesh,
    const std::string solverType,
    const bool usePoint)
: ForceBase(mesh, solverType, "Stress", dimensionSet(1, -1, -2, 0, 0, 0, 0), usePoint)
{}

void preciceAdapter::FSI::Stress::read(double* buffer, const unsigned int dim)
{
    this->readFromBuffer(buffer);
}

std::string preciceAdapter::FSI::Stress::getDataName() const
{
    return "Stress";
}

Foam::tmp<Foam::vectorField> preciceAdapter::FSI::Stress::getFaceVectors(const unsigned int patchID) const
{
    // face normal vectors
    return mesh_.boundary()[patchID].nf();
}

preciceAdapter::FSI::Stress::~Stress()
{}

#ifndef FORCEBASE_H
#define FORCEBASE_H


#include "CouplingDataUser.H"

#include "fvCFD.H"

#include "pointFields.H"
#include "vectorField.H"
#include "compressibleMomentumTransportModel.H"
#include "kinematicMomentumTransportModel.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads force
class ForceBase : public CouplingDataUser
{
protected:
    //- Stress tensor (see the OpenFOAM "Forces" function object)
    Foam::tmp<Foam::volSymmTensorField> devRhoReff() const;

    Foam::tmp<Foam::volScalarField> rho() const;

    Foam::tmp<Foam::volScalarField> mu() const;

    //- OpenFOAM fvMesh object (we need to access the objects' registry multiple times)
    const Foam::fvMesh& mesh_;

    const std::string solverType_;

    Foam::volVectorField* Force_;
    Foam::pointVectorField* pointForce_;

    //- Logs
    mutable autoPtr<OFstream> osAreaPtr_;
    mutable autoPtr<OFstream> osForcePtr_;
    mutable autoPtr<OFstream> osPointForcePtr_;

public:
    //- Constructor
    ForceBase(
        const Foam::fvMesh& mesh,
        const std::string solverType,
        const std::string forceName,
        const dimensionSet forceDims,
        const bool usePoint);

    virtual ~ForceBase();

    //- Write the forces values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    //- Read the forces values from the buffer
    void readFromBuffer(double* buffer) const;

    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    virtual Foam::tmp<Foam::vectorField> getFaceVectors(const unsigned int patchID) const = 0;
};
}
}

#endif // FORCEBASE_H

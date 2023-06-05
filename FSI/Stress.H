#ifndef FSI_STRESS_H
#define FSI_STRESS_H

#include "ForceBase.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads stress [N/m^2]:
// This is essentially a force (in spatial coordinates) scaled by the deformed
// cell face. Thus, a consistent quantity. Calculation concept has been copied
// from the force module, but the scaled version here is commonly used in FEM
// applications.
class Stress : public ForceBase
{
public:
    //- Constructor
    Stress(
        const Foam::fvMesh& mesh,
        const std::string solverType,
        const bool usePoint);

    //- Read the stress values from the buffer
    void read(double* buffer, const unsigned int dim) override;

    //- Get the name of the current data field
    std::string getDataName() const override;

    //- Returns the face normal vectors (no multiplication by area)
    Foam::tmp<Foam::vectorField> getFaceVectors(const unsigned int patchID) const override;

    //- Destructor
    ~Stress();
};
}
}

#endif

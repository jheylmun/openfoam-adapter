#ifndef FSI_FORCE_H
#define FSI_FORCE_H

#include "ForceBase.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads force
class Force : public ForceBase
{
public:
    //- Constructor
    Force(
        const Foam::fvMesh& mesh,
        const std::string solverType,
        const bool usePoint);

    //- Read the forces values from the buffer
    void read(double* buffer, const unsigned int dim) override;

    //- Get the name of the current data field
    std::string getDataName() const override;

    //- Returns the normal vectors multiplied by the face area
    Foam::tmp<Foam::vectorField> getFaceVectors(const unsigned int patchID) const override;

    //- Destructor
    ~Force();
};

}
}

#endif

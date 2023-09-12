#ifndef FSI_NULLDISPLACEMENT_H
#define FSI_NULLDISPLACEMENT_H

#include "CouplingDataUser.H"

#include "fvCFD.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads displacement
class NullDisplacement : public CouplingDataUser
{

private:

public:
    //- Constructor
    NullDisplacement(
        const Foam::fvMesh& mesh,
        const std::string namePointDisplacement,
        const std::string nameCellDisplacement);

    //- Write the displacement values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    //- Read the displacement values from the buffer
    void read(double* buffer, const unsigned int dim) override;

    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    //- Get the name of the current data field
    std::string getDataName() const override;

    //- We need to initialize the cell-based vector and the interpolation object
    // in case we want to use the faceCenter location for the coupling
    void initialize() override;
};

}
}

#endif

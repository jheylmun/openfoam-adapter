#ifndef FSI_DISPLACEMENTDELTA_H
#define FSI_DISPLACEMENTDELTA_H

#include "CouplingDataUser.H"

#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"
#include "primitivePatchInterpolation.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads displacement
class DisplacementDelta : public CouplingDataUser
{

private:
    // DisplacementDelta pointVectorField
    Foam::pointVectorField* pointDisplacementDelta_;

    // cellDisplacementDelta field
    Foam::volVectorField* cellDisplacementDelta_;

    const Foam::fvMesh& mesh_;
    // NOTE: this allocation could be avoided in case we directly write to the
    // Needs to be a pointer since the class disables assignment and copy constructors
    std::vector<const Foam::primitivePatchInterpolation*> interpolationObjects_;

    // log file
    mutable autoPtr<OFstream> osPtr_;

public:
    //- Constructor
    DisplacementDelta(
        const Foam::fvMesh& mesh,
        const std::string namePointDisplacementDelta,
        const std::string nameCellDisplacementDelta);

    //- Set the time step
    void setDeltaT(const double deltaT, const bool complete) override;

    //- Update
    virtual void update() override;

    //- Write the displacementDelta values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    //- Read the displacementDelta values from the buffer
    void read(double* buffer, const unsigned int dim) override;

    bool isLocationTypeSupported(const bool meshConnectivity) const override;

    //- Get the name of the current data field
    virtual std::string getDataName() const;

    //- We need to initialize the cell-based vector and the interpolation object
    // in case we want to use the faceCenter location for the coupling
    void initialize() override;
};

}
}

#endif

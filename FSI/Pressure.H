#ifndef FSI_PRESSURE_H
#define FSI_PRESSURE_H


#include "CouplingDataUser.H"

#include "fvCFD.H"

#include "pointFields.H"
#include "scalarField.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads force
class Pressure : public CouplingDataUser
{
protected:

    Foam::tmp<Foam::volScalarField> rho() const;


    //- OpenFOAM fvMesh object (we need to access the objects' registry multiple times)
    const Foam::fvMesh& mesh_;

    const std::string solverType_;

    //- Logs
    mutable autoPtr<OFstream> osAreaPtr_;
    mutable autoPtr<OFstream> osPressurePtr_;
    mutable autoPtr<OFstream> osPointPressurePtr_;

    const bool negPressure_;
    const bool pointPressure_;

public:
    //- Constructor
    Pressure(
        const Foam::fvMesh& mesh,
        const std::string solverType,
        const bool negPressure,
        const bool usePoint);

    virtual ~Pressure();

    //- Get the name of the current data field
    std::string getDataName() const override;

    //- Read the forces values from the buffer
    void read(double* buffer, const unsigned int dim) override;

    //- Write the forces values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim) override;

    bool isLocationTypeSupported(const bool meshConnectivity) const override;
};
}
}

#endif // FSI_PRESSURE_H

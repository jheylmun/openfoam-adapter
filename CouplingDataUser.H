#ifndef COUPLINGDATAUSER_H
#define COUPLINGDATAUSER_H

#include "Utilities.H"
#include <vector>
#include <string>

namespace preciceAdapter
{
// A small enum to deal with the different locations of the
// coupling interface
enum class LocationType
{
    none,
    faceCenters,
    faceNodes
};

class CouplingDataUser
{

protected:
    enum DataType
    {
        scalar,
        vector
    };

    //- Type of the coupling data (scalar or vector)
    DataType dataType_ = scalar;

    //- OpenFOAM patches that form the interface
    std::vector<int> patchIDs_;

    //- preCICE data ID
    int dataID_;

    //- location type of the interface
    LocationType locationType_ = LocationType::none;

public:
    //- Constructor
    CouplingDataUser();

    //- Returns true if the data are scalar
    bool hasScalarData();

    //- Returns true if the data are vector
    bool hasVectorData();

    //- Set the preCICE data ID
    void setDataID(int dataID);

    //- Get the preCICE data ID
    int dataID();

    //- Set the patch IDs that form the interface
    void setPatchIDs(std::vector<int> patchIDs);

    //- Set the locations type of the interface
    void setLocationsType(LocationType locationsType);

    // Check if the dataset supports this interface nodes location
    void checkDataLocation(const bool meshConnectivity) const;

    //- option to initialize data in derived data classes
    virtual void initialize();

    //- Write the coupling data to the buffer
    virtual void write(double* dataBuffer, bool meshConnectivity, const unsigned int dim) = 0;

    //- Read the coupling data from the buffer
    virtual void read(double* dataBuffer, const unsigned int dim) = 0;

    //- Given the meshConnectivity, return if the underlying loactionType of the
    //- interface nodes is supported by the data set
    virtual bool isLocationTypeSupported(const bool meshConnectivity) const = 0;

    //- Get the name of the current data field
    virtual std::string getDataName() const = 0;

    //- Destructor
    virtual ~CouplingDataUser() {}
};

}

#endif

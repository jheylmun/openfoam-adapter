#ifndef PRECICEADAPTER_H
#define PRECICEADAPTER_H

#include "Interface.H"

// Conjugate Heat Transfer module
#include "CHT/CHT.H"

// Fluid-Structure Interaction module
#include "FSI/FSI.H"

// Fluid-Fluid module
#include "FF/FF.H"

// NOTE: If you want to couple a new variable, include your module's header here.
// You also need to include it in the Make/files file.
// In case you use additional OpenFOAM symbols, you may also need to specify
// the respective libraries in the Make/options.

// OpenFOAM header files
#include "Time.H"
#include "fvMesh.H"

// preCICE Solver Interface
#include "precice/SolverInterface.hpp"

namespace preciceAdapter
{

class Adapter
{

private:
    //- Structure of the configuration of each coupling interface.
    //  Every interface needs to know the coupling mesh, the OpenFOAM
    //  patches that consist the coupling surface and the kinds
    //  of data that are exchanged.
    struct InterfaceConfig
    {
        std::string meshName;
        std::string locationsType;
        bool meshConnectivity;
        std::vector<std::string> patchNames;
        std::vector<std::string> writeData;
        std::vector<std::string> readData;
    };

    //- Configuration interfaces
    std::vector<struct InterfaceConfig> interfacesConfig_;

    //- OpenFOAM runTime object
    const Foam::Time& runTime_;

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // Configuration parameters used in the Adapter

    //- Remember if there were errors in the read() method
    bool errorsInConfigure = false;

    //- preCICE configuration file name
    Foam::word preciceConfigFilename_;

    //- preCICE participant name
    Foam::word participantName_;

    //- Switch for checkpointing
    bool checkpointing_ = false;

    //- Type of timestep (fixed, adjustable)
    bool adjustableTimestep_;

    //- Should the (fixed) timestep be stored and used?
    bool useStoredTimestep_ = false;

    //- Switch to enable the ConjugateHeatTransfer module
    bool CHTenabled_ = false;

    //- Switch to enable the FluidStructureInteraction module
    bool FSIenabled_ = false;

    //- Switch to enable the FluidFluid module
    bool FFenabled_ = false;

    // NOTE: Add a switch for your new module here

    //- Interfaces
    std::vector<Interface*> interfaces_;

    //- preCICE solver interface
    precice::SolverInterface* precice_ = NULL;

    //- preCICE solver interface initialized
    bool preciceInitialized_ = false;

    //- Conjugate Heat Transfer module object
    CHT::ConjugateHeatTransfer* CHT_ = NULL;

    //- Fluid-Structure Interaction module object
    FSI::FluidStructureInteraction* FSI_ = NULL;

    //- Fluid-Fluid module object
    FF::FluidFluid* FF_ = NULL;

    // NOTE: Add here a pointer for your new module object

    // Timesteps

    //- Timestep dictated by preCICE
    double timestepPrecice_;

    //- Timestep used by the solver
    double timestepSolver_;

    //- Stored (fixed) timestep
    double timestepStored_;

    // Checkpointing

    //- Checkpointed time (value)
    Foam::scalar couplingIterationTimeValue_;

    //- Checkpointed time (index)
    Foam::label couplingIterationTimeIndex_;

    //- Checkpointed mesh points
    Foam::pointField meshPoints_;
    Foam::pointField oldMeshPoints_;
    bool meshCheckPointed = false;

    // TODO: Currently unused, see storeMeshPoints().
    //- Checkpointed mesh volume
    // bool oldVolsStored = false;
    // Foam::volScalarField::Internal * oldVols_;
    // Foam::volScalarField::Internal * oldOldVols_;
    // int curTimeIndex_ = 0;

    // Vectors of pointers to the checkpointed mesh fields and their copies


    //- Checkpointed vol mesh fields
    std::vector<Foam::volScalarField*> volScalarFields_;
    std::vector<Foam::volScalarField*> volScalarFieldCopies_;

    std::vector<Foam::volVectorField*> volVectorFields_;
    std::vector<Foam::volVectorField*> volVectorFieldCopies_;

    std::vector<Foam::volSymmTensorField*> volSymmTensorFields_;
    std::vector<Foam::volSymmTensorField*> volSymmTensorFieldCopies_;

    std::vector<Foam::volSphericalTensorField*>
        volSphericalTensorFields_;
    std::vector<Foam::volSphericalTensorField*>
        volSphericalTensorFieldCopies_;

    std::vector<Foam::volTensorField*> volTensorFields_;
    std::vector<Foam::volTensorField*> volTensorFieldCopies_;


    //- Checkpointed internal vol mesh fields
    std::vector<Foam::volScalarField::Internal*>
        volScalarInternalFields_;
    std::vector<Foam::volScalarField::Internal*>
        volScalarInternalFieldCopies_;

    std::vector<Foam::volVectorField::Internal*>
        volVectorInternalFields_;
    std::vector<Foam::volVectorField::Internal*>
        volVectorInternalFieldCopies_;

    std::vector<Foam::volSymmTensorField::Internal*>
        volSymmTensorInternalFields_;
    std::vector<Foam::volSymmTensorField::Internal*>
        volSymmTensorInternalFieldCopies_;

    std::vector<Foam::volSphericalTensorField::Internal*>
        volSphericalTensorInternalFields_;
    std::vector<Foam::volSphericalTensorField::Internal*>
        volSphericalTensorInternalFieldCopies_;

    std::vector<Foam::volTensorField::Internal*>
        volTensorInternalFields_;
    std::vector<Foam::volTensorField::Internal*>
        volTensorInternalFieldCopies_;


    //- Checkpointed surface mesh fields
    std::vector<Foam::surfaceScalarField*> surfaceScalarFields_;
    std::vector<Foam::surfaceScalarField*> surfaceScalarFieldCopies_;

    std::vector<Foam::surfaceVectorField*> surfaceVectorFields_;
    std::vector<Foam::surfaceVectorField*> surfaceVectorFieldCopies_;

    std::vector<Foam::surfaceSymmTensorField*> surfaceSymmTensorFields_;
    std::vector<Foam::surfaceSymmTensorField*> surfaceSymmTensorFieldCopies_;

    std::vector<Foam::surfaceSphericalTensorField*>
        surfaceSphericalTensorFields_;
    std::vector<Foam::surfaceSphericalTensorField*>
        surfaceSphericalTensorFieldCopies_;

    std::vector<Foam::surfaceTensorField*> surfaceTensorFields_;
    std::vector<Foam::surfaceTensorField*> surfaceTensorFieldCopies_;


    //- Checkpointed point fields
    std::vector<Foam::pointScalarField*> pointScalarFields_;
    std::vector<Foam::pointScalarField*> pointScalarFieldCopies_;

    std::vector<Foam::pointVectorField*> pointVectorFields_;
    std::vector<Foam::pointVectorField*> pointVectorFieldCopies_;

    std::vector<Foam::pointSymmTensorField*> pointSymmTensorFields_;
    std::vector<Foam::pointSymmTensorField*> pointSymmTensorFieldCopies_;

    std::vector<Foam::pointSphericalTensorField*>
        pointSphericalTensorFields_;
    std::vector<Foam::pointSphericalTensorField*>
        pointSphericalTensorFieldCopies_;

    std::vector<Foam::pointTensorField*> pointTensorFields_;
    std::vector<Foam::pointTensorField*> pointTensorFieldCopies_;


    // NOTE: Declare additional vectors for any other types required.

    // Configuration

    //- Read the adapter's configuration file
    bool configFileRead();

    //- Check the adapter's configuration file
    bool configFileCheck(const std::string adapterConfigFileName);

    // Methods communicating with preCICE

    //- Initialize preCICE and exchange the first data
    void initialize();

    //- Finalize and destroy preCICE
    void finalize();

    //- Advance preCICE
    void advance();

    //- Read the coupling data at each interface
    void readCouplingData();

    //- Write the coupling data at each interface
    void writeCouplingData();

    //- Adjust the timestep of the solver according to preCICE
    void adjustSolverTimeStep();

    //- Determine if the coupling is still happening
    bool isCouplingOngoing();

    //- Determine if the coupling timestep has been completed
    bool isCouplingTimeWindowComplete();

    //- Determine if a checkpoint must be read
    bool isReadCheckpointRequired();

    //- Determine if a checkpoint must be written
    bool isWriteCheckpointRequired();

    //- Tell preCICE that the checkpoint has been read
    void fulfilledReadCheckpoint();

    //- Tell preCICE that the checkpoint has been written
    void fulfilledWriteCheckpoint();

    // Methods for checkpointing

    //- Configure the mesh checkpointing
    void setupMeshCheckpointing();

    //- Configure the mesh checkpointing
    void setupMeshVolCheckpointing();

    //- Configure the checkpointing
    void setupCheckpointing();

    //- Make a copy of the runTime object
    void storeCheckpointTime();

    //- Restore the copy of the runTime object
    void reloadCheckpointTime();

    //- Store the locations of the mesh points
    void storeMeshPoints();

    //- Restore the locations of the mesh points
    void reloadMeshPoints();

    // Add mesh checkpoint fields, depending on the type
    template<class GeoField>
    void addCheckpointFields
    (
        std::vector<GeoField*>& fields,
        std::vector<GeoField*>& fieldCopies
    )
    {
        wordList names(mesh_.names(GeoField::typeName));
        forAll(names, i)
        {
            GeoField& fld = mesh_.lookupObjectRef<GeoField>(names[i]);
            fields.push_back(&fld);
            fieldCopies.push_back
            (
                new GeoField
                (
                    "preCICE:" + fld.name(),
                    fld
                )
            );

            adapterInfo("Checkpointing: " + fld.name());
        }
    }

    // NOTE: Read checkpoint fields
    template<class GeoField>
    void readCheckpointFields
    (
        std::vector<GeoField*>& fields,
        const std::vector<GeoField*>& fieldCopies
    )
    {
        // Reload all the fields of type GeoField
        for (uint i = 0; i < fields.size(); i++)
        {
            const GeoField& fldCopy(*fieldCopies.at(i));
            GeoField& fld(*fields.at(i));

            // Load the volume field
            fld == fldCopy;

            int nOldTimes(fld.nOldTimes());
            int oldI = 0;
            GeoField* fldOld = &fld;
            const GeoField* fldCopyOld = &fldCopy;
            while (++oldI <= nOldTimes)
            {
                fldOld = &fldOld->oldTime();
                fldCopyOld = &fldCopyOld->oldTime();

                *fldOld == *fldCopyOld;
            }
        }
    }

    // NOTE: Read checkpoint fields
    template<class GeoField>
    void readInternalCheckpointFields
    (
        std::vector<GeoField*>& fields,
        const std::vector<GeoField*>& fieldCopies
    )
    {
        // Reload all the fields of type GeoField
        for (uint i = 0; i < fields.size(); i++)
        {
            const GeoField& fldCpy(*fieldCopies.at(i));
            GeoField& fld(*fields.at(i));

            // Load the volume field
            fld == fldCpy;
        }
    }

    // Update checkpoint fields
    template<class GeoField>
    void writeCheckpointFields
    (
        const std::vector<GeoField*>& fields,
        std::vector<GeoField*>& fieldCopies
    )
    {
        // Reload all the fields of type GeoField
        for (uint i = 0; i < fields.size(); i++)
        {
            // Store the field
            GeoField& fldCpy(*fieldCopies.at(i));
            const GeoField& fld(*fields.at(i));

            // Load the volume field
            fldCpy == fld;
        }
    }

    // Update checkpoint fields
    template<class GeoField>
    void clearCheckpointFields
    (
        std::vector<GeoField*>& fieldCopies
    )
    {
        // Clear all copied fields
        for (uint i = 0; i < fieldCopies.size(); i++)
        {
            delete fieldCopies.at(i);
        }
    }

    // TODO Probably these can be included to the mesh checkpoints.
    //- Read the volume checkpoint - restore the mesh volume fields
    void readCheckpoint();

    //- Write the volume checkpoint to a buffer - restore the mesh volume fields
    void writeCheckpoint();

    //- Destroy the preCICE interface and delete the allocated
    //  memory in a proper way. Called by the destructor.
    void teardown();

public:
    // Methods called by the functionObject

    //- Constructor
    Adapter(const Foam::Time& runTime, const Foam::fvMesh& mesh);

    //- Setup the adapter's configuration
    void configure();

    //- Called by the functionObject's execute()
    void execute();

    //- Called by the functionObject's setTimeStep()
    void setTimeStep();

    //- Called by the functionObject's end()
    void end();

    double timestepPrecice() const
    {
        return timestepPrecice_;
    }

    //- Destructor
    ~Adapter();
};

}

#endif

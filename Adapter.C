#include "Adapter.H"
#include "Interface.H"
#include "Utilities.H"

#include "IOstreams.H"

using namespace Foam;

preciceAdapter::Adapter::Adapter(const Time& runTime, const fvMesh& mesh)
: runTime_(runTime),
  mesh_(mesh)
{
    adapterInfo("Loaded the OpenFOAM-preCICE adapter - v1.1.0.", "info");

    return;
}

bool preciceAdapter::Adapter::configFileRead()
{

    // We need a try-catch here, as if reading preciceDict fails,
    // the respective exception will be reduced to a warning.
    // See also comment in preciceAdapter::Adapter::configure().
    try
    {
        adapterInfo("Reading preciceDict...", "info");

        // TODO: static is just a quick workaround to be able
        // to find the dictionary also out of scope (e.g. in KappaEffective).
        // We need a better solution.
        static IOdictionary preciceDict(
            IOobject(
                "preciceDict",
                runTime_.system(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE));

        // Read and display the preCICE configuration file name
        preciceConfigFilename_ = static_cast<fileName>(preciceDict.lookup("preciceConfig"));
        DEBUG(adapterInfo("  precice-config-file : " + preciceConfigFilename_));

        // Read and display the participant name
        participantName_ = static_cast<word>(preciceDict.lookup("participant"));
        DEBUG(adapterInfo("  participant name    : " + participantName_));

        // Read and display the list of modules
        DEBUG(adapterInfo("  modules requested   : "));
        auto modules_ = static_cast<wordList>(preciceDict.lookup("modules"));
        for (const auto& module : modules_)
        {
            DEBUG(adapterInfo("  - " + module + "\n"));

            // Set the modules switches
            if (module == "CHT")
            {
                CHTenabled_ = true;
            }
            if (module == "FSI")
            {
                FSIenabled_ = true;
            }
            if (module == "FF")
            {
                FFenabled_ = true;
            }
        }

        // Every interface is a subdictionary of "interfaces",
        // each with an arbitrary name. Read all of them and create
        // a list (here: pointer) of dictionaries.
        const auto* interfaceDictPtr = preciceDict.subDictPtr("interfaces");
        DEBUG(adapterInfo("  interfaces : "));

        // Check if we found any interfaces
        // and get the details of each interface
        if (!interfaceDictPtr)
        {
            adapterInfo("  Empty list of interfaces", "warning");
            return false;
        }
        else
        {
            for (const entry& interfaceDictEntry : *interfaceDictPtr)
            {
                if (interfaceDictEntry.isDict())
                {
                    const dictionary& interfaceDict = interfaceDictEntry.dict();
                    struct InterfaceConfig interfaceConfig;

                    interfaceConfig.meshName = static_cast<word>(interfaceDict.lookup("mesh"));
                    DEBUG(adapterInfo("  - mesh         : " + interfaceConfig.meshName));

                    // By default, assume "faceCenters" as locationsType
                    interfaceConfig.locationsType = interfaceDict.lookupOrDefault<word>("locations", "faceCenters");
                    DEBUG(adapterInfo("    locations    : " + interfaceConfig.locationsType));

                    // By default, assume that no mesh connectivity is required (i.e. no nearest-projection mapping)
                    interfaceConfig.meshConnectivity = interfaceDict.lookupOrDefault<bool>("connectivity", false);
                    // Mesh connectivity only makes sense in case of faceNodes, check and raise a warning otherwise
                    if (interfaceConfig.meshConnectivity && interfaceConfig.locationsType == "faceCenters")
                    {
                        DEBUG(adapterInfo("Mesh connectivity is not supported for faceCenters. \n"
                                          "Please configure the desired interface with the locationsType faceNodes. \n"
                                          "Have a look in the adapter documentation for detailed information.",
                                          "warning"));
                        return false;
                    }
                    DEBUG(adapterInfo("    connectivity : " + std::to_string(interfaceConfig.meshConnectivity)));

                    DEBUG(adapterInfo("    patches      : "));
                    auto patches = static_cast<wordList>(interfaceDict.lookup("patches"));
                    for (auto patch : patches)
                    {
                        interfaceConfig.patchNames.push_back(patch);
                        DEBUG(adapterInfo("      - " + patch));
                    }

                    DEBUG(adapterInfo("    writeData    : "));
                    auto writeData = static_cast<wordList>(interfaceDict.lookup("writeData"));
                    for (auto writeDatum : writeData)
                    {
                        interfaceConfig.writeData.push_back(writeDatum);
                        DEBUG(adapterInfo("      - " + writeDatum));
                    }

                    DEBUG(adapterInfo("    readData     : "));
                    auto readData = static_cast<wordList>(interfaceDict.lookup("readData"));
                    for (auto readDatum : readData)
                    {
                        interfaceConfig.readData.push_back(readDatum);
                        DEBUG(adapterInfo("      - " + readDatum));
                    }
                    interfacesConfig_.push_back(interfaceConfig);
                }
            }
        }

        // NOTE: set the switch for your new module here

        // If the CHT module is enabled, create it, read the
        // CHT-specific options and configure it.
        if (CHTenabled_)
        {
            CHT_ = new CHT::ConjugateHeatTransfer(mesh_);
            if (!CHT_->configure(preciceDict))
            {
                return false;
            }
        }

        // If the FSI module is enabled, create it, read the
        // FSI-specific options and configure it.
        if (FSIenabled_)
        {
            // Check for unsupported FSI with meshConnectivity
            for (uint i = 0; i < interfacesConfig_.size(); i++)
            {
                if (interfacesConfig_.at(i).meshConnectivity == true)
                {
                    adapterInfo(
                        "Mesh connectivity is not supported for FSI, as, usually, "
                        "the Solid participant needs to provide the connectivity information. "
                        "Therefore, set provideMeshConnectivity = false. "
                        "Have a look in the adapter documentation for more information. ",
                        "warning");
                    return false;
                }
            }

            FSI_ = new FSI::FluidStructureInteraction(mesh_, runTime_);
            if (!FSI_->configure(preciceDict))
            {
                return false;
            }
        }

        if (FFenabled_)
        {
            FF_ = new FF::FluidFluid(mesh_);
            if (!FF_->configure(preciceDict))
            {
                return false;
            }
        }

        // NOTE: Create your module and read any options specific to it here

        if (!CHTenabled_ && !FSIenabled_ && !FFenabled_) // NOTE: Add your new switch here
        {
            adapterInfo("No module is enabled.", "error-deferred");
            return false;
        }

        // TODO: Loading modules should be implemented in more general way,
        // in order to avoid code duplication. See issue #16 on GitHub.
    }
    catch (const Foam::error& e)
    {
        adapterInfo(e.message(), "error-deferred");
        return false;
    }

    return true;
}

void preciceAdapter::Adapter::configure()
{
    // Read the adapter's configuration file
    if (!configFileRead())
    {
        // This method is called from the functionObject's read() method,
        // which is called by the Foam::functionObjectList::read() method.
        // All the exceptions triggered in this method are caught as
        // warnings and the simulation continues simply without the
        // functionObject. However, we want the simulation to exit with an
        // error in case something is wrong. We store the information that
        // there was an error and it will be handled by the first call to
        // the functionObject's execute(), which can throw errors normally.
        errorsInConfigure = true;

        return;
    }

    try
    {
        // Check the timestep type (fixed vs adjustable)
        DEBUG(adapterInfo("Checking the timestep type (fixed vs adjustable)..."));
        adjustableTimestep_ = runTime_.controlDict().lookupOrDefault("adjustTimeStep", false);

        if (adjustableTimestep_)
        {
            DEBUG(adapterInfo("  Timestep type: adjustable."));
        }
        else
        {
            DEBUG(adapterInfo("  Timestep type: fixed."));
        }

        // Initialize preCICE
        DEBUG(adapterInfo("Creating the preCICE solver interface..."));
        DEBUG(adapterInfo("  Number of processes: " + std::to_string(Pstream::nProcs())));
        DEBUG(adapterInfo("  MPI rank: " + std::to_string(Pstream::myProcNo())));
        precice_ = new precice::SolverInterface(participantName_, preciceConfigFilename_, Pstream::myProcNo(), Pstream::nProcs());
        DEBUG(adapterInfo("  preCICE solver interface was created."));

        // Create interfaces
        DEBUG(adapterInfo("Creating interfaces..."));
        for (uint i = 0; i < interfacesConfig_.size(); i++)
        {
            Interface* interface = new Interface(*precice_, mesh_, interfacesConfig_.at(i).meshName, interfacesConfig_.at(i).locationsType, interfacesConfig_.at(i).patchNames, interfacesConfig_.at(i).meshConnectivity);
            interfaces_.push_back(interface);
            DEBUG(adapterInfo("Interface created on mesh " + interfacesConfig_.at(i).meshName));

            DEBUG(adapterInfo("Adding coupling data writers..."));
            for (uint j = 0; j < interfacesConfig_.at(i).writeData.size(); j++)
            {
                std::string dataName = interfacesConfig_.at(i).writeData.at(j);

                unsigned int inModules = 0;

                // Add CHT-related coupling data writers
                if (CHTenabled_ && CHT_->addWriters(dataName, interface))
                {
                    inModules++;
                }

                // Add FSI-related coupling data writers
                if (FSIenabled_ && FSI_->addWriters(dataName, interface))
                {
                    inModules++;
                }

                // Add FF-related coupling data writers
                if (FFenabled_ && FF_->addWriters(dataName, interface))
                {
                    inModules++;
                }

                if (inModules == 0)
                {
                    adapterInfo("I don't know how to write \"" + dataName
                                    + "\". Maybe this is a typo or maybe you need to enable some adapter module?",
                                "error-deferred");
                }
                else if (inModules > 1)
                {
                    adapterInfo("It looks like more than one modules can write \"" + dataName
                                    + "\" and I don't know how to choose. Try disabling one of the modules.",
                                "error-deferred");
                }

                // NOTE: Add any coupling data writers for your module here.
            } // end add coupling data writers

            DEBUG(adapterInfo("Adding coupling data readers..."));
            for (uint j = 0; j < interfacesConfig_.at(i).readData.size(); j++)
            {
                std::string dataName = interfacesConfig_.at(i).readData.at(j);

                unsigned int inModules = 0;

                // Add CHT-related coupling data readers
                if (CHTenabled_ && CHT_->addReaders(dataName, interface)) inModules++;

                // Add FSI-related coupling data readers
                if (FSIenabled_ && FSI_->addReaders(dataName, interface)) inModules++;

                // Add FF-related coupling data readers
                if (FFenabled_ && FF_->addReaders(dataName, interface)) inModules++;

                if (inModules == 0)
                {
                    adapterInfo("I don't know how to read \"" + dataName
                                    + "\". Maybe this is a typo or maybe you need to enable some adapter module?",
                                "error-deferred");
                }
                else if (inModules > 1)
                {
                    adapterInfo("It looks like more than one modules can read \"" + dataName
                                    + "\" and I don't know how to choose. Try disabling one of the modules.",
                                "error-deferred");
                }

                // NOTE: Add any coupling data readers for your module here.
            } // end add coupling data readers

            // Create the interface's data buffer
            interface->createBuffer();
        }

        // Initialize preCICE and exchange the first coupling data
        initialize();

        // Read the received coupling data
        readCouplingData();

        // If checkpointing is required, specify the checkpointed fields
        // and write the first checkpoint
        if (isWriteCheckpointRequired())
        {
            checkpointing_ = true;

            // Setup the checkpointing (find and add fields to checkpoint)
            setupCheckpointing();

            // Write checkpoint (for the first iteration)
            writeCheckpoint();
            fulfilledWriteCheckpoint();
        }

        // Adjust the timestep for the first iteration
        adjustSolverTimeStep();


        // If the solver tries to end before the coupling is complete,
        // e.g. because the solver's endTime was smaller or (in implicit
        // coupling) equal with the max-time specified in preCICE,
        // problems may occur near the end of the simulation,
        // as the function object may be called only once near the end.
        // See the implementation of Foam::Time::run() for more details.
        // To prevent this, we set the solver's endTime to "infinity"
        // and let only preCICE control the end of the simulation.
        // This has the side-effect of not triggering the end() method
        // in any function object normally. Therefore, we trigger it
        // when preCICE dictates to stop the coupling.
        adapterInfo(
            "Setting the solver's endTime to infinity to prevent early exits. "
            "Only preCICE will control the simulation's endTime. "
            "Any functionObject's end() method will be triggered by the adapter. "
            "You may disable this behavior in the adapter's configuration.",
            "info");
        const_cast<Time&>(runTime_).setEndTime(GREAT);
    }
    catch (const Foam::error& e)
    {
        adapterInfo(e.message(), "error-deferred");
        errorsInConfigure = true;
    }

    return;
}

void preciceAdapter::Adapter::execute()
{
    if (errorsInConfigure)
    {
        // Handle any errors during configure().
        // See the comments in configure() for details.
        adapterInfo(
            "There was a problem while configuring the adapter. "
            "See the log for details.",
            "error");
    }

    // The solver has already solved the equations for this timestep.
    // Now call the adapter's methods to perform the coupling.

    // TODO add a function which checks if all fields are checkpointed.
    // if (ncheckpointed is nregisterdobjects. )

    // Write the coupling data in the buffer
    writeCouplingData();

    // Advance preCICE
    advance();

    // Read checkpoint if required
    if (isReadCheckpointRequired())
    {
        readCheckpoint();
        fulfilledReadCheckpoint();
    }

    // Adjust the timestep, if it is fixed
//     if (!adjustableTimestep_)
//     {
//         adjustSolverTimeStep();
//     }

    // Write checkpoint if required
    if (isWriteCheckpointRequired())
    {
        writeCheckpoint();
        fulfilledWriteCheckpoint();
    }

    // As soon as OpenFOAM writes the results, it will not try to write again
    // if the time takes the same value again. Therefore, during an implicit
    // coupling, we write again when the coupling timestep is complete.
    // Check the behavior e.g. by using watch on a result file:
    //     watch -n 0.1 -d ls --full-time Fluid/0.01/T.gz
    if (checkpointing_ && isCouplingTimeWindowComplete())
    {
        // Check if the time directory already exists
        // (i.e. the solver wrote results that need to be updated)
        if (runTime_.timePath().type() == fileType::directory)
        {
            adapterInfo(
                "The coupling timestep completed. "
                "Writing the updated results.",
                "info");
            const_cast<Time&>(runTime_).writeNow();
        }
    }

    // Read the received coupling data from the buffer
    readCouplingData();

    // If the coupling is not going to continue, tear down everything
    // and stop the simulation.
    if (!isCouplingOngoing())
    {
        adapterInfo("The coupling completed.", "info");

        // Finalize the preCICE solver interface and delete data
        finalize();

        // Tell OpenFOAM to stop the simulation.
        // Set the solver's endTime to now. The next evaluation of
        // runTime.run() will be false and the solver will exit.
        const_cast<Time&>(runTime_).setEndTime(runTime_.value());
        adapterInfo(
            "The simulation was ended by preCICE. "
            "Calling the end() methods of any functionObject explicitly.",
            "info");
        adapterInfo("Great that you are using the OpenFOAM-preCICE adapter! "
                    "Next to the preCICE library and any other components, please also cite this adapter. "
                    "Find how on https://precice.org/adapter-openfoam-overview.html.",
                    "info");
        const_cast<Time&>(runTime_).functionObjects().end();
    }

    return;
}

void preciceAdapter::Adapter::setTimeStep()
{
    adjustSolverTimeStep();

    return;
}

void preciceAdapter::Adapter::readCouplingData()
{
    DEBUG(adapterInfo("Reading coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->readCouplingData();
    }

    return;
}

void preciceAdapter::Adapter::writeCouplingData()
{
    DEBUG(adapterInfo("Writing coupling data..."));

    for (uint i = 0; i < interfaces_.size(); i++)
    {
        interfaces_.at(i)->writeCouplingData();
    }

    return;
}

void preciceAdapter::Adapter::initialize()
{
    DEBUG(adapterInfo("Initalizing the preCICE solver interface..."));
    timestepPrecice_ = precice_->initialize();

    preciceInitialized_ = true;

    if (precice_->isActionRequired(precice::constants::actionWriteInitialData()))
    {
        writeCouplingData();
        precice_->markActionFulfilled(precice::constants::actionWriteInitialData());
    }

    DEBUG(adapterInfo("Initializing preCICE data..."));
    precice_->initializeData();

    adapterInfo("preCICE was configured and initialized", "info");

    return;
}

void preciceAdapter::Adapter::finalize()
{
    if (NULL != precice_ && preciceInitialized_ && !isCouplingOngoing())
    {
        DEBUG(adapterInfo("Finalizing the preCICE solver interface..."));

        // Finalize the preCICE solver interface
        precice_->finalize();

        preciceInitialized_ = false;

        // Delete the solver interface and all the related data
        teardown();
    }
    else
    {
        adapterInfo("Could not finalize preCICE.", "error");
    }

    return;
}

void preciceAdapter::Adapter::advance()
{
    DEBUG(adapterInfo("Advancing preCICE..."));
    timestepPrecice_ = precice_->advance(runTime_.deltaTValue());
    return;
}

void preciceAdapter::Adapter::adjustSolverTimeStep()
{
    DEBUG(adapterInfo("Adjusting the solver's timestep..."));

    // The timestep size that the solver has determined that it wants to use
    double timestepSolverDetermined;

    /* In this method, the adapter overwrites the timestep used by OpenFOAM.
       If the timestep is not adjustable, OpenFOAM will not try to re-estimate
       the timestep or read it again from the controlDict. Therefore, store
       the value that the timestep has is the beginning and try again to use this
       in every iteration.
       // TODO Treat also the case where the user modifies the timestep
       // in the controlDict during the simulation.
    */

    // Is the timestep adjustable or fixed?
    if (!adjustableTimestep_)
    {
        // Have we already stored the timestep?
        if (!useStoredTimestep_)
        {
            // Show a warning if runTimeModifiable is set
            if (runTime_.runTimeModifiable())
            {
                adapterInfo(
                    "You have enabled 'runTimeModifiable' in the "
                    "controlDict. The preciceAdapter does not yet "
                    "fully support this functionality when "
                    "'adjustableTimestep' is not enabled. "
                    "If you modify the 'deltaT' in the controlDict "
                    "during the simulation, it will not be updated.",
                    "warning");
            }

            // Store the value
            timestepStored_ = runTime_.deltaT().value();

            // Ok, we stored it once, we will use this from now on
            useStoredTimestep_ = true;

            if (timestepSolver_ > timestepPrecice_)
            {
                // Reduce the solver time step to the precice time step size
                adapterInfo(
                    "The solver's timestep is greater than the "
                    "coupling timestep.");

                timestepStored_ = timestepPrecice_;
                timestepSolver_ = timestepPrecice_;

                const_cast<Time&>(runTime_).setDeltaTNoAdjust(timestepSolver_);
            }
        }

        // Use the stored timestep as the determined solver's timestep
        timestepSolver_ = timestepStored_;
    }
    else
    {
        // The timestep is adjustable, so OpenFOAM will modify it
        // and therefore we can use the updated value
        // The actual time step will be set by the timeToNextWrite
        // in the function object
        timestepSolver_ = runTime_.deltaTValue();
    }

    return;
}

bool preciceAdapter::Adapter::isCouplingOngoing()
{
    bool isCouplingOngoing = false;

    // If the coupling ends before the solver ends,
    // the solver would try to access this method again,
    // giving a segmentation fault if precice_
    // was not available.
    if (NULL != precice_)
    {
        isCouplingOngoing = precice_->isCouplingOngoing();
    }

    return isCouplingOngoing;
}

bool preciceAdapter::Adapter::isCouplingTimeWindowComplete()
{
    return precice_->isTimeWindowComplete();
}

bool preciceAdapter::Adapter::isReadCheckpointRequired()
{
    return precice_->isActionRequired(precice::constants::actionReadIterationCheckpoint());
}

bool preciceAdapter::Adapter::isWriteCheckpointRequired()
{
    return precice_->isActionRequired(precice::constants::actionWriteIterationCheckpoint());
}

void preciceAdapter::Adapter::fulfilledReadCheckpoint()
{
    precice_->markActionFulfilled(precice::constants::actionReadIterationCheckpoint());

    return;
}

void preciceAdapter::Adapter::fulfilledWriteCheckpoint()
{
    precice_->markActionFulfilled(precice::constants::actionWriteIterationCheckpoint());

    return;
}

void preciceAdapter::Adapter::storeCheckpointTime()
{
    couplingIterationTimeIndex_ = runTime_.timeIndex();
    couplingIterationTimeValue_ = runTime_.value();
    DEBUG(adapterInfo("Stored time value t = " + std::to_string(runTime_.value())));

    return;
}

void preciceAdapter::Adapter::reloadCheckpointTime()
{
    const_cast<Time&>(runTime_).setTime(couplingIterationTimeValue_, couplingIterationTimeIndex_);
    // TODO also reset the current iteration?!
    DEBUG(adapterInfo("Reloaded time value t = " + std::to_string(runTime_.value())));

    return;
}

void preciceAdapter::Adapter::storeMeshPoints()
{
    DEBUG(adapterInfo("Storing mesh points..."));
    // TODO: In foam-extend, we would need "allPoints()". Check if this gives the same data.
    meshPoints_ = mesh_.points();
    oldMeshPoints_ = mesh_.oldPoints();

    /*
    // TODO  This is only required for subcycling. It should not be called when not subcycling!!
    // Add a bool 'subcycling' which can be evaluated every timestep.
    if ( !oldVolsStored && mesh_.foundObject<volScalarField::Internal>("V00") ) // For Ddt schemes which use one previous timestep
    {
        setupMeshVolCheckpointing();
        oldVolsStored = true;
    }
    // Update any volume fields from the buffer to the checkpointed values (if already exists.)
    */

    DEBUG(adapterInfo("Stored mesh points."));
    if (mesh_.moving())
    {
        if (!meshCheckPointed)
        {
            // Set up the checkpoint for the mesh flux: meshPhi
            setupMeshCheckpointing();
            meshCheckPointed = true;
        }
    }
}

void preciceAdapter::Adapter::reloadMeshPoints()
{
    // In Foam::polyMesh::movePoints.
    // TODO: The function movePoints overwrites the pointer to the old mesh.
    // Therefore, if you revert the mesh, the oldpointer will be set to the points, which are the new values.
    DEBUG(adapterInfo("Moving mesh points to their previous locations..."));

    // TODO
    // Switch oldpoints on for pure physics. (is this required?). Switch off for better mesh deformation capabilities?
    const_cast<pointField&>(mesh_.oldPoints()) = oldMeshPoints_;
    const_cast<fvMesh&>(mesh_).movePoints(meshPoints_);

    DEBUG(adapterInfo("Moved mesh points to their previous locations."));
}

void preciceAdapter::Adapter::setupMeshCheckpointing()
{
    // The other mesh <type>Fields:
    //      C
    //      Cf
    //      Sf
    //      magSf
    //      delta
    // are updated by the function fvMesh::movePoints. Only the meshPhi needs checkpointing.
    DEBUG(adapterInfo("Creating a list of the mesh checkpointed fields..."));

    // Add meshPhi to the checkpointed fields
//     surfaceScalarFields_.push_back
//     (
//         &const_cast<surfaceScalarField&>(mesh_.phi())
//     );
//     surfaceScalarFieldCopies_.push_back
//     (
//         new surfaceScalarField
//         (
//             "preCICE:" + mesh_.phi().name(),
//             mesh_.phi()
//         )
//     );
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.phi().name() + " in the list of checkpointed fields.");
#endif
}

void preciceAdapter::Adapter::setupMeshVolCheckpointing()
{
    DEBUG(adapterInfo("Creating a list of the mesh volume checkpointed fields..."));
    // Add the V0 and the V00 to the list of checkpointed fields.
    // For V0
//     volScalarInternalFields_.push_back
//     (
//         &const_cast<volScalarField::Internal&>(mesh_.V0())
//     );
//     volScalarInternalFieldCopies_.push_back
//     (
//         new volScalarField::Internal
//         (
//             "preCICE:" + mesh_.V0().name(),
//             mesh_.V0()
//         )
//     );
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V0().name() + " in the list of checkpointed fields.");
#endif
    // For V00
//     volScalarInternalFields_.push_back
//     (
//         &const_cast<volScalarField::Internal&>(mesh_.V00())
//     );
//     volScalarInternalFieldCopies_.push_back
//     (
//         new volScalarField::Internal
//         (
//             "preCICE:" + mesh_.V00().name(),
//             mesh_.V00()
//         )
//     );
#ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Added " + mesh_.V00().name() + " in the list of checkpointed fields.");
#endif
}


void preciceAdapter::Adapter::setupCheckpointing()
{
    // Add fields in the checkpointing list - sorted for parallel consistency
    DEBUG(adapterInfo("Adding in checkpointed fields..."));

    // Internal vol fields
    addCheckpointFields
    (
        volScalarInternalFields_,
        volScalarInternalFieldCopies_
    );
    addCheckpointFields
    (
        volVectorInternalFields_,
        volVectorInternalFieldCopies_
    );
    addCheckpointFields
    (
        volSymmTensorInternalFields_,
        volSymmTensorInternalFieldCopies_
    );
    addCheckpointFields
    (
        volSphericalTensorInternalFields_,
        volSphericalTensorInternalFieldCopies_
    );
    addCheckpointFields
    (
        volTensorInternalFields_,
        volTensorInternalFieldCopies_
    );

    // Vol fields
    addCheckpointFields
    (
        volScalarFields_,
        volScalarFieldCopies_
    );
    addCheckpointFields
    (
        volVectorFields_,
        volVectorFieldCopies_
    );
    addCheckpointFields
    (
        volSymmTensorFields_,
        volSymmTensorFieldCopies_
    );
    addCheckpointFields
    (
        volSphericalTensorFields_,
        volSphericalTensorFieldCopies_
    );
    addCheckpointFields
    (
        volTensorFields_,
        volTensorFieldCopies_
    );

    // Surface fields
    addCheckpointFields
    (
        surfaceScalarFields_,
        surfaceScalarFieldCopies_
    );
    addCheckpointFields
    (
        surfaceVectorFields_,
        surfaceVectorFieldCopies_
    );
    addCheckpointFields
    (
        surfaceSymmTensorFields_,
        surfaceSymmTensorFieldCopies_
    );
    addCheckpointFields
    (
        surfaceSphericalTensorFields_,
        surfaceSphericalTensorFieldCopies_
    );
    addCheckpointFields
    (
        surfaceTensorFields_,
        surfaceTensorFieldCopies_
    );

    // Point fields
    addCheckpointFields
    (
        pointScalarFields_,
        pointScalarFieldCopies_
    );
    addCheckpointFields
    (
        pointVectorFields_,
        pointVectorFieldCopies_
    );
    addCheckpointFields
    (
        pointSymmTensorFields_,
        pointSymmTensorFieldCopies_
    );
    addCheckpointFields
    (
        pointSphericalTensorFields_,
        pointSphericalTensorFieldCopies_
    );
    addCheckpointFields
    (
        pointTensorFields_,
        pointTensorFieldCopies_
    );
}


// NOTE: Add here methods to add other object types to checkpoint, if needed.
void preciceAdapter::Adapter::readCheckpoint()
{
    // TODO: To increase efficiency: only the oldTime() fields of the quantities which are used in the time
    //  derivative are necessary. (In general this is only the velocity). Also old information of the mesh
    //  is required.
    //  Therefore, loading the oldTime() and oldTime().oldTime() fields for the other fields can be excluded
    //  for efficiency.
    adapterInfo("Reading a checkpoint...");

    // Reload the runTime
    reloadCheckpointTime();

    // Reload the meshPoints (if FSI is enabled)
    if (FSIenabled_)
    {
        reloadMeshPoints();
    }

    // Vol fields
    readCheckpointFields
    (
        volScalarFields_,
        volScalarFieldCopies_
    );
    readCheckpointFields
    (
        volVectorFields_,
        volVectorFieldCopies_
    );
    readCheckpointFields
    (
        volSymmTensorFields_,
        volSymmTensorFieldCopies_
    );
    readCheckpointFields
    (
        volSphericalTensorFields_,
        volSphericalTensorFieldCopies_
    );
    readCheckpointFields
    (
        volTensorFields_,
        volTensorFieldCopies_
    );

    // Internal vol fields
    readInternalCheckpointFields
    (
        volScalarInternalFields_,
        volScalarInternalFieldCopies_
    );
    readInternalCheckpointFields
    (
        volVectorInternalFields_,
        volVectorInternalFieldCopies_
    );
    readInternalCheckpointFields
    (
        volSymmTensorInternalFields_,
        volSymmTensorInternalFieldCopies_
    );
    readInternalCheckpointFields
    (
        volSphericalTensorInternalFields_,
        volSphericalTensorInternalFieldCopies_
    );
    readInternalCheckpointFields
    (
        volTensorInternalFields_,
        volTensorInternalFieldCopies_
    );

    // Surface fields
    readCheckpointFields
    (
        surfaceScalarFields_,
        surfaceScalarFieldCopies_
    );
    readCheckpointFields
    (
        surfaceVectorFields_,
        surfaceVectorFieldCopies_
    );
    readCheckpointFields
    (
        surfaceSymmTensorFields_,
        surfaceSymmTensorFieldCopies_
    );
    readCheckpointFields
    (
        surfaceSphericalTensorFields_,
        surfaceSphericalTensorFieldCopies_
    );
    readCheckpointFields
    (
        surfaceTensorFields_,
        surfaceTensorFieldCopies_
    );

    // Point fields
    readCheckpointFields
    (
        pointScalarFields_,
        pointScalarFieldCopies_
    );
    readCheckpointFields
    (
        pointVectorFields_,
        pointVectorFieldCopies_
    );
    readCheckpointFields
    (
        pointSymmTensorFields_,
        pointSymmTensorFieldCopies_
    );
    readCheckpointFields
    (
        pointSphericalTensorFields_,
        pointSphericalTensorFieldCopies_
    );
    readCheckpointFields
    (
        pointTensorFields_,
        pointTensorFieldCopies_
    );

    // NOTE: Add here other field types to read, if needed.

// #ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint was read. Time = " + std::to_string(runTime_.value()));
// #endif

    return;
}


void preciceAdapter::Adapter::writeCheckpoint()
{
    adapterInfo("Writing a checkpoint...");

    // Store the runTime
    storeCheckpointTime();

    // Store the meshPoints (if FSI is enabled)
    if (FSIenabled_)
    {
        storeMeshPoints();
    }

    // Vol fields
    writeCheckpointFields
    (
        volScalarFields_,
        volScalarFieldCopies_
    );
    writeCheckpointFields
    (
        volVectorFields_,
        volVectorFieldCopies_
    );
    writeCheckpointFields
    (
        volSymmTensorFields_,
        volSymmTensorFieldCopies_
    );
    writeCheckpointFields
    (
        volSphericalTensorFields_,
        volSphericalTensorFieldCopies_
    );
    writeCheckpointFields
    (
        volTensorFields_,
        volTensorFieldCopies_
    );

    // Internal vol fields
    writeCheckpointFields
    (
        volScalarInternalFields_,
        volScalarInternalFieldCopies_
    );
    writeCheckpointFields
    (
        volVectorInternalFields_,
        volVectorInternalFieldCopies_
    );
    writeCheckpointFields
    (
        volSymmTensorInternalFields_,
        volSymmTensorInternalFieldCopies_
    );
    writeCheckpointFields
    (
        volSphericalTensorInternalFields_,
        volSphericalTensorInternalFieldCopies_
    );
    writeCheckpointFields
    (
        volTensorInternalFields_,
        volTensorInternalFieldCopies_
    );

    // Surface fields
    writeCheckpointFields
    (
        surfaceScalarFields_,
        surfaceScalarFieldCopies_
    );
    writeCheckpointFields
    (
        surfaceVectorFields_,
        surfaceVectorFieldCopies_
    );
    writeCheckpointFields
    (
        surfaceSymmTensorFields_,
        surfaceSymmTensorFieldCopies_
    );
    writeCheckpointFields
    (
        surfaceSphericalTensorFields_,
        surfaceSphericalTensorFieldCopies_
    );
    writeCheckpointFields
    (
        surfaceTensorFields_,
        surfaceTensorFieldCopies_
    );

    // Point fields
    writeCheckpointFields
    (
        pointScalarFields_,
        pointScalarFieldCopies_
    );
    writeCheckpointFields
    (
        pointVectorFields_,
        pointVectorFieldCopies_
    );
    writeCheckpointFields
    (
        pointSymmTensorFields_,
        pointSymmTensorFieldCopies_
    );
    writeCheckpointFields
    (
        pointSphericalTensorFields_,
        pointSphericalTensorFieldCopies_
    );
    writeCheckpointFields
    (
        pointTensorFields_,
        pointTensorFieldCopies_
    );

// #ifdef ADAPTER_DEBUG_MODE
    adapterInfo(
        "Checkpoint for time t = " + std::to_string(runTime_.value()) + " was stored.");
// #endif

    return;
}


void preciceAdapter::Adapter::end()
{
    // Throw a warning if the simulation exited before the coupling was complete
    if (NULL != precice_ && isCouplingOngoing())
    {
        adapterInfo("The solver exited before the coupling was complete.", "warning");
    }

    return;
}

void preciceAdapter::Adapter::teardown()
{
    // If the solver interface was not deleted before, delete it now.
    // Normally it should be deleted when isCouplingOngoing() becomes false.
    if (NULL != precice_)
    {
        DEBUG(adapterInfo("Destroying the preCICE solver interface..."));
        delete precice_;
        precice_ = NULL;
    }

    // Delete the preCICE solver interfaces
    if (interfaces_.size() > 0)
    {
        DEBUG(adapterInfo("Deleting the interfaces..."));
        for (uint i = 0; i < interfaces_.size(); i++)
        {
            delete interfaces_.at(i);
        }
        interfaces_.clear();
    }

    // Delete the copied fields for checkpointing
    if (checkpointing_)
    {
        DEBUG(adapterInfo("Deleting the checkpoints... "));

        // Vol fields
        clearCheckpointFields(volScalarFieldCopies_);
        clearCheckpointFields(volVectorFieldCopies_);
        clearCheckpointFields(volSymmTensorFieldCopies_);
        clearCheckpointFields(volSphericalTensorFieldCopies_);
        clearCheckpointFields(volTensorFieldCopies_);

        // Internal vol fields
        clearCheckpointFields(volScalarInternalFieldCopies_);
        clearCheckpointFields(volVectorInternalFieldCopies_);
        clearCheckpointFields(volSymmTensorInternalFieldCopies_);
        clearCheckpointFields(volSphericalTensorInternalFieldCopies_);
        clearCheckpointFields(volTensorInternalFieldCopies_);

        // Surface fields
        clearCheckpointFields(surfaceScalarFieldCopies_);
        clearCheckpointFields(surfaceVectorFieldCopies_);
        clearCheckpointFields(surfaceSymmTensorFieldCopies_);
        clearCheckpointFields(surfaceSphericalTensorFieldCopies_);
        clearCheckpointFields(surfaceTensorFieldCopies_);

        // Point fields
        clearCheckpointFields(pointScalarFieldCopies_);
        clearCheckpointFields(pointVectorFieldCopies_);
        clearCheckpointFields(pointSymmTensorFieldCopies_);
        clearCheckpointFields(pointSphericalTensorFieldCopies_);
        clearCheckpointFields(pointTensorFieldCopies_);

        // NOTE: Add here delete for other types, if needed

        checkpointing_ = false;
    }

    // Delete the CHT module
    if (NULL != CHT_)
    {
        DEBUG(adapterInfo("Destroying the CHT module..."));
        delete CHT_;
        CHT_ = NULL;
    }

    // Delete the FSI module
    if (NULL != FSI_)
    {
        DEBUG(adapterInfo("Destroying the FSI module..."));
        delete FSI_;
        FSI_ = NULL;
    }

    // Delete the FF module
    if (NULL != FF_)
    {
        DEBUG(adapterInfo("Destroying the FF module..."));
        delete FF_;
        FF_ = NULL;
    }

    // NOTE: Delete your new module here

    return;
}

preciceAdapter::Adapter::~Adapter()
{
    teardown();

    return;
}

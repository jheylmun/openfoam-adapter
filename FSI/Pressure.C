#include "Pressure.H"
#include "fluidThermo.H"
#include "facePointInterpolation.H"

using namespace Foam;


preciceAdapter::FSI::Pressure::Pressure(
    const Foam::fvMesh& mesh,
    const std::string solverType,
    const bool usePoint)
: mesh_(mesh),
  solverType_(solverType),
  osAreaPtr_(nullptr),
  osPressurePtr_(nullptr),
  osPointPressurePtr_(nullptr),
  pointPressure_(usePoint)
{
    //What about type "basic"?
    if
    (
        solverType_.compare("incompressible") != 0
     && solverType_.compare("compressible") != 0
    )
    {
        FatalErrorInFunction
            << "Force based calculations only support "
            << "compressible or incompressible solver types."
            << exit(FatalError);
    }

    dataType_ = scalar;

    if (logToFile)
    {
        osAreaPtr_.set(new OFstream("log.area"));
        osPressurePtr_.set(new OFstream("log.pressure"));
        if (usePoint)
        {
            osPointPressurePtr_.set(new OFstream("log.pointPressure"));
        }
        LogListHeader(osAreaPtr_(), false);
        LogListHeader(osPressurePtr_(), true);
        LogListHeader(osPointPressurePtr_(), false);
    }
}

std::string preciceAdapter::FSI::Pressure::getDataName() const
{
    return "Pressure";
}

//lookup correct rho
Foam::tmp<Foam::volScalarField> preciceAdapter::FSI::Pressure::rho() const
{
    // If volScalarField exists, read it from registry (for compressible cases)
    // interFoam is incompressible but has volScalarField rho

    if (mesh_.foundObject<volScalarField>("rho"))
    {
        return mesh_.lookupObject<volScalarField>("rho");
    }
    else if (solverType_.compare("incompressible") == 0)
    {
        const dictionary& FSIDict =
            mesh_.lookupObject<IOdictionary>("preciceDict").subOrEmptyDict("FSI");

        return tmp<volScalarField>(
            new volScalarField(
                IOobject(
                    "rho",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh_,
                dimensionedScalar(static_cast<dimensionedScalar>(FSIDict.lookup("rho")))));
    }
    else
    {
        FatalErrorInFunction
            << "Did not find the correct rho."
            << exit(FatalError);

        return volScalarField::null();
    }
}

void preciceAdapter::FSI::Pressure::write
(
    double* buffer,
    bool meshConnectivity,
    const unsigned int dim
)
{
    // Pressure boundary field
    tmp<volScalarField> tpressure;
    if (solverType_.compare("incompressible") == 0)
    {
        tpressure = mesh_.lookupObject<volScalarField>("p")*rho();
    }
    else
    {
        tpressure = tmp<volScalarField>
        (
            mesh_.lookupObject<volScalarField>("p")
        );
    }
    const volScalarField& pressure = tpressure();


    DynamicList<Foam::scalar> mappedArea;
    DynamicList<Foam::scalar> mappedPressure;

    if (log())
    {
        for (const label patchID : patchIDs_)
        {
            mappedArea.append(mesh_.magSf().boundaryField()[patchID]);
            mappedPressure.append(pressure.boundaryField()[patchID]);
        }
    }

    if (logToTerminal)
    {
        printListStatistics("Area", mappedArea);
        printListStatistics("Pressure", mappedPressure, mappedArea, true);
    }

    if (logToFile)
    {
        LogListStatistics
        (
            osAreaPtr_(),
            mesh_.time(),
            mappedArea
        );
        LogListStatistics
        (
            osPressurePtr_(),
            mesh_.time(),
            mappedPressure,
            mappedArea,
            true
        );
    }

    if (pointPressure_)
    {
        DynamicList<Foam::scalar> mappedPointPressure;
        const facePointInterpolation& interp =
            facePointInterpolation::New(mesh_);
        for (const label patchID : patchIDs_)
        {
            Field<Foam::scalar> pPressure
            (
                interp.interpolate(patchID, pressure.boundaryField()[patchID])
            );
            forAll(pPressure, i)
            {
                buffer[i] = pPressure[i];
            }
            if (log())
            {
                mappedPointPressure.append(pPressure);
            }
        }
        if (logToTerminal)
                printListStatistics("PointPressure", mappedPointPressure);
        if (logToFile)
        {
            LogListStatistics
            (
                osPointPressurePtr_(),
                mesh_.time(),
                mappedPointPressure
            );
        }
    }
    else
    {
        for (const label patchID : patchIDs_)
        {
            const scalarField& pPressure = pressure.boundaryField()[patchID];
            forAll(pPressure, i)
            {
                buffer[i] = pPressure[i];
            }
        }
    }
}

void preciceAdapter::FSI::Pressure::read
(
    double* buffer,
    const unsigned int dim
)
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE readBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Reading pressures is not supported."
        << exit(FatalError);
}

bool preciceAdapter::FSI::Pressure::isLocationTypeSupported(const bool meshConnectivity) const
{
    return
        this->locationType_ == LocationType::faceCenters
     || this->locationType_ == LocationType::faceNodes;
}

preciceAdapter::FSI::Pressure::~Pressure()
{}

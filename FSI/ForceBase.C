#include "ForceBase.H"
#include "fluidThermo.H"
#include "volPointInterpolation.H"

using namespace Foam;


preciceAdapter::FSI::ForceBase::ForceBase(
    const Foam::fvMesh& mesh,
    const std::string solverType,
    const std::string forceName,
    const dimensionSet forceDims,
    const bool usePoint)
: mesh_(mesh),
  solverType_(solverType),
  Force_(nullptr),
  pointForce_(nullptr)
{
    //What about type "basic"?
    if
    (
        solverType_.compare("incompressible") != 0
     && solverType_.compare("compressible") != 0
#ifdef SOLID
     && solverType_.compare("solid") != 0
#endif
    )
    {
        FatalErrorInFunction
            << "Force based calculations only support "
#ifdef SOLID
            << "solid, compressible, or incompressible solver types."
#else
            << "compressible or incompressible solver types."
#endif
            << exit(FatalError);
    }

    dataType_ = vector;

    Force_ = new volVectorField(
        IOobject(
            forceName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector(
            "fdim",
            forceDims,
            Foam::vector::zero));
    if (usePoint)
    {
        pointForce_ = new pointVectorField(
            IOobject(
                "point" + forceName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            pointMesh::New(mesh_),
            dimensionedVector(
                "fdim",
                forceDims,
                Foam::vector::zero));
    }
}

//Calculate viscous force
Foam::tmp<Foam::volSymmTensorField> preciceAdapter::FSI::ForceBase::devRhoReff() const
{
    //For turbulent flows
    typedef compressibleMomentumTransportModel cmpTurbModel;
    typedef incompressibleMomentumTransportModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::typeName))
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::typeName);

        return turb.devTau();
    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::typeName))
    {
        const icoTurbModel& turb =
            mesh_.lookupObject<icoTurbModel>(icoTurbModel::typeName);

        return rho() * turb.devSigma();
    }
    else
    {
        // For laminar flows get the velocity
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        return -mu() * dev(twoSymm(fvc::grad(U)));
    }
}

//lookup correct rho
Foam::tmp<Foam::volScalarField> preciceAdapter::FSI::ForceBase::rho() const
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

//lookup correct mu
Foam::tmp<Foam::volScalarField> preciceAdapter::FSI::ForceBase::mu() const
{

    if (solverType_.compare("incompressible") == 0)
    {
        if (mesh_.foundObject<fluidThermo>(basicThermo::dictName))
        {
            const fluidThermo& thermo =
                mesh_.lookupObject<fluidThermo>(basicThermo::dictName);

            return thermo.mu();
        }
        else
        {

            const dictionary& FSIDict =
                mesh_.lookupObject<IOdictionary>("preciceDict").subOrEmptyDict("FSI");

            dimensionedScalar nu(static_cast<dimensionedScalar>(FSIDict.lookup("nu")));

            return tmp<volScalarField>(
                new volScalarField(
                    nu * rho()));
        }
    }
    else if (solverType_.compare("compressible") == 0)
    {
        return mesh_.lookupObject<volScalarField>("thermo:mu");
    }
    else
    {
        FatalErrorInFunction
            << "Did not find the correct mu."
            << exit(FatalError);

        return volScalarField::null();
    }
}

void preciceAdapter::FSI::ForceBase::write
(
    double* buffer,
    bool meshConnectivity,
    const unsigned int dim
)
{
    // Compute forces. See the Forces function object.
    // Stress tensor boundary field
    tmp<volSymmTensorField> tdevRhoReff(devRhoReff());
    const volSymmTensorField::Boundary& devRhoReffb(
        tdevRhoReff().boundaryField());

    // Density boundary field
    tmp<volScalarField> trho(rho());
    const volScalarField::Boundary& rhob =
        trho().boundaryField();

    // Pressure boundary field
    const auto& pb = mesh_.lookupObject<volScalarField>("p").boundaryField();

    volVectorField::Boundary& bforce = Force_->boundaryFieldRef();
    // For every boundary patch of the interface
    for (const label patchID : patchIDs_)
    {
        tmp<vectorField> tsurface = getFaceVectors(patchID);
        const auto& surface = tsurface();

        // Pressure forces
        // FIXME: We need to substract the reference pressure for incompressible calculations
        if (solverType_.compare("incompressible") == 0)
        {
            bforce[patchID] = surface * pb[patchID] * rhob[patchID];
        }
        else if (solverType_.compare("compressible") == 0)
        {
            bforce[patchID] = surface * pb[patchID];
        }
        else
        {
            FatalErrorInFunction
                << "Forces calculation does only support "
                << "compressible or incompressible solver type."
                << exit(FatalError);
        }

        // Viscous forces
        bforce[patchID] += surface & devRhoReffb[patchID];
    }

    if (pointForce_)
    {
        volPointInterpolation::New(mesh_).interpolateBoundaryField
        (
            *Force_,
            *pointForce_
        );
        const pointVectorField::Boundary& bForce = pointForce_->boundaryField();
        for (const label patchID : patchIDs_)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchID];
            const labelList& meshPoints = patch.meshPoints();
            forAll(meshPoints, i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    buffer[i * dim + d] =
                        (*pointForce_)[meshPoints[i]][d];
                        // pForce[i][d];
            }
        }
    }
    else
    {
        for (const label patchID : patchIDs_)
        {
            forAll(bforce[patchID], i)
            {
                for (unsigned int d = 0; d < dim; ++d)
                    buffer[i * dim + d] =
                        Force_->boundaryField()[patchID][i][d];
            }
        }
    }
}

void preciceAdapter::FSI::ForceBase::readFromBuffer(double* buffer) const
{
    /* TODO: Implement
    * We need two nested for-loops for each patch,
    * the outer for the locations and the inner for the dimensions.
    * See the preCICE readBlockVectorData() implementation.
    */
    FatalErrorInFunction
        << "Reading forces is not supported."
        << exit(FatalError);
}

bool preciceAdapter::FSI::ForceBase::isLocationTypeSupported(const bool meshConnectivity) const
{
    return
        this->locationType_ == LocationType::faceCenters
     || this->locationType_ == LocationType::faceNodes;
}

preciceAdapter::FSI::ForceBase::~ForceBase()
{
    deleteDemandDrivenData(Force_);
    deleteDemandDrivenData(pointForce_);
}

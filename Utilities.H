#ifndef UTILITIES_H
#define UTILITIES_H

// Filter-out debug messages when not required
#ifdef ADAPTER_DEBUG_MODE
#define DEBUG(x) x
#else
#define DEBUG(x)
#endif

// String added in the beginning of every printed message
#define INFO_STR_ADAPTER "---[preciceAdapter] "

#include "IOstreams.H"
#include "Time.H"
#include "scalarField.H"
#include "vector.H"
#include "ops.H"

void adapterInfo(const std::string message, const std::string level = "debug");

void LogListHeader
(
    Foam::Ostream& os,
    const bool integrate = false
);

void printListStatistics
(
    const std::string& fieldName,
    const Foam::UList<Foam::scalar>& fld,
    const Foam::UList<Foam::scalar>& w = Foam::UList<Foam::scalar>(),
    const bool integrate = false,
    const std::string level = "info"
);
void printListStatistics
(
    const std::string& fieldName,
    const Foam::UList<Foam::vector>& fld,
    const Foam::UList<Foam::scalar>& w = Foam::UList<Foam::scalar>(),
    const bool integrate = false,
    const std::string level = "info"
);

void LogListStatistics
(
    Foam::Ostream& os,
    const Foam::Time& time,
    const Foam::UList<Foam::scalar>& fld,
    const Foam::UList<Foam::scalar>& w = Foam::UList<Foam::scalar>(),
    const bool integrate = false
);

void LogListStatistics
(
    Foam::Ostream& os,
    const Foam::Time& time,
    const Foam::UList<Foam::vector>& fld,
    const Foam::UList<Foam::scalar>& w = Foam::UList<Foam::scalar>(),
    const bool integrate = false
);

#endif

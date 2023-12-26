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
#include "ops.H"

void adapterInfo(const std::string message, const std::string level = "debug");

template<class T>
void printListStatistics
(
    const std::string& fieldName,
    const Foam::UList<T>& fld,
    const Foam::UList<Foam::scalar>& w = Foam::UList<Foam::scalar>(),
    const bool integrate = false,
    const std::string level = "debug"
)
{

    T fldAvg = Foam::Zero;
    T fldVar = Foam::Zero;
    Foam::scalar sumW = 0.0;
    bool weighted = false;
    if (Foam::returnReduce(w.size(), Foam::sumOp<Foam::scalar>()))
    {
        sumW = gSum(w);
        fldAvg = Foam::gSum(fld*w)/sumW;
        forAll(fld, i)
        {
            fldVar += Foam::cmptSqr(fld[i] - fldAvg)*w[i];
        }
        reduce(fldVar, Foam::sumOp<T>());
        fldVar /= sumW;
        weighted = true;
    }
    else
    {
        sumW = Foam::returnReduce(fld.size(), Foam::sumOp<Foam::label>());
        fldAvg = Foam::gSum(fld);
        if (sumW)
        {
            fldAvg /= sumW;
        }

        forAll(fld, i)
        {
            fldVar += Foam::cmptSqr(fld[i] - fldAvg);
        }
        Foam::reduce(fldVar, Foam::sumOp<T>());
        if (sumW)
        {
            fldVar /= sumW;
        }
    }

    Foam::OStringStream oss;
    oss <<  (
                weighted
              ? "Weighted statistics"
              : "Statistics"
            )
        << " for " << fieldName << '\n'
        << "    size (local/global) = " << fld.size() << " / "
        << Foam::returnReduce(fld.size(), Foam::sumOp<Foam::label>()) << '\n'
        << "    min = " << Foam::gMin(fld) << '\n'
        << "    minMagSqr = " << Foam::gMinMagSqr(fld) << '\n'
        << "    max = " << Foam::gMax(fld) << '\n'
        << "    maxMagSqr = " << Foam::gMaxMagSqr(fld) << '\n'
        << "    mean = " << fldAvg << '\n'
        << "    variance/std = " << fldVar << " / "
        << Foam::cmptPow(fldVar, Foam::pTraits<T>::one*0.5) << '\n';

    if (integrate && weighted)
    {
        oss << "    Integrated value = " << sumW*fldAvg << '\n';
    }

    if (fld.size())
    {
        adapterInfo(oss.str(), level);
    }
//     else
//     {
//         adapterInfo("Null", level);
//     }

    Foam::returnReduce(fld.size(), Foam::sumOp<Foam::label>());
}

inline void LogListHeader
(
    Foam::Ostream& os,
    const bool integrate = false
)
{
    if (Foam::Pstream::master())
    {
        os  << "# index" << '\n'
            << "# time" << '\n'
            << "# deltaT" << '\n'
            << "# min" << '\n'
            << "# minMagSqr" << '\n'
            << "# max" << '\n'
            << "# maxMagSqr" << '\n'
            << "# average" << '\n'
            << "# variance" << '\n';
        if (integrate)
        {
            os  << "# integral" << '\n';
        }
        os  << Foam::flush;
    }
}


template<class T>
void LogListStatistics
(
    Foam::Ostream& os,
    const Foam::Time& time,
    const Foam::UList<T>& fld,
    const Foam::UList<Foam::scalar>& w = Foam::UList<Foam::scalar>(),
    const bool integrate = false
)
{
    T fldMin = Foam::gMin(fld);
    T fldMinMagSqr = Foam::gMinMagSqr(fld);
    T fldMax = Foam::gMax(fld);
    T fldMaxMagSqr = Foam::gMaxMagSqr(fld);
    T fldAvg = Foam::Zero;
    T fldVar = Foam::Zero;
    Foam::scalar sumW = 1.0;
    bool validW = Foam::returnReduce(w.size(), Foam::sumOp<Foam::scalar>());
    if (validW)
    {
        sumW = gSum(w);
        fldAvg = Foam::gSum(fld*w)/sumW;
        forAll(fld, i)
        {
            fldVar += Foam::cmptSqr(fld[i] - fldAvg)*w[i];
        }
        reduce(fldVar, Foam::sumOp<T>());
        fldVar /= sumW;
    }
    else
    {
        sumW = Foam::returnReduce(fld.size(), Foam::sumOp<Foam::label>());
        fldAvg = Foam::gSum(fld);
        if (sumW)
        {
            fldAvg /= sumW;
        }

        forAll(fld, i)
        {
            fldVar += Foam::cmptSqr(fld[i] - fldAvg);
        }
        Foam::reduce(fldVar, Foam::sumOp<T>());
        if (sumW)
        {
            fldVar /= sumW;
        }
    }

    if (Foam::Pstream::master())
    {
        os  << time.timeIndex() << ' '
            << time.value() << ' '
            << time.deltaTValue() << ' '
            << fldMin << ' '
            << fldMinMagSqr << ' '
            << fldMax << ' '
            << fldMaxMagSqr << ' '
            << fldAvg << ' '
            << fldVar;
        if (integrate && validW)
        {
            os << ' ' << sumW*fldAvg;
        }
        os  << Foam::endl;
    }
}

#endif

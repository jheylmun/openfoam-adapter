#include "Utilities.H"

using namespace Foam;

void adapterInfo(const std::string message, const std::string level)
{
    if (level.compare("info") == 0)
    {
        // Prepend the message with a string
        Info << INFO_STR_ADAPTER
             << message.c_str()
             << nl;
    }
    else if (level.compare("warning") == 0)
    {
        // Produce a warning message with cyan header
        WarningInFunction
            << "\033[36m" // cyan color
            << "Warning in the preCICE adapter: "
            << "\033[0m" // restore color
            << nl
            << message.c_str()
            << nl
            << nl;
    }
    else if (level.compare("error") == 0)
    {
        // Produce an error message with red header
        // and exit the functionObject.
        // It will also exit the simulation, unless it
        // is called inside the functionObject's read().
        FatalErrorInFunction
            << "\033[31m" // red color
            << "Error in the preCICE adapter: "
            << "\033[0m" // restore color
            << nl
            << message.c_str()
            << nl
            << abort(FatalError);
    }
    else if (level.compare("error-deferred") == 0)
    {
        FatalError<<abort(FatalError);
        // Produce an warning message with red header.
        // OpenFOAM degrades errors inside read()
        // to warnings, stops the function object, but does
        // not exit. We throw a warning which is described
        // as an error, so that OpenFOAM does not exit,
        // but the user still sees that this is the actual
        // problem. We catch these errors and exit later.
        WarningInFunction
            << "\033[31m" // red color
            << "Error (deferred - will exit later) in the preCICE adapter: "
            << "\033[0m" // restore color
            << nl
            << message.c_str()
            << nl
            << nl;

    }
    else if (level.compare("debug") == 0)
    {
        Pout << INFO_STR_ADAPTER
             << "[DEBUG] "
             << message.c_str()
             << nl;
    }
    else if (level.compare("dev") == 0)
    {
        Info << "\033[35m" // cyan color
             << INFO_STR_ADAPTER
             << "[under development] "
             << "\033[0m " // restore color
             << message.c_str()
             << nl;
    }
    else
    {
        Info << INFO_STR_ADAPTER
             << "[unknown info level] "
             << message.c_str()
             << nl;
    }

    return;
}

void LogListHeader(Ostream& os, const bool integrate)
{
    if (Pstream::master())
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
        os  << flush;
    }
}

template<class T>
void printListStatistics
(
    const std::string& fieldName,
    const UList<T>& fld,
    const UList<scalar>& w,
    const bool integrate,
    const std::string level
)
{
    // Initialize statistics
    T fldAvg(Zero);
    T fldVar(Zero);

    scalar sumW = returnReduce(w.size(), sumOp<label>());

    // Weighted average and variance
    if (sumW)
    {
        sumW = gSum(w);
        fldAvg = gSum(fld*w)/sumW;
        forAll(fld, i)
        {
            fldVar += cmptSqr(fld[i] - fldAvg)*w[i];
        }
        reduce(fldVar, sumOp<T>());
        fldVar /= sumW;
    }

    // Number average and variance (w = 1/fld.size());
    else
    {
        sumW = returnReduce(fld.size(), sumOp<label>());
        fldAvg = gSum(fld)/sumW;

        forAll(fld, i)
        {
            fldVar += cmptSqr(fld[i] - fldAvg);
        }
        reduce(fldVar, sumOp<T>());
        fldVar /= sumW;
    }

    OStringStream oss;
    oss <<  (
                sumW
              ? "Weighted statistics"
              : "Statistics"
            )
        << " for " << fieldName << '\n'
        << "    size (local/global) = " << fld.size() << " / "
        << returnReduce(fld.size(), sumOp<label>()) << '\n'
        << "    min = " << gMin(fld) << '\n'
        << "    minMagSqr = " << gMinMagSqr(fld) << '\n'
        << "    max = " << gMax(fld) << '\n'
        << "    maxMagSqr = " << gMaxMagSqr(fld) << '\n'
        << "    mean = " << fldAvg << '\n'
        << "    variance/std = " << fldVar << " / "
        << cmptPow(fldVar, pTraits<T>::one*0.5) << '\n';

    if (integrate && sumW)
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
}

void printListStatistics
(
    const std::string& fieldName,
    const UList<scalar>& fld,
    const UList<scalar>& w,
    const bool integrate,
    const std::string level
)
{
    printListStatistics<scalar>(fieldName, fld, w, integrate, level);
}


void printListStatistics
(
    const std::string& fieldName,
    const UList<vector>& fld,
    const UList<scalar>& w,
    const bool integrate,
    const std::string level
)
{
    printListStatistics<vector>(fieldName, fld, w, integrate, level);
}


template<class T>
void LogListStatistics
(
    Ostream& os,
    const Time& time,
    const UList<T>& fld,
    const UList<scalar>& w,
    const bool integrate
)
{
    // Initialize statistics
    T fldMin(gMin(fld));
    T fldMinMagSqr(gMinMagSqr(fld));
    T fldMax(gMax(fld));
    T fldMaxMagSqr(gMaxMagSqr(fld));
    T fldAvg(Zero);
    T fldVar(Zero);

    // Check if weights are valid
    scalar sumW = returnReduce(w.size(), sumOp<label>());

    // Weighted average and variance
    if (sumW)
    {
        sumW = gSum(w);
        fldAvg = gSum(fld*w)/sumW;
        forAll(fld, i)
        {
            fldVar += cmptSqr(fld[i] - fldAvg)*w[i];
        }
        reduce(fldVar, sumOp<T>());
        fldVar /= sumW;
    }

    // Number average and variance (w = 1/fld.size())
    else
    {
        sumW = returnReduce(fld.size(), sumOp<label>());
        fldAvg = gSum(fld)/sumW;

        forAll(fld, i)
        {
            fldVar += cmptSqr(fld[i] - fldAvg);
        }
        reduce(fldVar, sumOp<T>());
        fldVar /= sumW;
    }

    if (Pstream::master())
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
        if (integrate && sumW)
        {
            os << ' ' << sumW*fldAvg;
        }
        os  << endl;
    }
}

void LogListStatistics
(
    Ostream& os,
    const Time& time,
    const UList<scalar>& fld,
    const UList<scalar>& w,
    const bool integrate
)
{
    LogListStatistics<scalar>(os, time, fld, w, integrate);
}

void LogListStatistics
(
    Ostream& os,
    const Time& time,
    const UList<vector>& fld,
    const UList<scalar>& w,
    const bool integrate
)
{
    LogListStatistics<vector>(os, time, fld, w, integrate);
}

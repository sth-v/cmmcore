//
// Created by Andrew Astakhov on 14.10.24.
//

#ifndef CMMCORE_HASH_H
#define CMMCORE_HASH_H
#include <unordered_set>
#ifdef CYTHON_ABI
#include "cm_limits.h"
#else
#include "cmmcore/cm_limits.h"
#endif

namespace cmmcore
{
    // Quantize to CMMCORE_DECIMALS decimal places
    static std::hash<double> doubleHash;

    inline int quantize(const double value, const int decimals)
    {
        return static_cast<int>(value * decimals);
    }

    namespace _hashimpl
    {
        struct PairDoubleHash
        {
            size_t operator()(const std::pair<double, double>& v) const
            {
                int qx = quantize(v.first, CMMCORE_DECIMALS);
                int qy = quantize(v.first,CMMCORE_DECIMALS);
                size_t h1 = std::hash<int>()(qx);
                size_t h2 = std::hash<int>()(qy);
                return h1 ^ (h2 << 1);
            }
        };

        struct PairDoubleEqual
        {
            bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b) const
            {
                return a.first == b.first && a.second == b.second;
            }
        };
    };

    using PairDoubleHash_t = _hashimpl::PairDoubleHash;
    using PairDoubleEqual_t = _hashimpl::PairDoubleEqual;

}
#endif //CMMCORE_HASH_H

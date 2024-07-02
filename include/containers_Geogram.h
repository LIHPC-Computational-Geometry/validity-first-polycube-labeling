#pragma once

#include <geogram/basic/numeric.h>      // for index_t
#include <geogram/basic/attributes.h>   // for Attribute
#include <geogram/basic/vecg.h>         // for vecng<>
#include <geogram/basic/geometry.h>     // for mat2

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <limits>  // for std::numeric_limits<>
#include <iomanip> // for std::setw()

using namespace GEO;

index_t min(const Attribute<index_t>& container);

// find the min element in a Geogram vector (vecng)
template <index_t DIM, class T>
T min(const vecng<DIM, T>& vector) {
    T min = std::numeric_limits<T>::min();
    FOR(d,DIM) {
        min = std::min(min,vector[d]);
    }
    return min;
}

// find the max element in a Geogram vector (vecng)
template <index_t DIM, class T>
T max(const vecng<DIM, T>& vector) {
    T max = std::numeric_limits<T>::min();
    FOR(d,DIM) {
        max = std::max(max,vector[d]);
    }
    return max;
}

// define how to print a GEO::mat2
std::ostream& operator<< (std::ostream &out, const GEO::mat2& data);
// Use the operator<< overloading with {fmt}
// https://fmt.dev/latest/api.html#std-ostream-support
template <> struct fmt::formatter<GEO::mat2> : ostream_formatter {};
#pragma once

#include <geogram/basic/numeric.h> // for index_t
#include <geogram/basic/attributes.h> // for Attribute

#include <algorithm>    // for std::find(), std::min_element()
#include <map>          // for std::map::find()

#define VECTOR_CONTAINS(container,value) (std::find(container.cbegin(),container.cend(),value) != container.cend())
#define MAP_CONTAINS(container,value) (container.find(value) != container.end())

#define VECTOR_MIN(container) (std::min_element(container.begin(),container.end())) // do not work on mesh attributes because they don't have iterators...

using namespace GEO;

index_t min(const Attribute<index_t>& container);

// to concatenate vectors
// https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
template<class T>
inline void operator +=(T& a, const T& b) {
    a.insert( a.end(), b.begin(), b.end() );
}

// to get the index of the last element in a vector
template<class T>
inline std::size_t index_of_last(const T& a) {
    return (a.size()-1);
}
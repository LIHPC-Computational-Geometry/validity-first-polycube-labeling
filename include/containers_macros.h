#pragma once

#include <algorithm>    // for std::find(), std::min_element(), std::max_element(), std::max()
#include <map>          // for std::map::find()
#include <iterator>     // for std::distance()
#include <numeric>      // for std::accumulate()

#define VECTOR_CONTAINS(container,value) (std::find(container.cbegin(),container.cend(),value) != container.cend())
#define MAP_CONTAINS(container,value) (container.find(value) != container.end())
#define INIT_LIST_CONTAINS(container,value) (std::find(container.begin(),container.end(),value) != container.end()) // there are no .cbegin()/.cend() for std::initializer_list

#define VECTOR_MIN(container) (std::min_element(container.begin(),container.end())) // do not work on mesh attributes because they don't have iterators...
#define VECTOR_MAX(container) (std::max_element(container.begin(),container.end())) // do not work on mesh attributes because they don't have iterators...
#define VECTOR_MIN_INDEX(container) (std::distance(container.begin(), VECTOR_MIN(container))) // get index of the min value
#define VECTOR_MAX_INDEX(container) (std::distance(container.begin(), VECTOR_MAX(container))) // get index of the max value

#define VECTOR_SUM(container) (std::accumulate(container.begin(),container.end(),0.0))
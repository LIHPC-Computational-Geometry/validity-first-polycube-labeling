#pragma once

#include <geogram/basic/numeric.h> // for index_t
#include <geogram/basic/attributes.h> // for Attribute
#include <geogram/basic/vecg.h> // for vecng<>

#include <algorithm>    // for std::find(), std::min_element(), std::max()
#include <map>          // for std::map::find()
#include <iterator>     // for std::iterator_traits, std::distance()
#include <numeric>      // for std::accumulate()
#include <limits>       // for std::numeric_limits<>

#define VECTOR_CONTAINS(container,value) (std::find(container.cbegin(),container.cend(),value) != container.cend())
#define MAP_CONTAINS(container,value) (container.find(value) != container.end())
#define INIT_LIST_CONTAINS(container,value) (std::find(container.begin(),container.end(),value) != container.end()) // there are no .cbegin()/.cend() for std::initializer_list

#define VECTOR_MIN(container) (std::min_element(container.begin(),container.end())) // do not work on mesh attributes because they don't have iterators...
#define VECTOR_MAX(container) (std::max_element(container.begin(),container.end())) // do not work on mesh attributes because they don't have iterators...
#define VECTOR_MIN_INDEX(container) (std::distance(container.begin(), VECTOR_MIN(container))) // get index of the min value
#define VECTOR_MAX_INDEX(container) (std::distance(container.begin(), VECTOR_MAX(container))) // get index of the max value

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

// to check is a vector has duplicate values
template<class T>
bool has_duplicates(const std::vector<T>& container) {
    std::set<T> values_encontered;
    for(const auto& value : container) {
        if(values_encontered.contains(value)) {
            return true;
        }
        values_encontered.insert(value);
    }
    return false;
}

// to compute the standard deviation
// thank you https://codereview.stackexchange.com/a/123278
template <typename It, 
    typename E = typename std::iterator_traits<It>::value_type, 
    typename R = typename std::common_type<double, E>::type>
R std_dev(It b, It e)
{
    R N          = (R) std::distance(b, e);
    R const mean = std::accumulate(b, e, R{}) / N;
    R variance   = std::accumulate(b, e, R{}, [mean](R a, E v)-> R { return a + (v-mean)*(v-mean); });
    return std::sqrt(variance / N);
}

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

// creation of a set of one element
// thank you Walter https://stackoverflow.com/a/37564479
template <typename T>
inline std::set<T> make_set(const T& x)
{ return {x}; }

// get map key where value is max
template <typename K, typename V>
K key_at_max_value(const std::map<K,V>& map) {
    // thanks Janek_Kozicki and cigien https://stackoverflow.com/a/54690905
    return std::max_element(map.begin(),map.end(),[] (const std::pair<K,V>& a, const std::pair<K,V>& b) -> bool{ return a.second < b.second; } )->first;
}

template <typename K, typename V>
void fill_set_with_map_keys(const std::map<K,V>& map, std::set<K>& keys) {
    keys.clear();
    for(auto [key,value] : map) {
        keys.insert(key);
    }
}

template <typename T>
bool no_item_in_common(const std::set<T>& a, const std::set<T>& b) {
    for(const auto& i : a) {
        if(b.contains(i)) {
            return false;
        }
    }
    for(const auto& j : b) {
        if(a.contains(j)) {
            return false;
        }
    }
    return true;
}
#pragma once

#include <iterator> // for std::iterator_traits
#include <vector>
#include <set>
#include <map>

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
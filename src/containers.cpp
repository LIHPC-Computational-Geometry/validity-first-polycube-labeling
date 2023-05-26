#include <algorithm>    // for std::min()

#include "containers.h"

index_t min(const Attribute<index_t>& container) {
    index_t min = max_index_t();
    FOR(index,container.size()) {
        min = std::min(min,container[index]);
    }
    return min;
}
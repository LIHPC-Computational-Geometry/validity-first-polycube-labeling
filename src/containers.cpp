#include <algorithm>    // for std::min()

#include "containers.h"

index_t min(const Attribute<index_t>& container) {
    index_t min = max_index_t();
    FOR(index,container.size()) {
        min = std::min(min,container[index]);
    }
    return min;
}

std::ostream& operator<< (std::ostream &out, const GEO::mat2& data) {
    out << std::setw(8) << data(0,0) << " " << data(0,1) << '\n' << 
           std::setw(8) << data(1,0) << " " << data(1,1) << std::endl;
    return out;
}
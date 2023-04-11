#pragma once

#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>

#include <iostream>
#include <vector>
#include <set>

using namespace GEO;

// A set of adjacent triangles having the same label
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L61
struct Chart {
    std::set<int> neighbors_ids; // IDs of adjacent charts
    int label = -1;

    void print() {
        std::cout << "\t label  " << label << std::endl;
        std::cout << "\t neighbors ";
        for (int i: neighbors_ids) std::cout << i << " ";
        std::cout << std::endl;
    }
};

// A labeling stored with charts, boundaries and corners
// based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/include/flagging.h#L135
class StaticLabelingGraph {

public:
    StaticLabelingGraph(Mesh& mesh, const char* facet_attribute);

    GEO::index_t get_nb_charts() const;

    const index_t* get_facet2chart_ptr() const;

private:
    std::vector<Chart> charts;
    std::vector<index_t> facet2chart;
    index_t nb_charts; // TODO replace by charts.size() when charts will be used
};
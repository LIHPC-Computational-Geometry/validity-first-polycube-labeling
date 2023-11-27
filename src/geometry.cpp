#include <queue>

#include "geometry.h"

void per_facet_distance(const Mesh& mesh, std::map<index_t,unsigned int>& distance) {
    // thank you Trizalio https://stackoverflow.com/a/72437022
    // Compute the distance to facets for which distance[facet]==0, for all facets in the keys of 'distance'
    std::queue<index_t> facets_to_explore; // facets we have to visit to update the distance in their neighborhood. FIFO data structure
    index_t current_facet = index_t(-1),
            adjacent_facet = index_t(-1);
    unsigned int current_distance = 0; // will store the distance currently computed for the current facet
    for(const auto& kv : distance) {
        if(kv.second == 0) {
            facets_to_explore.push(kv.first);
        }
    }
    while(!facets_to_explore.empty()) { // while there still are facets to explore
        current_facet = facets_to_explore.front(); // pick the facet at the front of the FIFO
        facets_to_explore.pop(); // remove it from the FIFO
        current_distance = distance[current_facet];
        FOR(le,3) { // for each local edge of the current facet
            adjacent_facet = mesh.facets.adjacent(current_facet,le); // get the facet index of the neighbor at the other side of the current local edge
            if(!distance.contains(adjacent_facet)) {
                continue; // -> do not explore the adjacent facet, it is not inside the set of facets
            }
            if (current_distance + 1 < distance[adjacent_facet]) { // if the distance of the neighbor was too much (passing by the current facet is closer)
                distance[adjacent_facet] = current_distance + 1; // update the distance
                facets_to_explore.emplace(adjacent_facet);
            }
        }
    }
}

bool facet_normals_are_inwards(Mesh& mesh) {
    // https://forums.cgsociety.org/t/check-if-mesh-is-inside-out/1688290

    // Find the vertex with the smallest X coordinate
    double smallest_X_coordinate = Numeric::max_float64();
    index_t corresponding_vertex = index_t(-1);
    FOR(v,mesh.vertices.nb()) { // for each vertex
        if(smallest_X_coordinate > mesh.vertices.point(v).x) {
            smallest_X_coordinate = mesh.vertices.point(v).x;
            corresponding_vertex = v;
        }
    }

    // find a facet adjacent to the corresponding_vertex
    // and create a halfedge whose origin is corresponding_vertex
    // Is there a easier way that going through all facets?
    MeshHalfedges::Halfedge halfedge;
    CustomMeshHalfedges mesh_he(mesh);
    FOR(f,mesh.facets.nb()) { // for each facet
        if(
        (mesh.facets.vertex(f,0) == corresponding_vertex) ||
        (mesh.facets.vertex(f,1) == corresponding_vertex) ||
        (mesh.facets.vertex(f,2) == corresponding_vertex) ) {
            halfedge.facet = f;
            break;
        }
    }
    halfedge.corner = mesh.facets.corner(halfedge.facet,0); // try the halfedge at local vertex 0
    while (Geom::halfedge_vertex_index_from(mesh,halfedge) != corresponding_vertex) {
        mesh_he.move_to_next_around_facet(halfedge);
    }

    // compute the vertex normal of the corresponding_vertex
    vec3 normal_of_corresponding_vertex;
    MeshHalfedges::Halfedge init_halfedge = halfedge;
    do {
        normal_of_corresponding_vertex += mesh_facet_normal(mesh,halfedge.facet);
        mesh_he.move_to_next_around_vertex(halfedge);
    } while (halfedge != init_halfedge);

    // compute the dot product of the vertex normal and a vector going outwards the mesh (towards -X)
    // if the dot product is negative, the normal is inwards
    return dot(normal_of_corresponding_vertex,vec3(-1.0,0.0,0.0)) < 0.0;
}

void flip_facet_normals(Mesh& mesh) {
    geo_assert(mesh.cells.nb() == 0); // must be a surface mesh
    geo_assert(mesh.facets.nb() != 0); // must have facets
    geo_assert(mesh.facets.are_simplices()); // must be a triangle mesh

    index_t tmp_vertex_index = index_t(-1);
    index_t tmp_facet_index = index_t(-1);
    FOR(f,mesh.facets.nb()) { // for each facet
        //
        // f = index of current facet
        // lv = local vertex in {0,1,2}
        // le = local edge in {0,1,2}
        // v = vertex index
        // af = adjacent facet
        //
        //                         vA
        //  ------------------------------------------------
        //   \                     / \                     /
        //    \                   /lv0\                   /
        //     \                 /     \                 /
        //      \     afK       /       \      afI      /
        //       \             /         \             /
        //        \           /           \           /
        //         \         /le2       le0\         /
        //          \       /      f        \       /
        //           \     /                 \     /
        //            \   /                   \   /
        //             \ /lv2     le1       lv1\ /
        //           vC ------------------------- vB
        //               \                     /
        //                \                   /
        //                 \                 /
        //                  \               /
        //                   \             /
        //                    \    afJ    /
        //                     \         /
        //                      \       /
        //                       \     /
        //                        \   /
        //                         \ /
        //
        // local vertices are clockwise -> right hand rule -> facet normal is inwards
        // swap lv1 and lv2 to flip the normal, and swap adjacent facet accordingly
        //
        //                         vA
        //  ------------------------------------------------
        //   \                     / \                     /
        //    \                   /lv0\                   /
        //     \                 /     \                 /
        //      \     afK       /       \      afI      /
        //       \             /         \             /
        //        \           /           \           /
        //         \         /le0       le2\         /
        //          \       /      f        \       /
        //           \     /                 \     /
        //            \   /                   \   /
        //             \ /lv1     le1       lv2\ /
        //           vC ------------------------- vB
        //               \                     /
        //                \                   /
        //                 \                 /
        //                  \               /
        //                   \             /
        //                    \    afJ    /
        //                     \         /
        //                      \       /
        //                       \     /
        //                        \   /
        //                         \ /
        //
        // local vertices are counterclockwise -> right hand rule -> facet normal is outwards
        //
        tmp_vertex_index = mesh.facets.vertex(f,1); // copy vB
        tmp_facet_index = mesh.facets.adjacent(f,0); // copy afI
        mesh.facets.set_vertex(f,1,mesh.facets.vertex(f,2));// at lv1 is no longer vB but vC
        mesh.facets.set_adjacent(f,0,mesh.facets.adjacent(f,2)); // beyond le0 is no longer afI but afK
        mesh.facets.set_vertex(f,2,tmp_vertex_index); // at lv2 is no longer vC but vB
        mesh.facets.set_adjacent(f,2,tmp_facet_index); // beyond le2 is no longer afK but afI
    }
}
#include <geogram/basic/memory.h> // for GEO::vector

#include <fmt/core.h>
#include <fmt/ostream.h>    // to use fmt::print() on ostreams

#include <queue>
#include <vector>
#include <utility>      // for std::pair, std::make_pair()
#include <algorithm>    // for std::min(), std::max()

#include "geometry.h"
#include "CustomMeshHalfedges.h"

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

bool facet_normals_are_inward(Mesh& mesh) {
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

    // compute the dot product of the vertex normal and a vector going outward the mesh (towards -X)
    // if the dot product is negative, the normal is inward
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
        // local vertices are clockwise -> right hand rule -> facet normals are inward
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
        // local vertices are counterclockwise -> right hand rule -> facet normals are outward
        //
        tmp_vertex_index = mesh.facets.vertex(f,1); // copy vB
        tmp_facet_index = mesh.facets.adjacent(f,0); // copy afI
        mesh.facets.set_vertex(f,1,mesh.facets.vertex(f,2));// at lv1 is no longer vB but vC
        mesh.facets.set_adjacent(f,0,mesh.facets.adjacent(f,2)); // beyond le0 is no longer afI but afK
        mesh.facets.set_vertex(f,2,tmp_vertex_index); // at lv2 is no longer vC but vB
        mesh.facets.set_adjacent(f,2,tmp_facet_index); // beyond le2 is no longer afK but afI
    }
}

void center_mesh(Mesh& mesh, bool normalize) {
    double xyzmin[3];
    double xyzmax[3];
    vec3 midpoint;
    double scale = Numeric::min_float64(); // will store width of widest dimension

    get_bbox(mesh, xyzmin, xyzmax);

    FOR(i,3) { // for each axis {0=X, 1=Y, 2=Z}
        if ( (xyzmax[i]-xyzmin[i]) > scale) {
            scale = xyzmax[i]-xyzmin[i];
        }
        midpoint[i] = (xyzmax[i]+xyzmin[i]) / 2.0;
    }

    if(normalize) {
        scale /= 2.0; // we need half of the bounding box length when rescaling
        FOR(v,mesh.vertices.nb()) {
            mesh.vertices.point(v) -= midpoint; // center
            mesh.vertices.point(v) /= scale; // scaling
        }
    }
    else {
        // center but no scaling
        FOR(v,mesh.vertices.nb()) {
            mesh.vertices.point(v) -= midpoint; // center
        }
    }
}

void compute_adjacent_facets_of_vertices(const Mesh& mesh, std::vector<std::vector<index_t>>& adj_facets) {
    adj_facets.clear();
    adj_facets.resize(mesh.vertices.nb());

    // link each facet to adjacent vertices
    FOR(f,mesh.facets.nb()) {
        FOR(lv,mesh.facets.nb_vertices(f)) { // for each local vertices of the facet f
            adj_facets[mesh.facets.vertex(f,lv)].push_back(f);
        }
    }
}

void remove_feature_edges_with_low_dihedral_angle(Mesh& mesh, std::vector<std::vector<index_t>>& adj_facets) {
    // Remove feature edges with low dihedral angle
    // We need facets at each side of edge e to compute angle between their normals
    // Issue : no adjecency between edges and facets
    // -> find an halfedge (oriented edge) which is on this edge. halfedges have adjacency with facets
    // issue : we need a facet index and a facet corner index to construct an halfedge
    // -> find the facet and the facet corner from the two vertices of the edge e
    // issue : no adjacency between vertices and facets
    // -> call compute_adjacent_facets_of_vertices() if not already done

    index_t nb_edges = mesh.edges.nb(),
            v0_index = index_t(-1),
            v1_index = index_t(-1),
            facet_corner_on_v0 = index_t(-1),
            facet_at_left = index_t(-1),
            facet_at_right = index_t(-1);
    GEO::vector<index_t> to_remove(nb_edges,0); // for each edge, store 0 = keep or 1 = remove
    CustomMeshHalfedges mesh_he(mesh);
    unsigned int nb_edges_removed = 0;
    CustomMeshHalfedges::Halfedge outgoing_halfedge;

    if (adj_facets.empty()) {
        compute_adjacent_facets_of_vertices(mesh,adj_facets);
    }

    FOR(e,nb_edges) {
        // get indices of the 2 vertices of the current edge
        v0_index = mesh.edges.vertex(e,0);
        v1_index = mesh.edges.vertex(e,1);
        // parse adjacent facets of v0 until we found an halfedge
        // whose origin is on v0 and tip on v1
        for(index_t f : adj_facets[v0_index]) {
            FOR(lv,mesh.facets.nb_vertices(f)) {
                facet_corner_on_v0 = mesh.facets.corner(f,lv);
                if(mesh.facet_corners.vertex(facet_corner_on_v0) == v0_index) {
                    break; // we found which facet corner is on v0
                }
            }
            // now that we have a facet and a facet cornet, we can have an halfedge going out of v0
            outgoing_halfedge.facet = f;
            outgoing_halfedge.corner = facet_corner_on_v0;
            geo_assert(halfedge_vertex_index_from(mesh,outgoing_halfedge) == v0_index)
            // check of the vertex at the tip is v1
            if(halfedge_vertex_index_to(mesh,outgoing_halfedge) == v1_index) {
                break;
            }
            // else : continue and check another adjacent facet
        }
        facet_at_left = halfedge_facet_left(mesh,outgoing_halfedge);
        facet_at_right = halfedge_facet_right(mesh,outgoing_halfedge);
        if(angle(mesh_facet_normal(mesh,facet_at_left),mesh_facet_normal(mesh,facet_at_right)) < FEATURE_EDGES_MIN_ANGLE) {
            to_remove[e] = 1;
            nb_edges_removed++;
        }
        // else : keep this feature edge because significant angle
    }
    
    mesh.edges.delete_elements(to_remove,false);

    fmt::println(
        Logger::out("feature edges"),
        "{} feature edge(s) removed ({:.1f}%) because dihedral angle < {}",
        nb_edges_removed,
        (nb_edges_removed * 100.0) / (double) nb_edges,
        FEATURE_EDGES_MIN_ANGLE);
    Logger::out("feature edges").flush();
}

void transfer_feature_edges(Mesh& mesh, std::set<std::pair<index_t,index_t>>& feature_edges) {
    feature_edges.clear();
    index_t v0 = index_t(-1),
            v1 = index_t(-1);
    FOR(e,mesh.edges.nb()) {
        v0 = mesh.edges.vertex(e,0);
        v1 = mesh.edges.vertex(e,1);
        feature_edges.insert(std::make_pair(
            std::min(v0,v1),
            std::max(v0,v1)
        ));
    }
    mesh.edges.clear();
}

bool halfedge_is_on_feature_edge(const Mesh& mesh, const MeshHalfedges::Halfedge& H, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    index_t v0 = halfedge_vertex_index_from(mesh,H),
            v1 = halfedge_vertex_index_to(mesh,H);
    return feature_edges.contains(std::make_pair(
        std::min(v0,v1),
        std::max(v0,v1)
    ));
}
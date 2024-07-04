#include <geogram/basic/memory.h>           // for GEO::vector
#include <geogram/points/principal_axes.h>  // for PrincipalAxes3d
#include <geogram/basic/vecg.h>             // for vec3, length()
#include <geogram/basic/matrix.h>           // for mat3
#include <geogram/mesh/mesh_geometry.h>     // for get_bbox()
#include <geogram/basic/numeric.h>          // for max_float64()

#include <fmt/core.h>
#include <fmt/ostream.h>    // to use fmt::print() on ostreams

#include <queue>
#include <vector>
#include <utility>          // for std::pair, std::make_pair()
#include <algorithm>        // for std::min(), std::max()
#include <initializer_list>
#include <set>
#include <cmath>            // for std::fabs()

#include "geometry.h"
#include "CustomMeshHalfedges.h"
#include "LabelingGraph.h"
#include "labeling.h" // for label2vector
#include "containers_macros.h" // for VECTOR_MAX_INDEX()

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

void compute_adjacent_facets_of_vertices(const Mesh& mesh, std::vector<std::vector<AdjacentFacetOfVertex>>& adj_facet_corners) {
    adj_facet_corners.clear();
    adj_facet_corners.resize(mesh.vertices.nb());

    // link each facet to adjacent vertices
    FOR(f,mesh.facets.nb()) {
        FOR(lv,mesh.facets.nb_vertices(f)) { // for each local vertices of the facet f
            adj_facet_corners[mesh.facets.vertex(f,lv)].push_back(AdjacentFacetOfVertex{f,lv});
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

bool local_edge_is_on_feature_edge(const Mesh& mesh, index_t facet_index, index_t local_edge, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    // local edge k is the one between local vertices k and (k+1)%3
    // see https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
    geo_assert(facet_index < mesh.facets.nb());
    geo_assert(local_edge < 3);
    index_t v0 = mesh.facets.vertex(facet_index,local_edge),
            v1 = mesh.facets.vertex(facet_index,(local_edge+1)%3);
    return feature_edges.contains(std::make_pair(
        std::min(v0,v1),
        std::max(v0,v1)
    ));
}

// return true if success
bool move_to_next_halfedge_on_feature_edge(const CustomMeshHalfedges& mesh_he, MeshHalfedges::Halfedge& H, const std::set<std::pair<index_t,index_t>>& feature_edges) {
    geo_assert(halfedge_is_on_feature_edge(mesh_he.mesh(),H,feature_edges));
    mesh_he.move_to_opposite(H);
    MeshHalfedges::Halfedge init_H = H;
    do {
        mesh_he.move_clockwise_around_vertex(H,true);
    } while(!halfedge_is_on_feature_edge(mesh_he.mesh(),H,feature_edges));
    // we found an halfedge on a feature edge !
    if(H == init_H) {
        // it's the same halfedge...
        mesh_he.move_to_opposite(H); // restore orientation
        return false;
    }
    // we found a *NEW* halfedge on a feature edge !
    return true;
}

void rotate_mesh_according_to_principal_axes(Mesh& mesh) {
    geo_assert(mesh.vertices.nb()!=0);

    ////////////////////////////////
    // Compute principal axes
    ////////////////////////////////
    
    PrincipalAxes3d principal_axes;
    principal_axes.begin();
    FOR(v,mesh.vertices.nb()) {
        principal_axes.add_point(mesh.vertices.point(v));
    }
    principal_axes.end();

    ////////////////////////////////
    // Compute vector toward (1,1,1) in the principal axes
    ////////////////////////////////

    vec3 rotated_vector = normalize(normalize(principal_axes.axis(0)) +
                                    normalize(principal_axes.axis(1)) + 
                                    normalize(principal_axes.axis(2)));
    
    ////////////////////////////////
    // Compute rotation matrix that align (1,1,1) (regarding init axes) to rotated_vector
    ////////////////////////////////

    vec3 reference(1.0,1.0,1.0);
    reference = normalize(reference);

    // thanks Jur van den Berg https://math.stackexchange.com/a/476311
    geo_assert(length(rotated_vector - vec3(-1.0,-1.0,-1.0)) > 10e-3); // formula not applicable if the 2 vectors are opposite
    vec3 v = cross(reference,rotated_vector); // = [v0 v1 v2]
    double c = dot(reference,rotated_vector); // <=> cosine of angle
    mat3 skew_symmetric_cross_product;
    /* 
     * = [[ 0  -v2  v1],
     *    [ v2  0  -v0],
     *    [-v1  v0  0 ]]
     */
    skew_symmetric_cross_product(0,0) = 0;
    skew_symmetric_cross_product(0,1) = -v[2];
    skew_symmetric_cross_product(0,2) = v[1];
    skew_symmetric_cross_product(1,0) = v[2];
    skew_symmetric_cross_product(1,1) = 0;
    skew_symmetric_cross_product(1,2) = -v[0];
    skew_symmetric_cross_product(2,0) = -v[1];
    skew_symmetric_cross_product(2,1) = v[0];
    skew_symmetric_cross_product(2,2) = 0;
    mat3 R; // initialized with identity
    R += skew_symmetric_cross_product + skew_symmetric_cross_product * skew_symmetric_cross_product * (1/(1+c));

    ////////////////////////////////
    // Rotate point cloud
    ////////////////////////////////

    FOR(v,mesh.vertices.nb()) {
        mesh.vertices.point(v) = mult(R,mesh.vertices.point(v));
    }

}

MeshHalfedges::Halfedge get_an_outgoing_halfedge_of_vertex(const Mesh& mesh, const std::vector<std::vector<index_t>>& adj_facets, index_t vertex_index) {
    MeshHalfedges::Halfedge output;
    output.facet = adj_facets[vertex_index][0];
    for(index_t facet_corner : mesh.facets.corners(output.facet)) { // for each facet corner of `output.facet`
        output.corner = facet_corner;
        if(halfedge_vertex_index_from(mesh,output) == vertex_index) {
            return output;
        }
    }
    geo_assert_not_reached;
}

MeshHalfedges::Halfedge get_halfedge_between_two_vertices(const CustomMeshHalfedges& mesh_he, const std::vector<std::vector<index_t>>& adj_facets, index_t origin_vertex, index_t extremity_vertex) {
    MeshHalfedges::Halfedge init_halfedge = get_an_outgoing_halfedge_of_vertex(mesh_he.mesh(),adj_facets,origin_vertex);
    MeshHalfedges::Halfedge current_halfedge = init_halfedge;
    index_t current_extremity_vertex = index_t(-1);

    do {
        current_extremity_vertex = halfedge_vertex_index_to(mesh_he.mesh(),current_halfedge);
        if(current_extremity_vertex == extremity_vertex) {
            return current_halfedge;
        }
        mesh_he.move_counterclockwise_around_vertex(current_halfedge,true); // ignore borders
    } while(current_halfedge != init_halfedge);
    geo_assert_not_reached;
}

MeshHalfedges::Halfedge get_most_aligned_halfedge_around_vertex(const MeshHalfedges::Halfedge& init_halfedge, const CustomMeshHalfedges& mesh_he, const vec3& reference) {
    vec3 normalized_reference = normalize(reference);
    MeshHalfedges::Halfedge current_halfedge = init_halfedge;
    MeshHalfedges::Halfedge most_aligned_halfedge; // init value is (NO_FACET,NO_CORNER)
    double current_cost = 0.0,
           min_cost = Numeric::max_float64();
    do {
        mesh_he.move_clockwise_around_vertex(current_halfedge,true);
        // dot product = 1  -> cost = 0
        // dot product = 0  -> cost = 0.5
        // dot product = -1 -> cost = 1
        // -> flip (*-1), halve range (*0.5) and shift (+0.5)
        current_cost = dot(normalize(halfedge_vector(mesh_he.mesh(),current_halfedge)),normalized_reference)*(-0.5)+0.5;
        if (current_cost < min_cost) {
            min_cost = current_cost;
            most_aligned_halfedge = current_halfedge;
        }
    } while (current_halfedge != init_halfedge);
    geo_assert(most_aligned_halfedge.facet != NO_FACET);
    geo_assert(most_aligned_halfedge.corner != NO_CORNER);
    return most_aligned_halfedge;
}

void get_adjacent_facets_conditional(const Mesh& mesh, index_t facet_index, index_t which_chart, const std::vector<index_t>& facet2chart, std::set<index_t>& out) {
    index_t neighbor_facet = index_t(-1);
    FOR(le,3) { // for each local edge of the facet
        neighbor_facet = mesh.facets.adjacent(facet_index,le);
        if(out.contains(neighbor_facet)) { // `neighbor_facet` found by another facet
            continue;
        }
        if(facet2chart.at(neighbor_facet) != which_chart) {
            continue;
        }
        out.insert(neighbor_facet);
    }
}

index_t get_nearest_vertex_of_coordinates(const CustomMeshHalfedges& mesh_he, const std::vector<std::vector<index_t>>& adj_facets, vec3 target_coordinates, index_t start_vertex, size_t max_dist) {
    const Mesh& mesh = mesh_he.mesh();
    geo_assert(mesh.vertices.double_precision());
    geo_assert(max_dist >= 1); // if max_dist==0, only `start_vertex` can be returned...

    index_t previous_vertex = index_t(-1);
    MeshHalfedges::Halfedge previous_halfedge;
    double previous_distance = nanf64("");
    index_t current_vertex = start_vertex;
    MeshHalfedges::Halfedge current_halfedge = get_an_outgoing_halfedge_of_vertex(mesh,adj_facets,start_vertex);
    double current_distance = distance(mesh.vertices.point(start_vertex),target_coordinates);
    size_t dist = 0;

    do {
        previous_halfedge = current_halfedge;
        previous_vertex = current_vertex;
        previous_distance = current_distance;
        dist++;
        current_halfedge = get_most_aligned_halfedge_around_vertex(previous_halfedge,mesh_he,target_coordinates-mesh.vertices.point(previous_vertex));
        mesh_he.move_to_opposite(current_halfedge);
        current_vertex = halfedge_vertex_index_from(mesh,current_halfedge);
        current_distance = distance(mesh.vertices.point(current_vertex),target_coordinates);
    } while ((previous_distance > current_distance) && (dist < max_dist));

    // now we are too far from the `target_coordinates`,
    // `previous_vertex` was closer than `current_vertex`
    // -> `previous_vertex` is the nearest vertex
    return previous_vertex;
}

index_t nearest_axis_of_edges(const Mesh& mesh, std::initializer_list<MeshHalfedges::Halfedge> edges, std::initializer_list<index_t> forbidden_axes) {
    std::set<index_t> forbidden_axes_as_set(forbidden_axes);
    geo_assert(forbidden_axes_as_set.size() < 3); // at least an axis must be allowed
    // find the axis to insert between the 2 part of this group
    // among the 2 others axes, find the one that is best suited for the current and next edge (= the second choice)
    std::array<double,3> per_axis_dot_product;
    per_axis_dot_product.fill(Numeric::min_float64());
    FOR(current_axis,3) {
        if(forbidden_axes_as_set.contains(current_axis)) {
            continue; // leave min_float64() value
        }
        // find the max dot product for this axis,
        //  for all edges in `edges` &
        //  for both the positive and negative axis direction (absolute dot product)
        for(const auto& edge : edges) {
            per_axis_dot_product[current_axis] = std::max(
                per_axis_dot_product[current_axis],
                std::fabs(dot(normalize(Geom::halfedge_vector(mesh,edge)),label2vector[current_axis*2]))
            );
        }
    }
    return (index_t) VECTOR_MAX_INDEX(per_axis_dot_product);
}

double dot_product_between_halfedge_and_axis(const Mesh& mesh, MeshHalfedges::Halfedge halfedge, index_t axis) {
    return std::fabs(dot(normalize(Geom::halfedge_vector(mesh,halfedge)),label2vector[axis*2]));
}

double angle_between_halfedge_and_axis(const Mesh& mesh, MeshHalfedges::Halfedge halfedge, index_t axis) {
    return std::min(
        angle(normalize(Geom::halfedge_vector(mesh,halfedge)),label2vector[axis*2]), // angle between `halfedge` and the positive direction of `axis` (eg +X)
        angle(normalize(Geom::halfedge_vector(mesh,halfedge)),label2vector[axis*2+1]) // angle between `halfedge` and the negative direction of `axis` (eg -X)
    );
}

void trace_path_on_chart(const CustomMeshHalfedges& mesh_he, const std::vector<std::vector<index_t>>& adj_facets, const std::vector<index_t>& facet2chart, const std::map<index_t,std::vector<std::pair<index_t,index_t>>>& turning_point_vertices, index_t start_vertex, vec3 direction, std::set<index_t>& facets_at_left, std::set<index_t>& facets_at_right, std::vector<MeshHalfedges::Halfedge>& halfedges) {
    // From `start_vertex`, move edge by edge in the given `direction`
    // To not drift away from this direction, we aim a target point, and not the same vector each time
    // To allow to go further than the init target point, we move the target point away at each step
    const Mesh& mesh = mesh_he.mesh();
    vec3 target_point = mesh_vertex(mesh,start_vertex) + direction;
    vec3 normalized_direction = normalize(direction);
    MeshHalfedges::Halfedge current_halfedge = get_an_outgoing_halfedge_of_vertex(mesh, adj_facets, start_vertex),
                            next_halfedge;
    next_halfedge = get_most_aligned_halfedge_around_vertex(current_halfedge,mesh_he,direction);
    // get which chart is around this halfedge
    index_t supporting_chart = index_t(-1); // will be computed at the 2nd halfedge
    // start walking
    do {
        current_halfedge = next_halfedge;
        facets_at_left.insert(halfedge_facet_left(mesh,current_halfedge));
        facets_at_right.insert(halfedge_facet_right(mesh,current_halfedge));
        halfedges.push_back(current_halfedge);
        // Stop if we found a turning point
        // see MAMBO B76 for example
        // https://gitlab.com/franck.ledoux/mambo
        if(turning_point_vertices.contains(halfedge_vertex_index_to(mesh,current_halfedge))) {
            break;
        }
        mesh_he.move_to_opposite(current_halfedge); // flip `current_halfedge` so that its origin vertex is the origin vertex of the next halfedge
        next_halfedge = get_most_aligned_halfedge_around_vertex(current_halfedge,mesh_he,target_point - halfedge_vertex_from(mesh,current_halfedge)); // reference vector = vector toward the `target_point`
        geo_assert(next_halfedge != current_halfedge); // assert the previous halfedge is not the one the best aligned. else we are backtracking
        target_point += normalized_direction * length(halfedge_vector(mesh,current_halfedge)); // move the target forward
        if(supporting_chart == index_t(-1)) {
            supporting_chart = facet2chart[Geom::halfedge_facet_left(mesh,next_halfedge)];
            geo_assert(facet2chart[Geom::halfedge_facet_right(mesh,next_halfedge)] == supporting_chart);
        }
    } while(
        facet2chart[Geom::halfedge_facet_left(mesh,next_halfedge)] == supporting_chart &&
        facet2chart[Geom::halfedge_facet_right(mesh,next_halfedge)] == supporting_chart
    );
    geo_assert(!facets_at_left.empty());
    geo_assert(!facets_at_right.empty());
}

vec3 average_facets_normal(const std::vector<vec3>& normals, const std::set<index_t>& facets) {
    geo_assert(!normals.empty());
    vec3 sum(0.0,0.0,0.0);
    for(const auto& f : facets) {
        sum += normals[f];
    }
    return sum / (double) facets.size();
}

double average_dot_product(const std::vector<vec3>& normals, const std::set<index_t>& facets, vec3 reference) {
    geo_assert(!normals.empty());
    vec3 normalized_reference = normalize(reference);
    double sum = 0.0;
    for(const auto& f : facets) {
        sum += dot(normals[f],normalized_reference);
    }
    return sum / (double) facets.size();
}

double average_angle(const std::vector<vec3>& normals, const std::set<index_t>& facets, vec3 reference) {
    geo_assert(!normals.empty());
    vec3 normalized_reference = normalize(reference);
    double sum = 0.0;
    for(const auto& f : facets) {
        sum += angle(normals[f],normalized_reference);
    }
    return sum / (double) facets.size();
}

bool vertex_has_lost_feature_edge_in_neighborhood(const CustomMeshHalfedges& mesh_he, const std::vector<std::vector<index_t>>& adj_facets, const std::set<std::pair<index_t,index_t>>& feature_edges, index_t vertex, MeshHalfedges::Halfedge& outgoing_halfedge_on_feature_edges) {
    geo_assert(mesh_he.is_using_facet_region()); // we need the labeling to be stored as facet region in `mesh_he`, so that halfedges on borders are defined
    if(feature_edges.empty()) {
        return false;
    }
    const Mesh& mesh = mesh_he.mesh();
    MeshHalfedges::Halfedge init_halfedge = get_an_outgoing_halfedge_of_vertex(mesh,adj_facets,vertex);
    MeshHalfedges::Halfedge current_halfedge = init_halfedge;
    do { // explore all outgoing halfedges of this vertex, starting with `init_halfedge`
        if(halfedge_is_on_feature_edge(mesh,current_halfedge,feature_edges) && !mesh_he.halfedge_is_border(current_halfedge)) {
            // `current_halfedge` is on a "lost" feature edge
            outgoing_halfedge_on_feature_edges = current_halfedge;
            return true;
        }
        mesh_he.move_clockwise_around_vertex(current_halfedge,true); // go to next halfedge clockwise & ignore borders
    } while(current_halfedge != init_halfedge);
    return false; // found no lost feature edge
}

MeshHalfedges::Halfedge follow_feature_edge_on_chart(const CustomMeshHalfedges& mesh_he, MeshHalfedges::Halfedge halfedge, const std::set<std::pair<index_t,index_t>>& feature_edges, const std::vector<index_t>& facet2chart, std::set<index_t>& facets_at_left, std::set<index_t>& facets_at_right) {
    const Mesh& mesh = mesh_he.mesh();
    geo_assert(halfedge_is_on_feature_edge(mesh,halfedge,feature_edges));
    MeshHalfedges::Halfedge previous_halfedge;
    index_t chart_on_which_the_lost_feature_edge_is = facet2chart[halfedge_facet_left(mesh,halfedge)];
    geo_assert(chart_on_which_the_lost_feature_edge_is == facet2chart[halfedge_facet_right(mesh,halfedge)]);
    bool move_successful = false;
    // walk along feature edge until we are no longer on the chart
    do {
        facets_at_left.insert(halfedge_facet_left(mesh,halfedge));
        facets_at_right.insert(halfedge_facet_right(mesh,halfedge));
        previous_halfedge = halfedge;
        move_successful = move_to_next_halfedge_on_feature_edge(mesh_he,halfedge,feature_edges);
        geo_assert(halfedge_is_on_feature_edge(mesh,halfedge,feature_edges));
    } while(
        move_successful && 
        (facet2chart[halfedge_facet_left(mesh,halfedge)] == chart_on_which_the_lost_feature_edge_is) &&
        (facet2chart[halfedge_facet_right(mesh,halfedge)] == chart_on_which_the_lost_feature_edge_is)
    );
    return previous_halfedge;
}

bool is_a_facet_to_tilt(const vec3& facet_normal, double sensitivity) {
    std::array<std::pair<double,index_t>,6> per_label_weights = {
        std::make_pair(facet_normal.x < 0.0 ? 0.0 : facet_normal.x,     0), // +X
        std::make_pair(facet_normal.x < 0.0 ? -facet_normal.x : 0.0,    1), // -X
        std::make_pair(facet_normal.y < 0.0 ? 0.0 : facet_normal.y,     2), // +Y
        std::make_pair(facet_normal.y < 0.0 ? -facet_normal.y : 0.0,    3), // -Y
        std::make_pair(facet_normal.z < 0.0 ? 0.0 : facet_normal.z,     4), // +Z
        std::make_pair(facet_normal.z < 0.0 ? -facet_normal.z : 0.0,    5)  // -Z
    };
    std::sort(per_label_weights.begin(),per_label_weights.end());
    return std::fabs(per_label_weights[5].first-per_label_weights[4].first) < sensitivity;
}

size_t get_facets_to_tilt(const Mesh& mesh, const std::vector<vec3>& normals, std::set<index_t>& facets_to_tilt, double sensitivity) {
    facets_to_tilt.clear();
    FOR(f,mesh.facets.nb()) { // for each facet
        if(is_a_facet_to_tilt(normals[f],sensitivity)) {
            // the 2 labels with the most weight are too close
            facets_to_tilt.insert(f);
        }
    }
    return facets_to_tilt.size();
}

double sd_adjacent_facets_area(const Mesh& mesh, const std::vector<std::vector<index_t>>& adj_facets, index_t vertex_index) {
    BasicStats stats;
    for(index_t f : adj_facets[vertex_index]) {
        stats.insert(mesh_facet_area(mesh,f));
    }
    return stats.sd();
}

void triangulate_facets(Mesh& M, std::vector<index_t>& triangle_index_to_old_facet_index, std::vector<index_t>& corner_index_to_old_corner_index) {
    // based on ext/geogram/src/lib/geogram/mesh/mesh.cpp MeshFacets::triangulate()
    // but returns the map between triangle indices and facet indices before the triangulation
    // allowing to transfer facet attributes
    // + also returns the map between facet corners before/after the triangulation
    if(M.facets.are_simplices()) {
        return;
    }
    triangle_index_to_old_facet_index.clear();
    corner_index_to_old_corner_index.clear();
    index_t nb_triangles = 0;
    for(index_t f = 0; f < M.facets.nb(); f++) {
        FOR(i,M.facets.nb_vertices(f) - 2) {
            triangle_index_to_old_facet_index.push_back(f);
            nb_triangles++;
        }
    }
    vector<index_t> new_corner_vertex_index;
    new_corner_vertex_index.reserve(nb_triangles * 3);
    corner_index_to_old_corner_index.reserve(nb_triangles * 3);
    for(index_t f = 0; f < M.facets.nb(); f++) {
        index_t v0 = M.facet_corners.vertex(M.facets.corners_begin(f));
        for(index_t c = M.facets.corners_begin(f) + 1;
            c + 1 < M.facets.corners_end(f); ++c
        ) {
            new_corner_vertex_index.push_back(v0);
            corner_index_to_old_corner_index.push_back(M.facets.corners_begin(f));
            new_corner_vertex_index.push_back(
                M.facet_corners.vertex(c)
            );
            corner_index_to_old_corner_index.push_back(c);
            new_corner_vertex_index.push_back(
                M.facet_corners.vertex(c + 1)
            );
            corner_index_to_old_corner_index.push_back(c+1);
        }
    }
    M.facets.assign_triangle_mesh(new_corner_vertex_index, true);
}

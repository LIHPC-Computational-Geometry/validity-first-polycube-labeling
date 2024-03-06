#include <geogram/mesh/mesh.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/assert.h>
#include <geogram/mesh/mesh_halfedges.h>    // for MeshHalfedges::Halfedge

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/os.h>

#include <DisjointSet.hpp>

#include <GCoptimization.h>

#include <utility>  // for std::pair
#include <fstream>
#include <vector>
#include <tuple>    // for std::tie()
#include <numbers>  // for std::numbers::pi

#include <nlohmann/json.hpp>

#include "LabelingGraph.h"
#include "geometry.h"               // for other_axis(), HalfedgeCompare
#include "containers.h"             // for VECTOR_CONTAINS(), MAP_CONTAINS(), index_of_last(), += on std::vector
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges
#include "labeling.h"               // for label2vector
#include "dump_mesh.h"              // for dump_vertex()

bool Chart::is_surrounded_by_feature_edges(const std::vector<Boundary>& all_boundaries) const {
    geo_assert(!boundaries.empty());
    for(index_t boundary_index : boundaries) { // for each boundary around this chart
        if(all_boundaries[boundary_index].on_feature_edge == false) {
            return false; // no, at least a subset of the contour is not on feature edges
        }
    }
    return true;
}

void Chart::counterclockwise_boundaries_order(
    const CustomMeshHalfedges& mesh_he,
    const std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare> halfedge2boundary,
    const std::vector<Boundary>& all_boundaries,
    std::vector<std::pair<index_t,bool>>& counterclockwise_order) const {

    // 1. get the fist halfedge of the fist element in the `boundaries` set
    // 2. if the current chart is at the right, we are not going counterclockwise -> register the boundary as in the opposite direction
    //    if the current chart is at the left -> register the boundary as the the same direction + go to the opposite corner
    // 3. Move clockwise until we are on a boundary halfedge (with the current chart on the left)
    // 4. Get the boundary associated to this halfedge, register its index and if it's in the same direction or not
    // ...
    // stop when we ordered all boundaries
    
    geo_assert(mesh_he.is_using_facet_region());
    geo_assert(!boundaries.empty());
    counterclockwise_order.clear();
    const Mesh& mesh = mesh_he.mesh();
    bool same_direction = false;

    index_t index_of_current_boundary = *boundaries.begin(); // first boundary in the set
    MeshHalfedges::Halfedge current_halfedge = all_boundaries[index_of_current_boundary].halfedges[0]; // first halfedge of first boundary
    if(facets.contains(halfedge_facet_left(mesh,current_halfedge))) {
        counterclockwise_order.push_back(std::make_pair(index_of_current_boundary,true)); // counterclockwise, the boundary is in the same direction
        current_halfedge = *all_boundaries[index_of_current_boundary].halfedges.rbegin(); // get last halfedge
        mesh_he.move_to_opposite(current_halfedge);
    }
    else if(facets.contains(halfedge_facet_right(mesh,current_halfedge))) {
        counterclockwise_order.push_back(std::make_pair(index_of_current_boundary,false)); // counterclockwise, the boundary is in the opposite direction
    }
    else {
        geo_assert_not_reached;
    }
    // now the base vertex of `current_halfedge` is on the next corner to explore

    unsigned int nb_iter = 0;
    do {
        MeshHalfedges::Halfedge previous_halfedge = current_halfedge;
        mesh_he.move_clockwise_around_vertex_until_on_border(current_halfedge);
        geo_assert(current_halfedge != previous_halfedge);
        geo_assert(facets.contains(halfedge_facet_left(mesh,current_halfedge)));
        std::tie(index_of_current_boundary, same_direction) = halfedge2boundary.at(current_halfedge);
        geo_assert(boundaries.contains(index_of_current_boundary));
        counterclockwise_order.push_back(std::make_pair(index_of_current_boundary,same_direction));
        if(same_direction) {
            current_halfedge = *all_boundaries[index_of_current_boundary].halfedges.rbegin(); // get last halfedge
            mesh_he.move_to_opposite(current_halfedge);
        }
        else {
            current_halfedge = *all_boundaries[index_of_current_boundary].halfedges.begin(); // get first halfedge
        }

        nb_iter++;
        if(nb_iter > 100) {
            fmt::println(Logger::err("fix_labeling"),"Infinite loop in Chart::counterclockwise_boundaries_order(), reached max iter"); Logger::err("fix_labeling").flush();
            break;
        }
    } while (counterclockwise_order.size() != boundaries.size());
}

std::ostream& operator<< (std::ostream &out, const Chart& data) {
    fmt::println(out,"\tlabel : {}",data.label);
    fmt::println(out,"\tfacets : {}",data.facets);
    fmt::println(out,"\tboundaries : {}",data.boundaries);
    return out;
}

std::size_t VertexRingWithBoundaries::valence() const {
    return boundary_edges.size();
}

bool VertexRingWithBoundaries::halfedge_is_in_boundary_edges(const MeshHalfedges::Halfedge& halfedge) const {
    if(VECTOR_CONTAINS(boundary_edges,halfedge)) {
        return true;
    }
    return false;
}

/**
 * \brief Explore the vertex at the base of \c initial_halfedge and fill the struct with boundary edges encountered
 * \param[in] initial_halfedge An half-edge going outward the vertex to explore
 * \param[in] mesh_halfedges The half-edges interface of the mesh
 * \param[in] labeling The facet-to-label association, to detect boundary edges (= edges between 2 labels)
 */
void VertexRingWithBoundaries::explore(const MeshHalfedges::Halfedge& initial_halfedge,
                                       const CustomMeshHalfedges& mesh_halfedges) {
    
    // prepare the vertex exploration
    MeshHalfedges::Halfedge current_halfedge = initial_halfedge; // create another Halfedge that we can modify (we need to keep the initial one)

    // go around the vertex, from boundary edge to boundary edge
    do {
        mesh_halfedges.move_counterclockwise_around_vertex(current_halfedge,true);
        if(mesh_halfedges.halfedge_is_border(current_halfedge)) {
            boundary_edges.push_back(current_halfedge); // register the current halfedge as outgoing boundary edge for this vertex ring
        }
    } while ((current_halfedge != initial_halfedge));
}

void VertexRingWithBoundaries::explore(index_t init_facet_corner, const CustomMeshHalfedges& mesh_halfedges) {
    MeshHalfedges::Halfedge initial_halfedge;
    initial_halfedge.facet = mesh_halfedges.mesh().facet_corners.adjacent_facet(init_facet_corner);
    index_t init_vertex = mesh_halfedges.mesh().facet_corners.vertex(init_facet_corner);
    // find which facet corner is init_vertex from initial_halfedge.facet
    FOR(lv,3) {
        if(mesh_halfedges.mesh().facet_corners.vertex(
            mesh_halfedges.mesh().facets.corner(initial_halfedge.facet,lv)
        ) == init_vertex) {
            initial_halfedge.corner = mesh_halfedges.mesh().facets.corner(initial_halfedge.facet,lv);
            break;
        }
    }
    geo_assert(mesh_halfedges.halfedge_is_valid(initial_halfedge));
    geo_assert(Geom::halfedge_vertex_index_from(mesh_halfedges.mesh(),initial_halfedge)==init_vertex);
    explore(initial_halfedge,mesh_halfedges);
}

void VertexRingWithBoundaries::check_boundary_edges(const CustomMeshHalfedges& mesh_halfedges) const {
    for(const auto& be: boundary_edges) {
        geo_assert(mesh_halfedges.halfedge_is_border(be));
    }
}

std::size_t VertexRingWithBoundaries::circular_previous(std::size_t index) const {
    return (index == 0) ? boundary_edges.size()-1 : --index;
}

std::size_t VertexRingWithBoundaries::circular_next(std::size_t index) const {
    return (index == boundary_edges.size()-1) ? 0 : ++index;
}

std::ostream& operator<< (std::ostream &out, const VertexRingWithBoundaries& data) {
    fmt::print(out,"{}",data.boundary_edges);
    return out;
}

std::size_t Corner::valence() const {
    // sum the valence of all vertex rings
    std::size_t valence = 0;
    for (const VertexRingWithBoundaries& vertex_ring_with_boundaries: vertex_rings_with_boundaries) {
        valence += vertex_ring_with_boundaries.valence();
    }
    return valence;
}

bool Corner::halfedge_is_in_boundary_edges(const MeshHalfedges::Halfedge& halfedge) const {
    for(const VertexRingWithBoundaries& vertex_ring_with_boundaries: vertex_rings_with_boundaries) { // parse boundary edges grouped by vertex ring
        if(vertex_ring_with_boundaries.halfedge_is_in_boundary_edges(halfedge)) {
            return true;
        }
    }
    return false;
}

bool Corner::compute_validity(bool allow_boundaries_between_opposite_labels, const std::vector<Boundary>& boundaries, const std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary) {
    int per_axis_counter[3] = {0,0,0}; // number of X, Y, and Z boundaries
    std::size_t nb_boundary_edges_in_vertex_ring;
    // going around all boundary edges of all vertex rings and count boundary axes
    for(const auto& vr : vertex_rings_with_boundaries) {
        nb_boundary_edges_in_vertex_ring = vr.boundary_edges.size();
        std::vector<int> boundary_edges_axis(nb_boundary_edges_in_vertex_ring);
        std::vector<bool> boundary_edges_masked(nb_boundary_edges_in_vertex_ring,false);
        // mask boundary edges around axis-less boundary edges
        FOR(lbe,nb_boundary_edges_in_vertex_ring) { // lbe = local boundary edge index
            boundary_edges_axis[lbe] = boundaries[halfedge2boundary.at(vr.boundary_edges[lbe]).first].axis; // get the boundary index, then the axis of this boundary
            if(boundary_edges_axis[lbe]==-1) {
                if(allow_boundaries_between_opposite_labels==false) {
                    is_valid=false;
                    return false;
                }
                // do not count the neighboring boundaries
                boundary_edges_masked[vr.circular_previous(lbe)] = true;
                boundary_edges_masked[vr.circular_next(lbe)] = true;
            }
        }
        // count axes
        FOR(lbe,nb_boundary_edges_in_vertex_ring) {
            if( (boundary_edges_masked[lbe]==false) && (boundary_edges_axis[lbe]!=-1) ) {
                per_axis_counter[boundary_edges_axis[lbe]]++;
            }
        }
    }

    // if after masking, the "corner" is now in the volume -> valid
    if ((per_axis_counter[0]==0) && (per_axis_counter[1]==0) && (per_axis_counter[2]==0)) {
        is_valid = true;
        return true;
    }

    // if after masking, the "corner" is just lying on a boundary -> valid
    if(
    ((per_axis_counter[0]==2) && (per_axis_counter[1]==0) && (per_axis_counter[2]==0)) ||  // only 2 X
    ((per_axis_counter[0]==0) && (per_axis_counter[1]==2) && (per_axis_counter[2]==0)) ||  // only 2 Y
    ((per_axis_counter[0]==0) && (per_axis_counter[1]==0) && (per_axis_counter[2]==2)) ) { // only 2 Z
        is_valid = true;
        return true;
    }

    // if it remains 1, 3 or more boundaries of the same axis -> invalid
    if(
    ((per_axis_counter[1]==0) && (per_axis_counter[2]==0)) || 
    ((per_axis_counter[0]==0) && (per_axis_counter[2]==0)) || 
    ((per_axis_counter[0]==0) && (per_axis_counter[1]==0)) ) {
        is_valid = false;
        return false;
    }
    
    // to be valid, boundary axes should match two by two, leaving nothing or an {X,Y,Z} group
    per_axis_counter[0] %= 2;
    per_axis_counter[1] %= 2;
    per_axis_counter[2] %= 2;
    is_valid = ((per_axis_counter[0]==0) && (per_axis_counter[1]==0) && (per_axis_counter[2]==0)) ||
               ((per_axis_counter[0]==1) && (per_axis_counter[1]==1) && (per_axis_counter[2]==1));
    return is_valid;
}

bool Corner::all_adjacent_boundary_edges_are_on_feature_edges(const Mesh& mesh, const std::set<std::pair<index_t,index_t>>& feature_edges) const {
    for(const auto& vr : vertex_rings_with_boundaries) {
        for(const auto& boundary_edge : vr.boundary_edges) {
            if(!halfedge_is_on_feature_edge(mesh,boundary_edge,feature_edges)) {
                return false;
            }
        }
    }
    return true;
}

vec3 Corner::average_coordinates_of_neighborhood(const Mesh& mesh, const StaticLabelingGraph& slg, bool include_itself, size_t max_dist) const {
    geo_assert(vertex != index_t(-1));
    geo_assert(max_dist >= 1);
    geo_assert(mesh.vertices.double_precision());

    vec3 sum(0.0,0.0,0.0);
    size_t count_vertices_in_sum = 0;
    index_t boundary_index = index_t(-1);
    bool same_direction = false;
    size_t dist = 0;

    if(include_itself) {
        // equivalent to a distance of 0
        sum += mesh.vertices.point(vertex);
        count_vertices_in_sum++;
    }

    for(const auto& vr : vertex_rings_with_boundaries) {
        for(auto first_boundary_halfedge : vr.boundary_edges) {
            std::tie(boundary_index,same_direction) = slg.halfedge2boundary.at(first_boundary_halfedge);
            const Boundary& boundary = slg.boundaries.at(boundary_index);
            if(same_direction) { // the corner is the start_corner of `boundary`
                dist = 1;
                for(auto boundary_halfedge : boundary.halfedges) { // go accross the boundary
                    sum += halfedge_vertex_to(mesh,boundary_halfedge) /* / (double) dist */; // the further the vertex is from the corner, the lesser its contribution to the sum
                    count_vertices_in_sum++;
                    dist++;
                    if(dist > max_dist) {
                        break;
                    }
                }
            }
            else { // the corner is the end_corner of `boundary`
                dist = 1;
                for(auto boundary_halfedge = boundary.halfedges.rbegin(); boundary_halfedge != boundary.halfedges.rend(); boundary_halfedge++) { // go accross the boundary in reverse order of halfedges
                    sum += halfedge_vertex_from(mesh,*boundary_halfedge) /* / (double) dist */; // the further the vertex is from the corner, the lesser its contribution to the sum
                    count_vertices_in_sum++;
                    dist++;
                    if(dist > max_dist) {
                        break;
                    }
                }
            }
        }
    }
    return sum / (double) count_vertices_in_sum;
}

MeshHalfedges::Halfedge Corner::get_most_aligned_boundary_halfedge(const Mesh& mesh, const vec3& reference) const {
    // based on geometry.cpp / get_most_aligned_halfedge_around_vertex()

    vec3 normalized_reference = normalize(reference);
    MeshHalfedges::Halfedge most_aligned_halfedge; // init value is (NO_FACET,NO_CORNER)
    double current_dot_product = 0.0,
           max_dot_product = Numeric::min_float64();

    for(const auto& vr : vertex_rings_with_boundaries) {
        for(const auto& boundary_halfedge : vr.boundary_edges) {
            current_dot_product = dot(normalize(halfedge_vector(mesh,boundary_halfedge)),normalized_reference);
            if (current_dot_product > max_dot_product) {
                max_dot_product = current_dot_product;
                most_aligned_halfedge = boundary_halfedge;
            }
        }
    }

    geo_assert(most_aligned_halfedge.facet != NO_FACET);
    geo_assert(most_aligned_halfedge.corner != NO_CORNER);
    return most_aligned_halfedge;    
}

std::ostream& operator<< (std::ostream &out, const Corner& data) {
    fmt::println(out,"\tvertex : {}",data.vertex);
    fmt::println(out,"\tis_valid : {}",data.is_valid);
    fmt::println(out,"\t{} vertex ring(s) with boundaries : {}",data.vertex_rings_with_boundaries.size(),data.vertex_rings_with_boundaries);
    return out;
}

void TurningPoint::fill_from(index_t outgoing_local_halfedge_index, const Boundary& boundary, const CustomMeshHalfedges& mesh_he) {
    geo_assert(mesh_he.is_using_facet_region()); // assert a facet region is defined in mesh_he
    outgoing_local_halfedge_index_ = outgoing_local_halfedge_index;
    // Explore clockwise : right side then left side
    // For each side, sum of angles to know towards with side the turning point is
    double right_side_sum_of_angles = 0.0;
    double left_side_sum_of_angles = 0.0;
    MeshHalfedges::Halfedge current_halfedge = boundary.halfedges[outgoing_local_halfedge_index], previous_halfedge = current_halfedge;
    mesh_he.move_clockwise_around_vertex(previous_halfedge,true);
    // explore the right side
    do {
        right_side_sum_of_angles += angle(halfedge_vector(mesh_he.mesh(),previous_halfedge),halfedge_vector(mesh_he.mesh(),current_halfedge));
        previous_halfedge = current_halfedge;
        mesh_he.move_clockwise_around_vertex(current_halfedge,true);
    } while (!mesh_he.halfedge_is_border(current_halfedge));
    // explore the left side
    do {
        left_side_sum_of_angles += angle(halfedge_vector(mesh_he.mesh(),previous_halfedge),halfedge_vector(mesh_he.mesh(),current_halfedge));
        previous_halfedge = current_halfedge;
        mesh_he.move_clockwise_around_vertex(current_halfedge,true);
    } while (!mesh_he.halfedge_is_border(current_halfedge));
    is_towards_right_ = (right_side_sum_of_angles > left_side_sum_of_angles);
}

index_t TurningPoint::get_closest_corner(const Boundary& boundary, const CustomMeshHalfedges& mesh_he) const {
    // the turning-point is at the base of the halfedge 'outgoing_local_halfedge_index_'
    // compute 2 distances :
    // - from boundary.start_corner to turning point (= first halfedge to outgoing_local_halfedge_index_-1 included)
    // - from turning point to boundary.end_corner (= outgoing_local_halfedge_index_ to last halfedge included)
    double distance_to_start_corner = 0.0;
    double distance_to_end_corner = 0.0;
    FOR(i,outgoing_local_halfedge_index_-1) {
        distance_to_start_corner += Geom::edge_length(mesh_he.mesh(),boundary.halfedges[i]);
    }
    for(index_t i = outgoing_local_halfedge_index_; i < boundary.halfedges.size(); i++) {
        distance_to_end_corner += Geom::edge_length(mesh_he.mesh(),boundary.halfedges[i]);
    }
    return (distance_to_start_corner < distance_to_end_corner) ? boundary.start_corner : boundary.end_corner;
}

index_t TurningPoint::vertex(const Boundary& boundary, const Mesh& mesh) const {
    return halfedge_vertex_index_from(mesh,boundary.halfedges[outgoing_local_halfedge_index_]);
}

std::ostream& operator<< (std::ostream &out, const TurningPoint& data) {
    fmt::println(out,"\t\toutgoing_local_halfedge_index : {}, direction : {}",data.outgoing_local_halfedge_index_,data.is_towards_left() ? "left" : "right");
    return out;
}

bool Boundary::empty() const {
    return (
        axis == -1 &&
        halfedges.empty() &&
        left_chart == index_t(-1) &&
        right_chart == index_t(-1) &&
        start_corner == index_t(-1) &&
        end_corner == index_t(-1) &&
        on_feature_edge == true
    );
}

void Boundary::explore(const MeshHalfedges::Halfedge& initial_halfedge,
                       const CustomMeshHalfedges& mesh_halfedges,
                       index_t index_of_self,
                       const std::set<std::pair<index_t,index_t>>& feature_edges,
                       const std::vector<index_t>& facet2chart,
                       std::vector<index_t>& vertex2corner,
                       std::vector<Chart>& charts,
                       std::vector<Corner>& corners,
                       std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary,
                       std::vector<MeshHalfedges::Halfedge>& boundary_edges_to_explore) {

    // only the start_corner should be filled
    geo_assert(axis == -1);
    geo_assert(turning_points.empty());
    geo_assert(halfedges.empty());
    geo_assert(left_chart == index_t(-1));
    geo_assert(right_chart == index_t(-1));
    geo_assert(end_corner == index_t(-1));
    geo_assert(on_feature_edge == true);

    // initialization
    index_t current_vertex = index_t(-1);
    const Mesh& mesh = mesh_halfedges.mesh();
    MeshHalfedges::Halfedge current_halfedge = initial_halfedge; // create a modifiable halfedge
    if(feature_edges.empty()) { // if this mesh doesn't have feature edges
        on_feature_edge = false; // this boundary is not on a feature edge
    }

    geo_assert(mesh_halfedges.halfedge_is_border(current_halfedge));
    halfedges.push_back(current_halfedge);
    halfedge2boundary[current_halfedge] = {index_of_self,true}; // mark this halfedge, link it to the current boundary
    mesh_halfedges.move_to_opposite(current_halfedge); // switch orientation
    halfedge2boundary[current_halfedge] = {index_of_self,false}; // mark this halfedge, link it to the current boundary (opposite orientation)
    mesh_halfedges.move_to_opposite(current_halfedge); // switch orientation
    if(!halfedge_is_on_feature_edge(mesh_halfedges.mesh(),current_halfedge,feature_edges)) {
        on_feature_edge = false;
    }

    // Step 1 : store information we already have with the first boundary edge:
    // the charts at the left and right, and the axis
    
    // fill left_chart field
    left_chart = facet2chart[Geom::halfedge_facet_left(mesh,initial_halfedge)]; // link the current boundary to the chart at its left
    charts[left_chart].boundaries.emplace(index_of_self); // link the left chart to the current boundary

    // fill right_chart field
    right_chart = facet2chart[Geom::halfedge_facet_right(mesh,initial_halfedge)]; // link the current boundary to the chart at its right
    charts[right_chart].boundaries.emplace(index_of_self); // link the left chart to the current boundary

    // fill axis field
    axis = other_axis(
        label2axis(charts[left_chart].label),
        label2axis(charts[right_chart].label)
    );

    // start computation of the average normal
    average_normal = halfedge_normal(mesh,current_halfedge);

    // Step 2 : go from boundary edge to boundary edge until we found a corner
    
    do {
        current_vertex = Geom::halfedge_vertex_index_to(mesh,current_halfedge);

        // compute the valence of the extremity_vertex
        //mesh_halfedges.move_to_opposite(current_halfedge); // flip current_halfedge so than extremity_vertex is at its base
        VertexRingWithBoundaries current_vertex_ring;
        current_vertex_ring.explore(Geom::halfedge_top_right_corner(mesh,current_halfedge),mesh_halfedges); // explore the halfedges around the extremity vertex

        if(current_vertex_ring.valence() < 3) {
            // Not a corner. At least, from the point of view of this vertex ring

            mesh_halfedges.move_to_next_around_border(current_halfedge); // cross current_vertex

            if(VECTOR_CONTAINS(halfedges,current_halfedge)) {
                // we went back on our steps
                // -> case of a boundary with no corner
                geo_assert(start_corner==index_t(-1));
                break;
            }

            halfedge2boundary[current_halfedge] = {index_of_self,true}; // mark this halfedge, link it to the current boundary
            halfedges.push_back(current_halfedge); // append to the list of halfedges composing the current boundary

            mesh_halfedges.move_to_opposite(current_halfedge); // switch orientation
            halfedge2boundary[current_halfedge] = {index_of_self,false}; // mark this halfedge, link it to the current boundary (opposite orientation)
            mesh_halfedges.move_to_opposite(current_halfedge); // switch orientation

            // update the averge normal of the boundary
            average_normal += halfedge_normal(mesh,current_halfedge);

            // check if we (still) are on a feature edge
            if(!halfedge_is_on_feature_edge(mesh_halfedges.mesh(),current_halfedge,feature_edges)) {
                on_feature_edge = false;
            }

            // the vertex at the base of halfedge will be explored in the next interation
        }
        else {
            // else : we found a corner !
            end_corner = vertex2corner[current_vertex];
            if(end_corner != index_t(-1)) {
                // this vertex has already been explored, maybe not by this vertex ring
                if(!corners[end_corner].halfedge_is_in_boundary_edges(current_vertex_ring.boundary_edges.front())) {
                    // this vertex has already been explored, but by another vertex ring
                    boundary_edges_to_explore += current_vertex_ring.boundary_edges;
                    corners[end_corner].vertex_rings_with_boundaries.push_back(current_vertex_ring);
                }
                // else : this vertex has already been explored by the same vertex ring
            }
            else {
                // we found an unexplored corner
                corners.push_back(Corner()); // create a new Corner
                corners.back().vertex = current_vertex; // link it to the current vertex
                corners.back().vertex_rings_with_boundaries.push_back(current_vertex_ring); // add it the just-explored vertex ring
                end_corner = (index_t) index_of_last(corners); // link the boundary to this corner
                vertex2corner[current_vertex] = end_corner; // link the current vertex to this corner
                boundary_edges_to_explore += current_vertex_ring.boundary_edges; // the boundary edges around this corner must be explored
            }
        }
    } while (end_corner == index_t(-1)); // continue the exploration until a corner is found
}

bool Boundary::contains_lower_than_180_degrees_angles(const CustomMeshHalfedges& mesh_halfedges) {
    for(auto be : halfedges) { // need to be mutable, but will be back on the original halfedge
        if(mesh_halfedges.is_on_lower_than_180_degrees_edge(be)) {
            return true;
        }
    }
    return false;
}

bool Boundary::compute_validity(bool allow_boundaries_between_opposite_labels, const CustomMeshHalfedges& mesh_halfedges) {
    if (axis==-1) {
        if(contains_lower_than_180_degrees_angles(mesh_halfedges)) {
            is_valid = false;
        }
        else if (allow_boundaries_between_opposite_labels) {
            is_valid = true;
        }
        else {
            is_valid = false;
        }
    }
    else {
        is_valid = true;
    }
    return is_valid;
}

bool Boundary::find_turning_points(const CustomMeshHalfedges& mesh_halfedges) {

    // based on https://github.com/LIHPC-Computational-Geometry/evocube/blob/master/src/graphcut_labeling.cpp#L142 graphcutTurningPoints()

    geo_assert(!halfedges.empty());
    std::size_t nb_elem = halfedges.size(); // there are one vertex less than the number of edges
    // elem n°0 is local vertex 0, which is at the base (origin) of halfedge 0, and so on

    // if boundary of 1 edge only -> no turning point
    if(nb_elem == 1) return false;

    // if axis==-1, no direction to compare with edges
    if (axis == -1) return false;

    std::vector<int> per_vertex_result(nb_elem);

    // unary costs
    std::vector<int> data_cost(nb_elem*2); // 2 "labels" : is a turning point, or is not a turning point
    vec3 desired_direction = label2vector[axis*2]; // vector corresponding to the axis

    FOR(i,nb_elem) {
        vec3 edge =  normalize(halfedge_vector(mesh_halfedges.mesh(),halfedges[i])); // edge from vertex i to vertex i+1 -> halfedge i
        double dot = (GEO::dot(desired_direction,edge))/0.9;
        double cost = 1.0 - std::exp(-(1./2.)*std::pow(dot,2));
        if (GEO::dot(desired_direction,edge) > 0){
            data_cost[i*2 + 0] = (int) (100.0 * cost);
        }
        else {
            data_cost[i*2 + 1] = (int) (100.0 * cost);
        }
    }

	// binary cost coefficients
	std::vector<int> smooth_cost(2*2);
	for(unsigned int l1 = 0; l1 < 2; l1++ ){
		for(unsigned int l2 = 0; l2 < 2; l2++ ){ 
			if (l1==l2) smooth_cost[l1+l2*2] = 0;
            else smooth_cost[l1+l2*2] = 1;
        }
    }

    try{
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph( (GCoptimization::SiteID) nb_elem,2);
        gc->setDataCost(data_cost.data());
        gc->setSmoothCost(smooth_cost.data());
		
        FOR(i,nb_elem-1) {
            vec3 edge1 = normalize(halfedge_vector(mesh_halfedges.mesh(),halfedges[i])); // edge from vertex i to vertex i+1 -> halfedge i
            vec3 edge2 = normalize(halfedge_vector(mesh_halfedges.mesh(),halfedges[i+1])); // edge from vertex i+1 to vertex i+2 -> halfedge i+1
            double dot = (GEO::dot(edge1,edge2) - 1.0)/FALLOFF_BINARY;
            double cost = std::exp(-(1./2.0)*std::pow(dot,2));
            gc->setNeighbors( (GCoptimization::SiteID) i, (GCoptimization::SiteID) i+1, (int) (100.0*cost));
        }

        gc->expansion(2);// run expansion for 2 iterations
		
        for ( std::size_t  i = 0; i < nb_elem; i++ )
			per_vertex_result[i] = gc->whatLabel( (GCoptimization::SiteID) i);

		delete gc;
	}
	catch (GCException e){
		e.Report();
	}

    // parse result for turning-points
    for (unsigned int i=1; i<per_vertex_result.size(); i++){
        if (per_vertex_result[i] != per_vertex_result[i-1]){
            turning_points.push_back(TurningPoint{});
            turning_points.back().fill_from((index_t) i, *this, mesh_halfedges);
        }
    }

    // handle loops
    if ( (start_corner == end_corner) &&
         ((*per_vertex_result.begin()) != (*per_vertex_result.rbegin())) ) {
        turning_points.push_back(TurningPoint{});
        turning_points.back().fill_from((index_t) 0, *this, mesh_halfedges); // the first (& last) vertex is a turning point (loop)
    }

    return !turning_points.empty(); // true if there are turning-points, else false
}

index_t Boundary::chart_at_other_side(index_t origin_chart) const {
    geo_assert(left_chart != index_t(-1));
    geo_assert(right_chart != index_t(-1));
    if(left_chart == origin_chart) {
        return right_chart;
    }
    else if(right_chart == origin_chart) {
        return left_chart;
    }
    else {
        geo_assert_not_reached;
    }
}

index_t Boundary::other_label(const std::vector<Chart>& charts, index_t label) {
    if(label == charts[left_chart].label) {
        return charts[right_chart].label;
    }
    else if(label == charts[right_chart].label) {
        return charts[left_chart].label;
    }
    else {
        geo_assert_not_reached;
    }
}

void Boundary::print_successive_halfedges(fmt::v9::ostream& out, Mesh& mesh) {
    // example:
    //          
    // 1022 --{124,1054}--> 54684 --{45348,4806}--> (7654) --{23347,1987}--> 8541
    // 
    // 1022         = vertex index
    // --{124,1054}-->  = halfedge (facet index then corner index)
    // (7654)       = vertex index on which there is a turning point
    //
    out.print("\t ");
    FOR(he_index,halfedges.size()) {
        // print vertex at the beginning of current halfedge
        if(halfedge_has_turning_point_at_base(he_index)) {
            // there is a turning point on this vertex
            out.print("({}) ",mesh.facet_corners.vertex(halfedges[he_index].corner));
        }
        else {
            out.print("{} ",mesh.facet_corners.vertex(halfedges[he_index].corner));
        }
        out.print("--{{{},{}}}--> ",halfedges[he_index].facet,halfedges[he_index].corner);
    }
    // print last vertex of the boundary
    CustomMeshHalfedges mesh_halfedges(mesh);
    out.print("{}\n",Geom::halfedge_vertex_index_to(mesh,*(halfedges.rbegin()))); // get vertex at extremity of last halfedge = last vertex
}

index_t Boundary::turning_point_vertex(index_t turning_point_index, const Mesh& mesh) const {
    return mesh.facet_corners.vertex(halfedges[turning_points[turning_point_index].outgoing_local_halfedge_index_].corner);
}

bool Boundary::halfedge_has_turning_point_at_base(index_t local_halfedge_index) const {
    for(auto& tp : turning_points) {
        if(tp.outgoing_local_halfedge_index_ == local_halfedge_index) {
            return true;
        }
    }
    return false;
}

index_t Boundary::get_closest_boundary_of_turning_point(const TurningPoint& turning_point, index_t closest_corner, const CustomMeshHalfedges& mesh_he, const std::map<MeshHalfedges::Halfedge,std::pair<index_t,bool>,HalfedgeCompare>& halfedge2boundary, const std::vector<Corner>& corners) const {
    geo_assert(mesh_he.is_using_facet_region()); // assert a facet region is defined in mesh_he

    // find the closest corner if not already done in parent function
    if(closest_corner == index_t(-1)) {
        closest_corner = turning_point.get_closest_corner(*this,mesh_he);
    }
    
    MeshHalfedges::Halfedge current_halfedge;
    // TODO check association between closest corner, turning point direction and direction around vertex (clockwise/counterclockwise)
    if(closest_corner == start_corner) {
        current_halfedge = halfedges[0]; // get first halfedge of the boundary. its origin vertex is on start_corner
        geo_assert(halfedge_vertex_index_from(mesh_he.mesh(),current_halfedge) == corners[closest_corner].vertex);
        if(turning_point.is_towards_left()) {
            mesh_he.move_counterclockwise_around_vertex_until_on_border(current_halfedge);
        }
        else { // turning point towards right
            mesh_he.move_clockwise_around_vertex_until_on_border(current_halfedge);
        }
    }
    else if(closest_corner == end_corner) {
        current_halfedge = *halfedges.rbegin(); // get last halfedge of the boundary. its extremity vertex is on end_corner
        mesh_he.move_to_opposite(current_halfedge); // flip it so that its origin vertex is on end_corner
        geo_assert(halfedge_vertex_index_from(mesh_he.mesh(),current_halfedge) == corners[closest_corner].vertex);
        if(turning_point.is_towards_right()) {
            mesh_he.move_counterclockwise_around_vertex_until_on_border(current_halfedge);
        }
        else { // turning point towards left
            mesh_he.move_clockwise_around_vertex_until_on_border(current_halfedge);
        }
    }
    else {
        geo_assert_not_reached;
    }
    geo_assert(halfedge_vertex_index_from(mesh_he.mesh(),current_halfedge) == corners[closest_corner].vertex);
    geo_assert(halfedge2boundary.contains(current_halfedge)); // assert `current_halfedge` is on a boundary
    geo_assert(halfedge2boundary.at(current_halfedge).first != halfedge2boundary.at(halfedges[0]).first); // assert the boundary found is not *this
    return halfedge2boundary.at(current_halfedge).first;
}

vec3 Boundary::vector_between_corners(const Mesh& mesh, const std::vector<Corner>& corners) const {
    geo_assert(start_corner != index_t(-1));
    geo_assert(end_corner != index_t(-1));
    // compute the vector going from the start corner to the end corner of the boundary
    return mesh_vertex(mesh,corners[end_corner].vertex) - mesh_vertex(mesh,corners[start_corner].vertex);
}

vec3 Boundary::average_vector_between_corners(const Mesh& mesh, const std::vector<Corner>& corners) const {
    return vector_between_corners(mesh,corners) / (double) halfedges.size();
}

void Boundary::get_adjacent_facets(const Mesh& mesh, std::set<index_t>& adjacent_facets, BoundarySide boundary_side_to_explore, const std::vector<index_t>& facet2chart, size_t max_dist) const {
    if(max_dist > 1) {
        fmt::println("Not implemented");
        geo_assert_not_reached;
    }
    index_t adjacent_facet = index_t(-1);
    for(const auto& current_halfedge : halfedges) {
        if( (boundary_side_to_explore == LeftAndRight) || (boundary_side_to_explore == OnlyLeft) ) {
            adjacent_facet = halfedge_facet_left(mesh,current_halfedge);
            adjacent_facets.insert(adjacent_facet);
            if(max_dist == 1) {
                get_adjacent_facets_conditional(mesh,adjacent_facet,left_chart,facet2chart,adjacent_facets);
            }
        }
        if( (boundary_side_to_explore == LeftAndRight) || (boundary_side_to_explore == OnlyRight) ) {
            adjacent_facet = halfedge_facet_right(mesh,current_halfedge);
            adjacent_facets.insert(adjacent_facet);
            if(max_dist == 1) {
                get_adjacent_facets_conditional(mesh,adjacent_facet,right_chart,facet2chart,adjacent_facets);
            }
        }
    }
}

void Boundary::get_flipped(const MeshHalfedges& mesh_he, Boundary& flipped_boundary) const {
    flipped_boundary.axis = axis;
    flipped_boundary.is_valid = is_valid;
    flipped_boundary.on_feature_edge = on_feature_edge;
    flipped_boundary.average_normal = average_normal;
    flipped_boundary.left_chart = right_chart;
    flipped_boundary.right_chart = left_chart;
    flipped_boundary.start_corner = end_corner;
    flipped_boundary.end_corner = start_corner;
    // go through halfedges in the reverse order
    flipped_boundary.halfedges.clear();
    for(auto it = halfedges.rbegin(); it != halfedges.rend(); ++it) {
        MeshHalfedges::Halfedge current_halfedge = *it;
        mesh_he.move_to_opposite(current_halfedge);
        flipped_boundary.halfedges.push_back(current_halfedge);
    }
    geo_assert(flipped_boundary.halfedges.size() == halfedges.size());
    // `turning_points` store indices in `halfedges` where a turning poin is on the origin vertex
    for(auto it = turning_points.rbegin(); it != turning_points.rend(); ++it) {
        TurningPoint current_turning_point = *it;
        geo_assert(current_turning_point.outgoing_local_halfedge_index_ < halfedges.size());
        // modify `current_turning_point` so it is relative to the flipped boundary
        current_turning_point.is_towards_right_ = !current_turning_point.is_towards_right_;
        current_turning_point.outgoing_local_halfedge_index_ = (index_t) halfedges.size()-current_turning_point.outgoing_local_halfedge_index_;
        geo_assert(current_turning_point.outgoing_local_halfedge_index_ < flipped_boundary.halfedges.size());
        flipped_boundary.turning_points.push_back(current_turning_point);
    }
}

void Boundary::split_at_turning_point(const MeshHalfedges& mesh_he, Boundary& downward_boundary, Boundary& upward_boundary) const {
    /* Given `this` being a boundary with a single turning-point:
     *
     *   start_corner         turning-point              end_corner        
     *       X ---> ---> ---> ---> o ---> ---> ---> ---> ---> X
     *
     * Returns :
     *  - `downward_boundary`:
     *       X <--- <--- <--- <--- X
     *  - `upward_boundary`:
     *                             X ---> ---> ---> ---> ---> X
     */
    geo_assert(turning_points.size() == 1); // cannot split boundary in 2 parts if there is 0, or more than 1 turning-points
    if(turning_points[0].outgoing_local_halfedge_index_ == 0) {
        fmt::println(Logger::err("labeling graph"),"Cannot use split_at_turning_point() on a boundary with a turning-point on the first vertex (cyclic boundary?)"); Logger::err("labeling graph").flush();
        geo_assert_not_reached;
    }

    downward_boundary.axis = axis;
    upward_boundary.axis = axis;
    downward_boundary.is_valid = is_valid;
    upward_boundary.is_valid = is_valid;
    downward_boundary.on_feature_edge = on_feature_edge; // should be recomputed on `downward_boundary` only
    upward_boundary.on_feature_edge = on_feature_edge; // should be recomputed on `upward_boundary` only
    downward_boundary.average_normal = average_normal; // should be recomputed on `downward_boundary` only
    upward_boundary.average_normal = average_normal; // should be recomputed on `upward_boundary` only
    downward_boundary.left_chart = right_chart;
    upward_boundary.left_chart = left_chart;
    downward_boundary.right_chart = left_chart;
    upward_boundary.right_chart = right_chart;
    downward_boundary.start_corner = index_t(-1); // no corner
    upward_boundary.start_corner = index_t(-1); // no corner
    downward_boundary.end_corner = start_corner;
    upward_boundary.end_corner = end_corner;
    downward_boundary.turning_points.clear();
    upward_boundary.turning_points.clear();
    downward_boundary.halfedges.clear();
    upward_boundary.halfedges.clear();
    for(signed_index_t he_index = (signed_index_t) turning_points[0].outgoing_local_halfedge_index_-1; he_index >= 0; --he_index) {
        downward_boundary.halfedges.push_back(halfedges[(index_t) he_index]);
        mesh_he.move_to_opposite(downward_boundary.halfedges.back());
    }
    for(index_t he_index = turning_points[0].outgoing_local_halfedge_index_; he_index < halfedges.size(); ++he_index) {
        upward_boundary.halfedges.push_back(halfedges[he_index]);
    }
    geo_assert(downward_boundary.halfedges.size() + upward_boundary.halfedges.size() == halfedges.size());
}

void Boundary::per_edges_axis_assignement_cost(const Mesh& mesh, index_t axis, std::vector<double>& costs) const {
    geo_assert(axis < 3); // must be either 0=X, 1=Y or 2=Z
    costs.resize(halfedges.size());
    FOR(he_index,halfedges.size()) {
        costs[he_index] = angle_between_halfedge_and_axis(mesh,halfedges[he_index],axis) / std::numbers::pi;
    }
}

void Boundary::per_edges_cumulative_axis_assignement_cost(const Mesh& mesh, index_t axis, std::vector<double>& costs, bool accumulation_from_start_corner) const {
    per_edges_axis_assignement_cost(mesh,axis,costs);
    if(accumulation_from_start_corner) {
        for(index_t he_index = 1; he_index < halfedges.size(); ++he_index) {
            costs[he_index] += costs[he_index-1];
        }
    }
    else { // accumulation from end corner
        for(signed_index_t he_index = (signed_index_t) halfedges.size()-2; he_index >= 0; --he_index) {
            costs[(index_t) he_index] += costs[(index_t) he_index+1];
        }
    }
}

std::ostream& operator<< (std::ostream &out, const Boundary& data) {
    fmt::println(out,"\taxis : {}",data.axis);
    fmt::println(out,"\tis_valid : {}",data.is_valid);
    fmt::println(out,"\tturning_points : ");
    for(const auto tp : data.turning_points) {
        fmt::print(out,"{}",tp);
    }
    fmt::println(out,"\ton_feature_edge : {}",data.on_feature_edge);
    fmt::println(out,"\thalfedges : {}",data.halfedges);
    fmt::println(out,"\tleft_chart : {}",OPTIONAL_TO_STRING(data.left_chart));
    fmt::println(out,"\tright_chart : {}",OPTIONAL_TO_STRING(data.right_chart));
    fmt::println(out,"\tstart_corner : {}",OPTIONAL_TO_STRING(data.start_corner));
    fmt::println(out,"\tend_corner : {}",OPTIONAL_TO_STRING(data.end_corner));
    
    return out;
}

void StaticLabelingGraph::fill_from(Mesh& mesh, std::string facet_attribute, bool allow_boundaries_between_opposite_labels, const std::set<std::pair<index_t,index_t>>& feature_edges) {

    // based on https://github.com/LIHPC-Computational-Geometry/genomesh/blob/main/src/flagging.cpp#L795
    // but here we use a disjoint-set for step 1

    // expecting a surface triangle mesh
    geo_assert(mesh.facets.are_simplices());
    geo_assert(mesh.cells.nb()==0);

    clear();
    allow_boundaries_between_opposite_labels_ = allow_boundaries_between_opposite_labels;

    facet2chart.resize(mesh.facets.nb()); // important: memory allocation allowing to call ds.getSetsMap() on the underlying array
    vertex2corner.resize(mesh.vertices.nb(),index_t(-1)); // contrary to facet2chart where all facets are associated to a chart, not all vertices are associated to a corner -> use UNDEFINED for vertices that are not on a corner
    CustomMeshHalfedges mesh_half_edges_(mesh); // Half edges API

    // STEP 1 : Aggregate adjacent triangles of same label as chart

    Attribute<index_t> labeling(mesh.facets.attributes(),facet_attribute); // read and store the facets attribute corresponding to the labeling

    DisjointSet<index_t> ds(labeling.size()); // create a disjoint-set data structure to group facets by labels
    for(index_t f: mesh.facets) { // for each facet of the mesh
        for(index_t le = 0; le < 3; ++le) { // for each local edge of facet f in {0,1,2}
            index_t adjacent_facet = mesh.facets.adjacent(f,le); // get the facet id of the facet beyond le
            if(labeling[f] == labeling[adjacent_facet]) { // if facets f and adjacent_facet have the same label
                ds.mergeSets(f,adjacent_facet); // merge facets f and adjacent_facet
            }
        }
    }
    std::size_t nb_charts = ds.getSetsMap(facet2chart.data()); // get the map (facet id -> chart id) and the number of charts
    charts.resize(nb_charts);

    // fill the Chart objects
    for(index_t f: mesh.facets) { // for each facet of the mesh
        Chart& current_chart = charts[facet2chart[f]]; // get the chart associated to this facet
        current_chart.label = labeling[f]; // (re)define its label
        current_chart.facets.emplace(f); // register the facet
    }

    // STEP 2 : Explore corners and boundaries
    // Ideally, we have to iterate over all vertices, and,
    // for each of them, explore adjacent facets to compute valence (number of boundary edges).
    // But I don't know how to get the adjacent facets of a given vertex.
    // What I can do, with half-edges, is to iterate over all facet corners, get the vertex at this corner,
    // and go around the vertex ring with CustomMeshHalfedges::move_to_next_around_vertex(). See VertexRingWithBoundaries::explore()
    // This works for shapes having several solids connected by a vertex only

    mesh_half_edges_.set_use_facet_region(facet_attribute); // indicate to Geogram the attribute with the charts (= regions for Geogram) from which we want to find the boundaries (= borders for Geogram)
    std::vector<MeshHalfedges::Halfedge> boundary_edges_to_explore; // vector of boundary edges we encountered, and that must be explored later

    index_t current_vertex = index_t(-1);
    for(index_t f: mesh.facets) { for(index_t c: mesh.facets.corners(f)) { // for each facet corner (f,c)
        
        MeshHalfedges::Halfedge halfedge(f,c); // halfedge on facet f having corner c as base
        current_vertex = mesh.facet_corners.vertex(halfedge.corner); // get the vertex at the base of halfedge
        
        // (re)explore this vertex, because maybe it was explored from another solid of the same shape (multiple vertex rings)
        VertexRingWithBoundaries current_vertex_ring;
        current_vertex_ring.explore(halfedge,mesh_half_edges_);
        if(current_vertex_ring.valence() < 3) {
            // Not a corner. At least, from the point of view of this vertex ring
            continue;            
        }
        // else : we found a corner !

        if(vertex2corner[current_vertex] != index_t(-1)) { // if a corner already exists on this vertex
            if(corners[vertex2corner[current_vertex]].halfedge_is_in_boundary_edges(halfedge)) { // if the current halfedge is already referenced in this corner
                continue; // nothing to do, skip to next vertex
            }
            else {
                // else : the corner on this vertex has never seen the current halfedge : add the boundary edges encountered to the boundary edges to explore
                boundary_edges_to_explore += current_vertex_ring.boundary_edges;
            }
        }
        else {
            // we found an unexplored corner
            corners.push_back(Corner()); // create a new Corner
            corners.back().vertex = current_vertex; // link it to the current vertex
            corners.back().vertex_rings_with_boundaries.push_back(current_vertex_ring); // add it the just-explored vertex ring
            vertex2corner[current_vertex] = (index_t) index_of_last(corners); // link the boundary to this corner
            boundary_edges_to_explore += current_vertex_ring.boundary_edges; // the boundary edges around this corner must be explored
        }        
        
        while(!boundary_edges_to_explore.empty()) { // while there still are boundaries to explore

            halfedge = boundary_edges_to_explore.back(); // pick an halfedge from the vector
            boundary_edges_to_explore.pop_back(); // remove this halfedge from the vector

            // Check if this halfedge has already been explored since its insertion in the vector
            if(MAP_CONTAINS(halfedge2boundary,halfedge)) {
                // this halfedge is already linked to a boundary, no need to re-explore the boundary
                continue;
            }

            boundaries.push_back(Boundary()); // create a new boundary
            boundaries.back().start_corner = vertex2corner[mesh.facet_corners.vertex(halfedge.corner)]; // link this boundary to the corner at the beginning
            // explore it, edge by edge
            boundaries.back().explore(halfedge,
                                       mesh_half_edges_,
                                       (index_t) index_of_last(boundaries), // get the index of this boundary
                                       feature_edges,
                                       facet2chart,
                                       vertex2corner,
                                       charts,
                                       corners,
                                       halfedge2boundary,
                                       boundary_edges_to_explore);
            if(boundaries.back().find_turning_points(mesh_half_edges_)) {
                non_monotone_boundaries.push_back((index_t) index_of_last(boundaries));
            }
            if(boundaries.back().compute_validity(allow_boundaries_between_opposite_labels,mesh_half_edges_)==false) {
                invalid_boundaries.push_back((index_t) index_of_last(boundaries));
            }
        }
    }}

    // STEP 3 : Find boundaries with no corners
    // explore all half edges
    // if not linked to a boundary, look at neighboring labels
    // if they are different, explore this boundary
    for(index_t f: mesh.facets) { for(index_t c: mesh.facets.corners(f)) {
        MeshHalfedges::Halfedge halfedge(f,c); // halfedge on facet f having corner c as base

        if (MAP_CONTAINS(halfedge2boundary,halfedge)) {
            // halfedge is on an already explorer boundary
            continue;
        }

        if (!mesh_half_edges_.halfedge_is_border(halfedge)) {
            // halfedge is inside a chart
            continue;
        }

        boundaries.push_back(Boundary()); // create a new boundary
        boundaries.back().start_corner = index_t(-1); // no corner at the beginning
        // explore it, edge by edge
        boundaries.back().explore(halfedge,
                                  mesh_half_edges_,
                                  (index_t) index_of_last(boundaries), // get the index of this boundary
                                  feature_edges,
                                  facet2chart,
                                  vertex2corner,
                                  charts,
                                  corners,
                                  halfedge2boundary,
                                  boundary_edges_to_explore);
        if(boundaries.back().find_turning_points(mesh_half_edges_)) {
            non_monotone_boundaries.push_back((index_t) index_of_last(boundaries));
        }
        if(boundaries.back().compute_validity(allow_boundaries_between_opposite_labels,mesh_half_edges_)==false) {
            invalid_boundaries.push_back((index_t) index_of_last(boundaries));
        }
    }}

    // STEP 4 : Find invalid charts

    Attribute<bool> facet_on_invalid_chart(mesh.facets.attributes(),"on_invalid_chart");
    facet_on_invalid_chart.fill(false);
    FOR(chart_index,charts.size()) {
        if(charts[chart_index].boundaries.size() <= 3) { // if this chart has 3 neighbors or less
            invalid_charts.push_back(chart_index); // register this chart

            for(index_t facet: charts[chart_index].facets) { // for all facets in this chart
                facet_on_invalid_chart[facet] = true; // mark them
            }
        }
    }

    // STEP 5 : Find invalid corners
    FOR(c,corners.size()) {
        if(corners[c].compute_validity(allow_boundaries_between_opposite_labels,boundaries,halfedge2boundary)==false) {
            invalid_corners.push_back(c);
        }
    }
}

void StaticLabelingGraph::clear() {
    charts.clear();
    boundaries.clear();
    corners.clear();
    facet2chart.clear();
    halfedge2boundary.clear();
    vertex2corner.clear();
    invalid_charts.clear();
    invalid_boundaries.clear();
    invalid_corners.clear();
    non_monotone_boundaries.clear();
    allow_boundaries_between_opposite_labels_ = false;
    // note : the mesh attributes will not be removed, because we no longer have a ref/pointer to the mesh
}

bool StaticLabelingGraph::is_valid() {
    return (nb_invalid_charts()==0) && 
		   (nb_invalid_boundaries()==0) && 
		   (nb_invalid_corners()==0);
}

std::size_t StaticLabelingGraph::nb_charts() const {
    return charts.size();
}

std::size_t StaticLabelingGraph::nb_boundaries() const {
    return boundaries.size();
}

std::size_t StaticLabelingGraph::nb_corners() const {
    return corners.size();
}

std::size_t StaticLabelingGraph::nb_facets() const {
    return facet2chart.size();
}

std::size_t StaticLabelingGraph::nb_vertices() const {
    return vertex2corner.size();
}

std::size_t StaticLabelingGraph::nb_invalid_charts() const {
    return invalid_charts.size();
}

std::size_t StaticLabelingGraph::nb_invalid_boundaries() const {
    return invalid_boundaries.size();
}

std::size_t StaticLabelingGraph::nb_invalid_corners() const {
    return invalid_corners.size();
}

std::size_t StaticLabelingGraph::nb_non_monotone_boundaries() const {
    return non_monotone_boundaries.size();
}

std::size_t StaticLabelingGraph::nb_turning_points() const {
    std::size_t count = 0;
    for(index_t i : non_monotone_boundaries) {
        count += boundaries[i].turning_points.size();
    }
    return count;
}

bool StaticLabelingGraph::is_allowing_boundaries_between_opposite_labels() const {
    return allow_boundaries_between_opposite_labels_;
}

bool StaticLabelingGraph::vertex_is_only_surrounded_by(index_t vertex_index, std::vector<index_t> expected_charts, const std::vector<std::vector<index_t>>& vertex_to_adj_facets) const {
    geo_assert(!vertex_to_adj_facets.empty());
    for(auto adj_facet : vertex_to_adj_facets.at(vertex_index)) { // for each adjacent facet of `vertex_index`
        if (!VECTOR_CONTAINS(expected_charts,facet2chart.at(adj_facet))) { // if the chart index of this facet is not among the `expected_charts`
            return false;
        }
    }
    return true;
}

void StaticLabelingGraph::get_adjacent_charts_of_vertex(index_t vertex_index, const std::vector<std::vector<index_t>>& vertex_to_adj_facets, std::set<index_t>& adjacent_charts) const {
    geo_assert(vertex_index < vertex_to_adj_facets.size());
    adjacent_charts.clear();
    for(index_t adj_facet : vertex_to_adj_facets[vertex_index]) {
        geo_assert(adj_facet < facet2chart.size());
        adjacent_charts.insert(facet2chart[adj_facet]);
    }
}

void StaticLabelingGraph::dump_to_text_file(const char* filename, Mesh& mesh) {
    auto out = fmt::output_file(filename);
    out.print("{}",(*this));

    // pretty print boundaries
    for(std::size_t boundary_index = 0; boundary_index < nb_boundaries(); ++boundary_index) {
        out.print("boundaries[{}]\n",boundary_index);
        boundaries[boundary_index].print_successive_halfedges(out,mesh);
    }
}

void StaticLabelingGraph::dump_to_D3_graph(const char* filename) {
    nlohmann::json graph;
    // export nodes. group 1 = chart, 2 = boundary, 3 = corner
    FOR(chart_index,nb_charts()) {
        graph["nodes"] += nlohmann::json::object({{"id", fmt::format("chart{}",chart_index)}, {"group", 1}});
    }
    FOR(boundary_index,nb_boundaries()) {
        graph["nodes"] += nlohmann::json::object({{"id", fmt::format("boundary{}",boundary_index)}, {"group", 2}});
    }
    FOR(corner_index,nb_corners()) {
        graph["nodes"] += nlohmann::json::object({{"id", fmt::format("corner{}",corner_index)}, {"group", 3}});
    }
    // export links
    for(const auto& b : boundaries) {
        // link this boundary to the left and right charts
        graph["links"] += nlohmann::json::object({{"source", fmt::format("chart{}",b.left_chart)}, {"target", fmt::format("chart{}",b.right_chart)}, {"value", 1}});
        // link this boundary to the start and end corner
        graph["links"] += nlohmann::json::object({{"source", fmt::format("corner{}",b.start_corner)}, {"target", fmt::format("corner{}",b.end_corner)}, {"value", 1}});
    }
    std::fstream ofs(filename,std::ios_base::out);
    if(ofs.good()) {
        ofs << std::setw(4) << graph << std::endl;
    }
    else {
        fmt::println(Logger::err("I/O"),"Cannot write into {}",filename); Logger::err("I/O").flush();
    }
}

bool StaticLabelingGraph::is_turning_point(const Mesh& mesh, index_t vertex_index) const {
    for(auto b : non_monotone_boundaries) {
        FOR(tp,boundaries[b].turning_points.size()) {
            if(boundaries[b].turning_point_vertex(tp,mesh) == vertex_index) {
                return true;
            }
        }
    }
    return false;
}

std::ostream& operator<< (std::ostream &out, const StaticLabelingGraph& data) {

    // write charts

    for(std::size_t chart_index = 0; chart_index < data.nb_charts(); ++chart_index) {
        fmt::println(out,"chart[{}]",chart_index);
        fmt::println(out,"{}",data.charts[chart_index]);
    }
    if(data.nb_charts()==0) fmt::println(out,"no charts");

    // write boundaries

    for(std::size_t boundary_index = 0; boundary_index < data.nb_boundaries(); ++boundary_index) {
        fmt::println(out,"boundaries[{}]",boundary_index);
        fmt::println(out,"{}",data.boundaries[boundary_index]);
    }
    if(data.nb_boundaries()==0) fmt::println(out,"no boundaries");

    // write corners

    for(std::size_t corner_index = 0; corner_index < data.nb_corners(); ++corner_index) {
        fmt::println(out,"corners[{}]",corner_index);
        fmt::println(out,"{}",data.corners[corner_index]);
    }
    if(data.nb_corners()==0) fmt::println(out,"no corners");

    // write facet2chart
    
    fmt::println(out,"facet2chart");
    for(std::size_t f = 0; f < data.nb_facets(); ++f) {
        fmt::println(out,"\t[{}] {}",f,data.facet2chart[f]);
    }

    // write halfedge2boundary

    fmt::println(out,"halfedge2boundary");
    {
        for(auto entry : data.halfedge2boundary) {
            fmt::println(out,"\t[{}] {}",entry.first,entry.second);
        }
    }

    // write vertex2corner
    
    fmt::println(out,"vertex2corner");
    for(std::size_t v = 0; v < data.nb_vertices(); ++v) {
        fmt::println(out,"\t[{}] {}",v,OPTIONAL_TO_STRING(data.vertex2corner[v]));
    }

    // write invalid_charts
    
    fmt::println(out,"invalid_charts");
    FOR(i,data.invalid_charts.size()) {
        fmt::println(out,"\t[{}] {}",i,data.invalid_charts[i]);
    }

    // write invalid_boundaries
    
    fmt::println(out,"invalid_boundaries");
    FOR(i,data.invalid_boundaries.size()) {
        fmt::println(out,"\t[{}] {}",i,data.invalid_boundaries[i]);
    }

    // write invalid_corners
    
    fmt::println(out,"invalid_corners");
    FOR(i,data.invalid_corners.size()) {
        fmt::println(out,"\t[{}] {}",i,data.invalid_corners[i]);
    }

    // write non_monotone_boundaries
    
    fmt::println(out,"non_monotone_boundaries");
    FOR(i,data.non_monotone_boundaries.size()) {
        fmt::println(out,"\t[{}] {}",i,data.non_monotone_boundaries[i]);
    }

    // write allow_boundaries_between_opposite_labels_

    fmt::println(out,"allow_boundaries_between_opposite_labels_ = {}",data.allow_boundaries_between_opposite_labels_);

    return out;
}
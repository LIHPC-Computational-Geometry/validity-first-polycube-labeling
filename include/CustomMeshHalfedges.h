// # Changelog - https://keepachangelog.com
//
// Modifications of Geogram's source files
//   ext/geogram/src/lib/geogram/mesh/mesh_halfedges.h
// & ext/geogram/src/lib/geogram/mesh/mesh_halfedges.cpp
//
// ## Modified move operators
//
// ### Added
//
// - new function : custom_move_to_next_around_border()
//
// ### Changed
//
// - move_to_next_around_vertex() and move_to_prev_around_vertex() are independant of facet_region_ if ignore_borders==true
//
// ## To have an independent copy of MeshHalfedges
//
// ### Added
//
// - using namespace GEO
//
// ### Changed
//
// - replace MeshHalfedges by CustomMeshHalfedges
// - Geogram's include guard replaced by #pragma once
//
// ### Removed
//
// - namespace GEO and Geom wrapping the code
// - inclusion of <geogram/mesh/mesh_halfedges.h>

#pragma once

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <iostream>

using namespace GEO;
using namespace GEO::Geom;

    /**
     * \brief Exposes a half-edge like API for
     * traversing a Mesh.
     */
    class GEOGRAM_API CustomMeshHalfedges {
    public:
        /**
         * \brief Stores a reference to a mesh corner and facet, and
         *  provides a halfedge-like API.
         */
        struct Halfedge {

            static const index_t NO_FACET  = index_t(-1);
            static const index_t NO_CORNER = index_t(-1);

            /**
             * \brief Constructs a new uninitialized Halfedge.
             */
            Halfedge() :
                facet(NO_FACET),
                corner(NO_CORNER) {
            }

            /**
             * \brief Constructs a new Halfedge from a facet and corner index.
             * \param[in] f the facet index
             * \param[in] c the corner index
             */
            Halfedge(index_t f, index_t c) :
                facet(f),
                corner(c) {
            }

            /**
             * \brief Clears this Halfedge.
             */
            void clear() {
                facet = NO_FACET;
                corner = NO_CORNER;
            }

            /**
             * \brief Tests whether this Halfedge is initialized.
             * \return true if this Halfedge is uninitialized, false otherwise
             */
            bool is_nil() const {
                return (facet == NO_FACET) && (corner == NO_CORNER);
            }

            /**
             * \brief Tests whether this Halfedge is the same as another one
             * \param[in] rhs the comparand
             * \return true if this Halfedge and \p rhs refer to the same
             *  facet and corner, false otherwise.
             */
            bool operator== (const Halfedge& rhs) const {
                return facet == rhs.facet && corner == rhs.corner;
            }

            /**
             * \brief Tests whether this Halfedge is different from another one
             * \param[in] rhs the comparand
             * \return true if this Halfedge and \p rhs refer to a different
             *  facet or corner, false otherwise.
             */
            bool operator!= (const Halfedge& rhs) const {
                return !(rhs == *this);
            }

            index_t facet;
            index_t corner;

        };

        /**
         * \brief Creates a new CustomMeshHalfedges
         * \param[in] mesh the Mesh
         */
        CustomMeshHalfedges(Mesh& mesh) : mesh_(mesh) {
        }

        /**
         * \brief Gets the mesh.
         * \return a reference to the mesh.
         */
        Mesh& mesh() {
            return mesh_;
        }

        /**
         * \brief Gets the mesh.
         * \return a const reference to the mesh.
         */
        const Mesh& mesh() const {
            return mesh_;
        }

        /**
         * \brief Sets whether facet regions determine borders.
         * \param[in] x if set, then an halfedge incident to two facets
         *  with different facet regions is considered to be a
         *  border
         */
        void set_use_facet_region(bool x) {
            if(x) {
                if(!facet_region_.is_bound()) {
                    facet_region_.bind(mesh_.facets.attributes(),"region");
                }
            } else {
                if(facet_region_.is_bound()) {
                    facet_region_.unbind();
                }
            }
        }

	/**
	 * \brief Sets a facet attribute name that determines borders.
	 * \param[in] attribute_name the name of the facet attribute to
	 *  be used to determine borders.
	 */
	void set_use_facet_region(const std::string& attribute_name) {
	    if(facet_region_.is_bound()) {
		facet_region_.unbind();
	    }
	    facet_region_.bind(mesh_.facets.attributes(),attribute_name);
	}

	
        /**
         * \brief Tests whether a Halfedge is valid.
         * \param[in] H the Halfedge to be tested
         * \return true if \p H refers to a halfedge that
         *  exists in the mesh, false otherwise
         * \note It only tests whether H.corner and H.facet are valid
         *  indices in the mesh, but does not test whether H.corner exists
         *  in H.facet.
         */
        bool halfedge_is_valid(const Halfedge& H) const {
            return
                H.facet != Halfedge::NO_FACET && 
                H.corner != Halfedge::NO_CORNER &&
                H.facet < mesh_.facets.nb() &&
                H.corner < mesh_.facet_corners.nb()
            ;
        }

        /**
         * \brief Tests whether a Halfedge is on the boder.
         * \details If set_use_facet_region() is set, then
         *  Halfedges incident to two different facet regions are
         *  considered as borders.
         * \param[in] H the Halfedge
         * \return true if \p H is on the border, false otherwise
         */
        bool halfedge_is_border(const Halfedge& H) const {
            geo_debug_assert(halfedge_is_valid(H));
            if(facet_region_.is_bound()) {
                index_t f = H.facet;
                index_t adj_f =
                    mesh_.facet_corners.adjacent_facet(H.corner);
                return
                    adj_f == NO_FACET ||
                    facet_region_[f] != facet_region_[adj_f]
                ;
            } 
            return mesh_.facet_corners.adjacent_facet(H.corner) == NO_FACET;
        }

        /**
         * \brief Replaces a Halfedge with the next one around the facet.
         * \param[in,out] H the Halfedge
         */
        void move_to_next_around_facet(Halfedge& H) const {
            geo_debug_assert(halfedge_is_valid(H));
            H.corner = mesh_.facets.next_corner_around_facet(H.facet, H.corner);
        }

        /**
         * \brief Replaces a Halfedge with the previous one around the facet.
         * \param[in,out] H the Halfedge
         */
        void move_to_prev_around_facet(Halfedge& H) const {
            geo_debug_assert(halfedge_is_valid(H));
            H.corner = mesh_.facets.prev_corner_around_facet(H.facet, H.corner);
        }
        
        /**
         * \brief Replaces a Halfedge with the next one around the vertex.
         * \param[in,out] H the Halfedge
         * \return true if the move was successful, false otherwise. On borders,
         *  the next halfedge around a vertex may not exist.
         */
        bool move_to_next_around_vertex(Halfedge& H, bool ignore_borders = false) const;

        /**
         * \brief Replaces a Halfedge with the previous one around the vertex.
         * \param[in,out] H the Halfedge
         * \return true if the move was successful, false otherwise. On borders,
         *  the previous halfedge around a vertex may not exist.
         */
        bool move_to_prev_around_vertex(Halfedge& H, bool ignore_borders = false) const;

        /**
         * \brief Replaces a Halfedge with the next one around the border.
         * \details If set_use_facet_region() is set, then
         *  Halfedges incident to two different facet regions are
         *  considered as borders.
         * \param[in,out] H the Halfedge
         */
        void move_to_next_around_border(Halfedge& H) const;

        // why move_to_next_around_border() does not have the same behavior ??
        void custom_move_to_next_around_border(Halfedge& H) const;

        /**
         * \brief Replaces a Halfedge with the previous one around the border.
         * \details If set_use_facet_region() is set, then
         *  Halfedges incident to two different facet regions are
         *  considered as borders.
         * \param[in,out] H the Halfedge
         */
        void move_to_prev_around_border(Halfedge& H) const;
        
        /**
         * \brief Replaces a Halfedge with the opposite one in the
         *  adjacent facet.
         * \param[in,out] H the Halfedge
         * \pre !is_on_border(H)
         */
        void move_to_opposite(Halfedge& H) const;

    private:
        Mesh& mesh_;
        Attribute<index_t> facet_region_;
    };

    /**
     * \brief Displays a Halfedge.
     * \param[out] out the stream where to print the Halfedge
     * \param[in] H the Halfedge
     * \return a reference to the stream \p out
     */
    inline std::ostream& operator<< (
        std::ostream& out, const CustomMeshHalfedges::Halfedge& H
    ) {
        return out << '(' << H.facet << ',' << H.corner << ')';
    }

        /**
         * \brief Gets the origin point of a Halfedge
         * \param[in] M the mesh
         * \param[in] H the Halfedge
         * \return a const reference to the origin of \p H
         */
        inline const vec3& halfedge_vertex_from(
            const Mesh& M, const CustomMeshHalfedges::Halfedge& H
        ) {
            return mesh_vertex(M, M.facet_corners.vertex(H.corner));
        }

        /**
         * \brief Gets the arrow extremity point of a Halfedge
         * \param[in] M the mesh
         * \param[in] H the Halfedge
         * \return a const reference to the arrow extremity of \p H
         */
        inline const vec3& halfedge_vertex_to(
            const Mesh& M, const CustomMeshHalfedges::Halfedge& H
        ) {
            index_t c = M.facets.next_corner_around_facet(H.facet, H.corner);
            return mesh_vertex(M, M.facet_corners.vertex(c));
        }

        /**
         * \brief Gets a 3d vector that connects the origin with the arrow
         *  extremity of a Halfedge.
         * \param[in] M the Mesh
         * \param[in] H the Halfedge
         * \return a 3d vector that connects the origin with the arrow
         *  extremity of \p H
         */
        inline vec3 halfedge_vector(
            const Mesh& M, const CustomMeshHalfedges::Halfedge& H
        ) {
            return halfedge_vertex_to(M, H) - halfedge_vertex_from(M, H);
        }

        /**
         * \brief Gets the length of a Halfedge
         * \param[in] M the Mesh
         * \param[in] H the Halfedge
         * \return the 3d length of \p H
         */
        inline double edge_length(
            const Mesh& M, const CustomMeshHalfedges::Halfedge& H
        ) {
            return length(halfedge_vector(M, H));
        }


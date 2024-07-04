#include <gtest/gtest.h>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <array>
#include <random>

#include <geogram/mesh/mesh.h>                  // for Mesh
#include <geogram/basic/vecg.h>                 // for vec3
#include <geogram/mesh/mesh_io.h>               // for mesh_save()
#include <geogram/basic/command_line.h>         // for CmdLine::initialize()
#include <geogram/basic/command_line_args.h>	// for CmdLine::import_arg_group()
#include <geogram/basic/file_system.h>          // for FileSystem::initialize()
#include <geogram/mesh/mesh_halfedges.h>        // for MeshHalfedges::Halfedge

#include "stats.h"                  // for BasicStats, std_dev()
#include "geometry_halfedges.h"     // for MeshHalfedgesExt
#include "labeling_generators.h"    // for naive_labeling()
#include "geometry.h"               // for comparison between vec3

#define DOUBLE_MAX_ABS_ERROR 10e5

using namespace GEO;

TEST(BasicStats, HardCoded) {
    BasicStats stats;
    std::array<double,10> values = {
        2.3,
        6.7,
        5.1,
        5.0,
        3.9,
        1.4,
        6.4,
        9.2,
        2.0,
        5.1
    };
    FOR(n,10) {
        stats.insert(values[n]);
    }
    EXPECT_EQ(10, stats.count());
    EXPECT_DOUBLE_EQ(47.1, stats.sum());
    EXPECT_DOUBLE_EQ(1.4, stats.min());
    EXPECT_DOUBLE_EQ(9.2, stats.max());
    EXPECT_DOUBLE_EQ(4.71, stats.avg());
    EXPECT_NEAR(5.2129, stats.variance(), DOUBLE_MAX_ABS_ERROR);
    EXPECT_NEAR(2.2831776102616, stats.sd(), DOUBLE_MAX_ABS_ERROR);
}

TEST(BasicStats, Random_1000_sd) { // test the standard deviation computation on 1000 random values
    BasicStats stats;
    std::array<double,1000> values;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution dis(0.0,100.0);
    FOR(n,1000) {
        values[n] = dis(gen);
        stats.insert(values[n]);
    }
    // compare standard deviation computation from all the values, and the iterative computation
    EXPECT_NEAR(std_dev(values.begin(),values.end()), stats.sd(), DOUBLE_MAX_ABS_ERROR);
}

// Halfedge operators tests
// http://google.github.io/googletest/primer.html#same-data-multiple-tests
class HalfedgesTest : public testing::Test {
public:

    HalfedgesTest() : testing::Test(),
                      mesh_halfedges(cube) {}; // link the halfedges API to the cube mesh

protected:
    void SetUp() override {

        // silent Geogram initialization
        // copy of what's inside GEO::initialize(), modified to have the logger in the "minimal" mode
        Environment* env = Environment::instance();
        env->set_value("version", "inaccessible"); // some code after expects the "version" environment variable to exist
        env->set_value("release_date", "inaccessible"); // idem
        env->set_value("SVN revision", "inaccessible"); // idem
        FileSystem::initialize();
        Logger::initialize();
        Logger::instance()->set_minimal(true);
        CmdLine::initialize();
        CmdLine::import_arg_group("sys"); // declares sys:compression_level, needed by mesh_save() for .geogram files
        mesh_io_initialize();

        // Create a simple triangle mesh of a cube

        // create the 8 vertices of the cube
        // Vertex ordering from the Geogram convention for hexahedra, see include/geometry_hexahedra.h
        //
        //                  vertex | X | Y | Z |
        //                  -------|---|---|---|
        //          4-------6   v0 | 0 | 0 | 1 |
        //         /|      /|   v1 | 0 | 0 | 0 |
        //        / |     / |   v2 | 1 | 0 | 1 |
        //       0-------2  |   v3 | 1 | 0 | 0 |
        //       |  5----|--7   v4 | 0 | 1 | 1 |
        // Z     | /     | /    v5 | 0 | 1 | 0 |
        // ^  Y  |/      |/     v6 | 1 | 1 | 1 |
        // | /   1-------3      v7 | 1 | 1 | 0 |
        // |/
        // o----> X
        index_t index_of_first_vertex = cube.vertices.create_vertices(8);
        geo_assert(index_of_first_vertex == 0);
        cube.vertices.point(0) = GEO::vec3(0.0, 0.0, 1.0);
        cube.vertices.point(1) = GEO::vec3(0.0, 0.0, 0.0);
        cube.vertices.point(2) = GEO::vec3(1.0, 0.0, 1.0);
        cube.vertices.point(3) = GEO::vec3(1.0, 0.0, 0.0);
        cube.vertices.point(4) = GEO::vec3(0.0, 1.0, 1.0);
        cube.vertices.point(5) = GEO::vec3(0.0, 1.0, 0.0);
        cube.vertices.point(6) = GEO::vec3(1.0, 1.0, 1.0);
        cube.vertices.point(7) = GEO::vec3(1.0, 1.0, 0.0);

        // create the 6*2 triangle facets
        // inside a given facet, local vertices will be ordered to have an outgoing normal (right hand rule -> counterclockwise)
        // starting with the vertex of smallest index
        //
        //                                | vertex (facet corner) |
        //                                | at each local vertex  |
        //                                |-------|-------|-------|
        //                                |   0   |   1   |   2   |
        //                                |-------|-------|-------|
        // 0-------2    0  0-------2      |       |       |       |
        // |       |    | \  \  f1 |   f0 | 0 (0) | 1 (1) | 3 (2) |
        // | front | => |   \  \   |      |       |       |       |
        // |       |    | f0  \  \ |   f1 | 0 (3) | 3 (4) | 2 (5) |
        // 1-------3    1-------3  3      |       |       |       |
        //                                |       |       |       |
        // 4-------0    4  4-------0      |       |       |       |
        // |       |    | \  \  f3 |   f2 | 1 (6) | 4 (7) | 5 (8) |
        // | left  | => |   \  \   |      |       |       |       |
        // |       |    | f2  \  \ |   f3 | 0 (9) | 4 (10)| 1 (11)|
        // 5-------1    5-------1  1      |       |       |       |
        //                                |       |       |       |
        // 6-------4    6  6-------4      |       |       |       |
        // |       |    | \  \  f5 |   f4 | 5 (12)| 6 (13)| 7 (14)|
        // | back  | => |   \  \   |      |       |       |       |
        // |       |    | f4  \  \ |   f5 | 4 (15)| 6 (16)| 5 (17)|
        // 7-------5    7-------5  5      |       |       |       |
        //                                |       |       |       |
        // 2-------6    2  2-------6      |       |       |       |
        // |       |    | \  \  f7 |   f6 | 2 (18)| 3 (19)| 7 (20)|
        // | right | => |   \  \   |      |       |       |       |
        // |       |    | f6  \  \ |   f7 | 2 (21)| 7 (22)| 6 (23)|
        // 3-------7    3-------7  7      |       |       |       |
        //                                |       |       |       |
        // 4-------6    4  4-------6      |       |       |       |
        // |       |    | \  \  f9 |   f8 | 0 (24)| 2 (25)| 4 (26)|
        // |  top  | => |   \  \   |      |       |       |       |
        // |       |    | f8  \  \ |   f9 | 2 (27)| 6 (28)| 4 (29)|
        // 0-------2    0-------2  2      |       |       |       |
        //                                |       |       |       |
        // 1-------3    1  1-------3      |       |       |       |
        // |       |    | \  \  f11|  f10 | 1 (30)| 5 (31)| 7 (32)|
        // |bottom | => |   \  \   |      |       |       |       |
        // |       |    |f10  \  \ |  f11 | 1 (33)| 7 (34)| 3 (35)|
        // 5-------7    5-------7  7      |       |       |       |
        //
        // About local edges (le) : local edge k is the one between local vertices k and (k+1)%3
        // See https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
        //
        index_t index_of_first_facet = cube.facets.create_triangles(12);
        geo_assert(index_of_first_facet == 0);
        // facet 0
        cube.facets.set_vertex(0,0,0);
        cube.facets.set_vertex(0,1,1);
        cube.facets.set_vertex(0,2,3);
        // facet 1
        cube.facets.set_vertex(1,0,0);
        cube.facets.set_vertex(1,1,3);
        cube.facets.set_vertex(1,2,2);
        // facet 2
        cube.facets.set_vertex(2,0,1);
        cube.facets.set_vertex(2,1,4);
        cube.facets.set_vertex(2,2,5);
        // facet 3
        cube.facets.set_vertex(3,0,0);
        cube.facets.set_vertex(3,1,4);
        cube.facets.set_vertex(3,2,1);
        // facet 4
        cube.facets.set_vertex(4,0,5);
        cube.facets.set_vertex(4,1,6);
        cube.facets.set_vertex(4,2,7);
        // facet 5
        cube.facets.set_vertex(5,0,4);
        cube.facets.set_vertex(5,1,6);
        cube.facets.set_vertex(5,2,5);
        // facet 6
        cube.facets.set_vertex(6,0,2);
        cube.facets.set_vertex(6,1,3);
        cube.facets.set_vertex(6,2,7);
        // facet 7
        cube.facets.set_vertex(7,0,2);
        cube.facets.set_vertex(7,1,7);
        cube.facets.set_vertex(7,2,6);
        // facet 8
        cube.facets.set_vertex(8,0,0);
        cube.facets.set_vertex(8,1,2);
        cube.facets.set_vertex(8,2,4);
        // facet 9
        cube.facets.set_vertex(9,0,2);
        cube.facets.set_vertex(9,1,6);
        cube.facets.set_vertex(9,2,4);
        // facet 10
        cube.facets.set_vertex(10,0,1);
        cube.facets.set_vertex(10,1,5);
        cube.facets.set_vertex(10,2,7);
        // facet 11
        cube.facets.set_vertex(11,0,1);
        cube.facets.set_vertex(11,1,7);
        cube.facets.set_vertex(11,2,3);

        // auto-compute facets adjacency
        cube.facets.connect();

        // initialize the halfedge with the one going from vertex 1 to vertex 3
        halfedge.facet = 0; // the facet at its left is facet 0
        halfedge.corner = 1; // the facet corner of facet 0 that is on the origin vertex (1) is 1
    }

    void print_facet_corners() {
        index_t facet_corner = index_t(-1);
        for(index_t f : cube.facets) { // for each facet
            FOR(lv,3) { // for each local vertex of the current facet
                facet_corner = cube.facets.corner(f,lv);
                fmt::println("At facet {} local vertex {} is facet corner {}, which is on (global) vertex {}. Adjacent facet is {}.",
                    f,
                    lv,
                    facet_corner,
                    cube.facet_corners.vertex(facet_corner),
                    OPTIONAL_TO_STRING(cube.facet_corners.adjacent_facet(facet_corner))
                );
            }
        }
    }

    void compute_normals() {
        normals.resize(cube.facets.nb());
        FOR(f,cube.facets.nb()) {
            normals[f] = normalize(Geom::mesh_facet_normal(cube,f));
        }
    }

    Mesh cube;
    MeshHalfedgesExt mesh_halfedges;
    MeshHalfedges::Halfedge halfedge;
    std::vector<vec3> normals;
};

TEST_F(HalfedgesTest, ExportCubeMesh) {
    EXPECT_TRUE(GEO::mesh_save(cube,"cube.obj"));
    // .geogram export doesn't seem to work properly
    // Or graphite doesn't properly render .geogram files with few elements
}

TEST_F(HalfedgesTest, InitializedHalfedge) {
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveToOpposite) {
    mesh_halfedges.move_to_opposite(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 3 to 1
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),35);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),1);

    mesh_halfedges.move_to_opposite(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveToPrevAroundFacet) {
    mesh_halfedges.move_to_prev_around_facet(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 0 to 1
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),11);

    mesh_halfedges.move_to_prev_around_facet(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 3 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),3);

    mesh_halfedges.move_to_prev_around_facet(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveToNextAroundFacet) {
    mesh_halfedges.move_to_next_around_facet(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 3 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),3);

    mesh_halfedges.move_to_next_around_facet(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 0 to 1
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),11);

    mesh_halfedges.move_to_next_around_facet(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveClockwiseAroundVertex) {
    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 7
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),7);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),34);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),32);

    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 5
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),31);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),8);

    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 4
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),7);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),10);

    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),0);

    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveCounterclockwiseAroundVertex) {
    EXPECT_TRUE(mesh_halfedges.move_counterclockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),0);

    EXPECT_TRUE(mesh_halfedges.move_counterclockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 4
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),7);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),10);

    EXPECT_TRUE(mesh_halfedges.move_counterclockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 5
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),31);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),8);

    EXPECT_TRUE(mesh_halfedges.move_counterclockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 7
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),7);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),34);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),32);

    EXPECT_TRUE(mesh_halfedges.move_counterclockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, ComputeNormals) {
    compute_normals();
    EXPECT_EQ(normals.size(),12);
    // facets 0 and 1 are on the front square, toward -Y
    EXPECT_EQ(normals[0],vec3(0.0,-1.0,0.0));
    EXPECT_EQ(normals[1],vec3(0.0,-1.0,0.0));
    // facets 2 and 3 are on the left square, toward -X
    EXPECT_EQ(normals[2],vec3(-1.0,0.0,0.0));
    EXPECT_EQ(normals[3],vec3(-1.0,0.0,0.0));
    // facets 4 and 5 are on the back square, toward +Y
    EXPECT_EQ(normals[4],vec3(0.0,1.0,0.0));
    EXPECT_EQ(normals[5],vec3(0.0,1.0,0.0));
    // facets 6 and 7 are on the right square, toward +X
    EXPECT_EQ(normals[6],vec3(1.0,0.0,0.0));
    EXPECT_EQ(normals[7],vec3(1.0,0.0,0.0));
    // facets 8 and 9 are on the top square, toward +Z
    EXPECT_EQ(normals[8],vec3(0.0,0.0,1.0));
    EXPECT_EQ(normals[9],vec3(0.0,0.0,1.0));
    // facets 10 and 11 are on the bottom square, toward -Z
    EXPECT_EQ(normals[10],vec3(0.0,0.0,-1.0));
    EXPECT_EQ(normals[11],vec3(0.0,0.0,-1.0));
}

TEST_F(HalfedgesTest, NaiveLabeling) {
    compute_normals();
    naive_labeling(cube,normals,"label");
    Attribute<index_t> label(cube.facets.attributes(), "label"); // retrieve the labeling
    EXPECT_EQ(label.size(),12);
    // facets 0 and 1 are on the front square, toward -Y = 3
    EXPECT_EQ(label[0],3);
    EXPECT_EQ(label[1],3);
    // facets 2 and 3 are on the left square, toward -X = 1
    EXPECT_EQ(label[2],1);
    EXPECT_EQ(label[3],1);
    // facets 4 and 5 are on the back square, toward +Y = 2
    EXPECT_EQ(label[4],2);
    EXPECT_EQ(label[5],2);
    // facets 6 and 7 are on the right square, toward +X = 0
    EXPECT_EQ(label[6],0);
    EXPECT_EQ(label[7],0);
    // facets 8 and 9 are on the top square, toward +Z = 4
    EXPECT_EQ(label[8],4);
    EXPECT_EQ(label[9],4);
    // facets 10 and 11 are on the bottom square, toward -Z = 5
    EXPECT_EQ(label[10],5);
    EXPECT_EQ(label[11],5);
}

TEST_F(HalfedgesTest, MoveToPrevAroundBorder) {
    compute_normals();
    naive_labeling(cube,normals,"label");
    mesh_halfedges.set_use_facet_region("label");

    mesh_halfedges.move_to_prev_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 0 to 1
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),11);

    mesh_halfedges.move_to_prev_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 2 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),8);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),25);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),24);

    mesh_halfedges.move_to_prev_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 3 to 2
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),19);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),18);

    mesh_halfedges.move_to_prev_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveToNextAroundBorder) {
    compute_normals();
    naive_labeling(cube,normals,"label");
    mesh_halfedges.set_use_facet_region("label");

    mesh_halfedges.move_to_next_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 3 to 2
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),19);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),18);

    mesh_halfedges.move_to_next_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 2 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),8);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),25);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),24);

    mesh_halfedges.move_to_next_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 0 to 1
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),11);

    mesh_halfedges.move_to_next_around_border(halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

TEST_F(HalfedgesTest, MoveClockwiseAroundVertexIgnoreBorders) {
    compute_normals();
    naive_labeling(cube,normals,"label");
    mesh_halfedges.set_use_facet_region("label");

    // move should fail because halfedge is on border
    EXPECT_FALSE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);

    // move should success because we ignore borders (2nd argument)
    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge,true));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_FALSE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 7
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),7);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),34);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),32);

    // move should success because the left facet is still on the same region
    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 5
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),31);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),8);

    // move should fail because halfedge is on border
    EXPECT_FALSE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 5
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),5);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),10);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),30);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),31);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),8);

    // move should success because we ignore borders (2nd argument)
    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge,true));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_FALSE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 4
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),4);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),6);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),7);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),10);

    // move should success because the left facet is still on the same region
    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),0);

    // move should fail because halfedge is on border
    EXPECT_FALSE(mesh_halfedges.move_clockwise_around_vertex(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 0
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),0);

    // move should success because we ignore borders (2nd argument)
    EXPECT_TRUE(mesh_halfedges.move_clockwise_around_vertex(halfedge,true));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(halfedge));
    EXPECT_TRUE(mesh_halfedges.halfedge_is_border(halfedge));
    // halfedge going from vertices 1 to 3
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_left(cube,halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,halfedge),35);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
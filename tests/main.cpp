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

#include "basic_stats.h"            // for BasicStats
#include "containers.h"             // for std_dev()
#include "CustomMeshHalfedges.h"    // for CustomMeshHalfedges

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
        // Vertex ordering from the Geogram convention for hexahedra, see include/hex_mesh.h
        //          4-------6   v0 = (0,0,1)
        //         /|      /|   v1 = (0,0,0)
        //        / |     / |   v2 = (1,0,1)
        //       0-------2  |   v3 = (1,0,0)
        //       |  5----|--7   v4 = (0,1,1)
        // Z     | /     | /    v5 = (0,1,0)
        // ^  Y  |/      |/     v6 = (1,1,1)
        // | /   1-------3      v7 = (1,1,0)
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
        // 0-------2    0  0-------2  
        // |       |    | \  \  f1 |  f0 = vertices 0,1,3 = facet corners 0,1,2
        // | front | => |   \  \   |
        // |       |    | f0  \  \ |  f1 = vertices 0,3,2 = facet corners 3,4,5
        // 1-------3    1-------3  3
        //
        // 4-------0    4  4-------0  
        // |       |    | \  \  f3 |  f2 = vertices 1,4,5 = facet corners 6,7,8
        // | left  | => |   \  \   |
        // |       |    | f2  \  \ |  f3 = vertices 0,4,1 = facet corners 9,10,11
        // 5-------1    5-------1  1
        //
        // 6-------4    6  6-------4  
        // |       |    | \  \  f5 |  f4 = vertices 5,6,7 = facet corners 12,13,14
        // | back  | => |   \  \   |
        // |       |    | f4  \  \ |  f5 = vertices 4,6,5 = facet corners 15,16,17
        // 7-------5    7-------5  5
        //
        // 2-------6    2  2-------6  
        // |       |    | \  \  f7 |  f6 = vertices 2,3,7 = facet corners 18,19,20
        // | right | => |   \  \   |
        // |       |    | f6  \  \ |  f7 = vertices 2,7,6 = facet corners 21,22,23
        // 3-------7    3-------7  7
        //
        // 4-------6    4  4-------6  
        // |       |    | \  \  f9 |  f8 = vertices 0,2,4 = facet corners 24,25,26
        // |  top  | => |   \  \   |
        // |       |    | f8  \  \ |  f9 = vertices 2,6,4 = facet corners 27,28,29
        // 0-------2    0-------2  2
        //
        // 1-------3    1  1-------3  
        // |       |    | \  \  f11|  f10 = vertices 1,5,7 = facet corners 30,31,32
        // |bottom | => |   \  \   |
        // |       |    |f10  \  \ |  f11 = vertices 1,7,3 = facet corners 33,34,35
        // 5-------7    5-------7  7
        //
        // About local edges (le) : local edge k is the one between k and (k+1)%3
        // See https://github.com/BrunoLevy/geogram/wiki/Mesh#triangulated-and-polygonal-meshes
        //
        index_t index_of_first_facet = cube.facets.create_triangles(12);
        geo_assert(index_of_first_facet == 0);
        // facet 0
        cube.facets.set_vertex(0,0,0);
        cube.facets.set_vertex(0,1,1);
        cube.facets.set_vertex(0,2,3);
        cube.facets.set_adjacent(0,0,3);
        cube.facets.set_adjacent(0,1,11);
        cube.facets.set_adjacent(0,2,1);
        // facet 1
        cube.facets.set_vertex(1,0,0);
        cube.facets.set_vertex(1,1,3);
        cube.facets.set_vertex(1,2,2);
        cube.facets.set_adjacent(1,0,0);
        cube.facets.set_adjacent(1,1,6);
        cube.facets.set_adjacent(1,2,8);
        // facet 2
        cube.facets.set_vertex(2,0,1);
        cube.facets.set_vertex(2,1,4);
        cube.facets.set_vertex(2,2,5);
        cube.facets.set_adjacent(2,0,3);
        cube.facets.set_adjacent(2,1,5);
        cube.facets.set_adjacent(2,2,10);
        // facet 3
        cube.facets.set_vertex(3,0,0);
        cube.facets.set_vertex(3,1,4);
        cube.facets.set_vertex(3,2,1);
        cube.facets.set_adjacent(3,0,8);
        cube.facets.set_adjacent(3,1,2);
        cube.facets.set_adjacent(3,2,0);
        // facet 4
        cube.facets.set_vertex(4,0,5);
        cube.facets.set_vertex(4,1,6);
        cube.facets.set_vertex(4,2,7);
        cube.facets.set_adjacent(4,0,5);
        cube.facets.set_adjacent(4,1,7);
        cube.facets.set_adjacent(4,2,10);
        // facet 5
        cube.facets.set_vertex(5,0,4);
        cube.facets.set_vertex(5,1,6);
        cube.facets.set_vertex(5,2,5);
        cube.facets.set_adjacent(5,0,9);
        cube.facets.set_adjacent(5,1,4);
        cube.facets.set_adjacent(5,2,2);
        // facet 6
        cube.facets.set_vertex(6,0,2);
        cube.facets.set_vertex(6,1,3);
        cube.facets.set_vertex(6,2,7);
        cube.facets.set_adjacent(6,0,1);
        cube.facets.set_adjacent(6,1,11);
        cube.facets.set_adjacent(6,2,7);
        // facet 7
        cube.facets.set_vertex(7,0,2);
        cube.facets.set_vertex(7,1,7);
        cube.facets.set_vertex(7,2,6);
        cube.facets.set_adjacent(7,0,6);
        cube.facets.set_adjacent(7,1,4);
        cube.facets.set_adjacent(7,2,9);
        // facet 8
        cube.facets.set_vertex(8,0,0);
        cube.facets.set_vertex(8,1,2);
        cube.facets.set_vertex(8,2,4);
        cube.facets.set_adjacent(8,0,1);
        cube.facets.set_adjacent(8,1,9);
        cube.facets.set_adjacent(8,2,3);
        // facet 9
        cube.facets.set_vertex(9,0,2);
        cube.facets.set_vertex(9,1,6);
        cube.facets.set_vertex(9,2,4);
        cube.facets.set_adjacent(9,0,7);
        cube.facets.set_adjacent(9,1,5);
        cube.facets.set_adjacent(9,2,8);
        // facet 10
        cube.facets.set_vertex(10,0,1);
        cube.facets.set_vertex(10,1,5);
        cube.facets.set_vertex(10,2,7);
        cube.facets.set_adjacent(10,0,2);
        cube.facets.set_adjacent(10,1,4);
        cube.facets.set_adjacent(10,2,11);
        // facet 11
        cube.facets.set_vertex(11,0,1);
        cube.facets.set_vertex(11,1,7);
        cube.facets.set_vertex(11,2,3);
        cube.facets.set_adjacent(11,0,10);
        cube.facets.set_adjacent(11,1,6);
        cube.facets.set_adjacent(11,2,0);

        // initialize the halfedge with the one going from vertex 1 to vertex 3
        init_halfedge.facet = 0; // the facet at its left is facet 11
        init_halfedge.corner = 1; // the facet corner of facet 11 that is on the origin vertex (1) is 1
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

    Mesh cube;
    CustomMeshHalfedges mesh_halfedges;
    MeshHalfedges::Halfedge init_halfedge;
};

TEST_F(HalfedgesTest, ExportCubeMesh) {
    EXPECT_TRUE(GEO::mesh_save(cube,"cube.obj"));
    // .geogram export doesn't seem to work properly
    // Or graphite doesn't properly render .geogram files with few elements
}

TEST_F(HalfedgesTest, InitializedHalfedge) {
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(init_halfedge));
    EXPECT_EQ(Geom::halfedge_facet_left(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,init_halfedge),11);
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,init_halfedge),3);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,init_halfedge),33);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,init_halfedge),2);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,init_halfedge),35);
}

TEST_F(HalfedgesTest, MoveToPrevAroundFacet) {
    mesh_halfedges.move_to_prev_around_facet(init_halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(init_halfedge));
    EXPECT_EQ(Geom::halfedge_facet_left(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,init_halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,init_halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,init_halfedge),11);
}

TEST_F(HalfedgesTest, MoveToNextAroundFacet) {
    mesh_halfedges.move_to_next_around_facet(init_halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(init_halfedge));
    EXPECT_EQ(Geom::halfedge_facet_left(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,init_halfedge),3);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,init_halfedge),2);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,init_halfedge),4);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,init_halfedge),3);
}

TEST_F(HalfedgesTest, MoveClockwiseAroundVertex) {
    mesh_halfedges.move_clockwise_around_vertex(init_halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(init_halfedge));
    EXPECT_EQ(Geom::halfedge_facet_left(cube,init_halfedge),11);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,init_halfedge),10);
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,init_halfedge),7);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,init_halfedge),33);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,init_halfedge),30);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,init_halfedge),34);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,init_halfedge),32);
}

TEST_F(HalfedgesTest, MoveCounterclockwiseAroundVertex) {
    mesh_halfedges.move_counterclockwise_around_vertex(init_halfedge);
    EXPECT_TRUE(mesh_halfedges.halfedge_is_valid(init_halfedge));
    EXPECT_EQ(Geom::halfedge_facet_left(cube,init_halfedge),3);
    EXPECT_EQ(Geom::halfedge_facet_right(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_vertex_index_from(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_vertex_index_to(cube,init_halfedge),0);
    EXPECT_EQ(Geom::halfedge_bottom_left_corner(cube,init_halfedge),11);
    EXPECT_EQ(Geom::halfedge_bottom_right_corner(cube,init_halfedge),1);
    EXPECT_EQ(Geom::halfedge_top_left_corner(cube,init_halfedge),9);
    EXPECT_EQ(Geom::halfedge_top_right_corner(cube,init_halfedge),0);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
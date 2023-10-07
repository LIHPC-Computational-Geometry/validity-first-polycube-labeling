#include <gtest/gtest.h>

#include <array>
#include <random>

#include "basic_stats.h"    // for BasicStats
#include "containers.h"     // for std_dev()

#define DOUBLE_MAX_ABS_ERROR 10e5

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

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
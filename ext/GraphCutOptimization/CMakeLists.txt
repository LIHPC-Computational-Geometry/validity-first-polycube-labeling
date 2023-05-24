cmake_minimum_required(VERSION 2.8)

project(gco)

SET(GCO_HEADERS
	src/LinkedBlockList.h
    src/graph.h
    src/block.h
    src/GCoptimization.h
)

SET(GCO_SOURCES
	src/LinkedBlockList.cpp
    src/graph.cpp
    src/maxflow.cpp
    src/GCoptimization.cpp
)

add_library(gco STATIC ${GCO_HEADERS} ${GCO_SOURCES})
target_include_directories(gco PUBLIC src)
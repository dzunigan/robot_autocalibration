project(sync CXX)

## CMake config
cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

## Set CXX optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wpedantic -Wall -Wextra")

## Set build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
else(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
  set(MRPT_DONT_USE_DBG_LIBS ON)
endif(CMAKE_BUILD_TYPE STREQUAL "Release")
message(STATUS "Flags: " ${CMAKE_CXX_FLAGS})

## Dependencies
set(GFLAGS_PREFER_EXPORTED_GFLAGS_CMAKE_CONFIGURATION FALSE CACHE BOOL "Prefer exported GFlags?")

find_package(Boost REQUIRED filesystem system)
#find_package(Ceres REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(MRPT REQUIRED base opengl gui)
#find_package(OpenCV REQUIRED)
find_package(Gflags REQUIRED)
find_package(Glog REQUIRED)
#find_package(SQLite3 REQUIRED)

#find_package(Git)

#include(GenerateVersionDefinitions)

# ${OpenCV_INCLUDE_DIRS} ${SQLITE3_INCLUDE_DIRS}
include_directories(include ${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${GFLAGS_INCLUDE_DIR} ${GLOG_INCLUDE_DIRS})

## Header files
file(GLOB_RECURSE HEADER_FILES include/*.hpp include/*.h)
add_custom_target(sync_header_files SOURCES ${HEADER_FILES})

## Link executables
add_executable(sync src/sync.cpp)
target_link_libraries(sync ${Boost_LIBRARIES} ${GLOG_LIBRARIES} ${GFLAGS_LIBRARY})


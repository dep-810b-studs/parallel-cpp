cmake_minimum_required(VERSION 3.8)
project(work_01)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_PREFIX_PATH TBB)
set(TDD_DIR libs/tbb)

add_executable(work_01 main.cpp vector.cpp vector.h struct.h optimizator.h optimizator.cpp)
 
option(JSON_BuildTests OFF)
option(JSON_Install OFF)

add_subdirectory(libs/json)
#add_subdirectory(libs/tbb/cmake)


#target_link_libraries(work_01 PRIVATE nlohmann_json::nlohmann_json)

include(./libs/tbb/cmake/TBBGet.cmake)
include(./libs/tbb/cmake/TBBBuild.cmake)
#tbb_get(TBB_ROOT tbb_root SOURCE_CODE)
tbb_build(TBB_ROOT ${tbb_root} CONFIG_DIR TBB_DIR)
find_package(TBB REQUIRED)

target_link_libraries(work_01 ${TBB_IMPORTED_TARGETS} nlohmann_json::nlohmann_json)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")



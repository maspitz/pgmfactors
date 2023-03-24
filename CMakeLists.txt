cmake_minimum_required(VERSION 3.16)

project(
    pgmfactors
    VERSION 0.1.0
    DESCRIPTION "Probability Factors using the XTensor library"
    HOMEPAGE_URL "https://github.com/maspitz/pgmfactors"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

include(CTest)

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)


# ---- Declare library ----

add_library(pgmfactors_pgmfactors OBJECT source/pgmfactors.cpp)
target_include_directories(pgmfactors_pgmfactors PUBLIC include /usr/include)
target_link_libraries(pgmfactors_pgmfactors xtensor)

set_target_properties(
    pgmfactors_pgmfactors PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME pgmfactors
    OUTPUT_NAME pgmfactors
)

add_executable(pgmexample apps/pgmexample.cpp)
target_link_libraries(pgmexample xtensor pgmfactors_pgmfactors)

enable_testing()
find_package(Catch2 REQUIRED)
include(Catch)

add_executable(pgmfactors_test test/pgmfactors_test.cpp)
target_link_libraries(pgmfactors_test PRIVATE pgmfactors_pgmfactors)
catch_discover_tests(pgmfactors_test)


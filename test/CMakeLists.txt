cmake_minimum_required(VERSION 3.14)

project(pgmfactorsTests LANGUAGES CXX)

include(../cmake/project-is-top-level.cmake)
include(../cmake/folders.cmake)

# ---- Dependencies ----

if(PROJECT_IS_TOP_LEVEL)
  find_package(pgmfactors REQUIRED)
  enable_testing()
endif()

find_package(Catch2 REQUIRED)
include(Catch)

# ---- Tests ----

add_executable(pgmfactors_test source/pgmfactors_test.cpp)
target_link_libraries(
    pgmfactors_test PRIVATE
    pgmfactors::pgmfactors
    Catch2::Catch2WithMain
)
target_compile_features(pgmfactors_test PRIVATE cxx_std_17)

catch_discover_tests(pgmfactors_test)

# ---- End-of-file commands ----

add_folders(Test)

cmake_minimum_required(VERSION 3.16)

project(
    pgm
    VERSION 0.1.0
    DESCRIPTION "Probability Factors and Graphical Models using the XTensor library"
    HOMEPAGE_URL "https://github.com/maspitz/pgm"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)

include(CTest)

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)


# ---- Declare library ----

add_library(pgm_pgm OBJECT source/factor.cpp)
target_include_directories(pgm_pgm PUBLIC include /usr/include)
target_link_libraries(pgm_pgm xtensor)

set_target_properties(
    pgm_pgm PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME pgm
    OUTPUT_NAME pgm
)

add_executable(pgmexample apps/pgmexample.cpp)
target_link_libraries(pgmexample xtensor pgm_pgm)

enable_testing()

Include(FetchContent)

FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v3.3.2
)

FetchContent_MakeAvailable(Catch2)
include(Catch)

add_executable(pgmtest test/test_factor.cpp)
target_link_libraries(pgmtest PRIVATE pgm_pgm Catch2::Catch2WithMain)
catch_discover_tests(pgmtest)


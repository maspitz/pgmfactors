cmake_minimum_required(VERSION 3.14)

include(cmake/prelude.cmake)

project(
    pgmfactors
    VERSION 0.1.0
    DESCRIPTION "Short description"
    HOMEPAGE_URL "https://example.com/"
    LANGUAGES CXX
)

include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)

# ---- Declare library ----

add_library(
    pgmfactors_pgmfactors
    source/pgmfactors.cpp
)
add_library(pgmfactors::pgmfactors ALIAS pgmfactors_pgmfactors)

include(GenerateExportHeader)
generate_export_header(
    pgmfactors_pgmfactors
    BASE_NAME pgmfactors
    EXPORT_FILE_NAME export/pgmfactors/pgmfactors_export.hpp
    CUSTOM_CONTENT_FROM_VARIABLE pragma_suppress_c4251
)

if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(pgmfactors_pgmfactors PUBLIC PGMFACTORS_STATIC_DEFINE)
endif()

set_target_properties(
    pgmfactors_pgmfactors PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION_MAJOR}"
    EXPORT_NAME pgmfactors
    OUTPUT_NAME pgmfactors
)

target_include_directories(
    pgmfactors_pgmfactors ${warning_guard}
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_include_directories(
    pgmfactors_pgmfactors SYSTEM
    PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/export>"
)

target_compile_features(pgmfactors_pgmfactors PUBLIC cxx_std_17)

find_package(fmt REQUIRED)
target_link_libraries(pgmfactors_pgmfactors PRIVATE fmt::fmt)

# ---- Install rules ----

if(NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/install-rules.cmake)
endif()

# ---- Developer mode ----

if(NOT pgmfactors_DEVELOPER_MODE)
  return()
elseif(NOT PROJECT_IS_TOP_LEVEL)
  message(
      AUTHOR_WARNING
      "Developer mode is intended for developers of pgmfactors"
  )
endif()

include(cmake/dev-mode.cmake)

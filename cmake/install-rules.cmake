if(PROJECT_IS_TOP_LEVEL)
  set(
      CMAKE_INSTALL_INCLUDEDIR "include/pgmfactors-${PROJECT_VERSION}"
      CACHE PATH ""
  )
endif()

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# find_package(<package>) call for consumers to find this project
set(package pgmfactors)

install(
    DIRECTORY
    include/
    "${PROJECT_BINARY_DIR}/export/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    COMPONENT pgmfactors_Development
)

install(
    TARGETS pgmfactors_pgmfactors
    EXPORT pgmfactorsTargets
    RUNTIME #
    COMPONENT pgmfactors_Runtime
    LIBRARY #
    COMPONENT pgmfactors_Runtime
    NAMELINK_COMPONENT pgmfactors_Development
    ARCHIVE #
    COMPONENT pgmfactors_Development
    INCLUDES #
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

write_basic_package_version_file(
    "${package}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

# Allow package maintainers to freely override the path for the configs
set(
    pgmfactors_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${package}"
    CACHE PATH "CMake package config location relative to the install prefix"
)
mark_as_advanced(pgmfactors_INSTALL_CMAKEDIR)

install(
    FILES cmake/install-config.cmake
    DESTINATION "${pgmfactors_INSTALL_CMAKEDIR}"
    RENAME "${package}Config.cmake"
    COMPONENT pgmfactors_Development
)

install(
    FILES "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
    DESTINATION "${pgmfactors_INSTALL_CMAKEDIR}"
    COMPONENT pgmfactors_Development
)

install(
    EXPORT pgmfactorsTargets
    NAMESPACE pgmfactors::
    DESTINATION "${pgmfactors_INSTALL_CMAKEDIR}"
    COMPONENT pgmfactors_Development
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()

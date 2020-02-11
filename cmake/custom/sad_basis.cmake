set(SAD_BASIS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/share/sad_basis CACHE STRING "Path to SAD basis and density files")
set(SAD_BASIS_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/sad_basis)
install(
  DIRECTORY ${SAD_BASIS_DIR}
  DESTINATION share/${PROJECT_NAME}/sad_basis
  )


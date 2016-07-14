option(MRCPP_INSTALL_FILTERS "Install MW filters" OFF)

if(MRCPP_INSTALL_FILTERS)
    install(DIRECTORY ${MW_FILTER_DIR}
        DESTINATION share
        )
    set(MW_FILTER_DIR ${CMAKE_INSTALL_PREFIX}/share/MultiwaveletFilters)
endif()

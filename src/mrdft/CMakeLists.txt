target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Factory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MRDFT.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Functional.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xc_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/xclib.cpp
    )

if(ENABLE_XCFUN)
    target_sources(mrchem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/xclib_XCFun.cpp
        )
endif()

if(ENABLE_LIBXC)
    target_sources(mrchem PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/xclib_Libxc.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/xc_func_alias.cpp
        )
endif()

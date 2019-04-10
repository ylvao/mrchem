find_package(nlohmann_json 3.5.0 CONFIG QUIET)

if(TARGET nlohmann_json::nlohmann_json)
  get_property(_loc TARGET nlohmann_json::nlohmann_json PROPERTY LOCATION)
  message(STATUS "Found nlohmann_json: ${nlohmann_json_INCLUDE_DIR} (found version ${nlohmann_json_VERSION})")
else()
  message(STATUS "Suitable nlohmann_json could not be located: downloading and building nlohmann_json instead.")
  include(FetchContent)                                                          
                                                                                 
  FetchContent_Populate(nlohmann_json_sources
    QUIET
    GIT_REPOSITORY
      https://github.com/nlohmann/json
    GIT_TAG
      v3.5.0
    LOG_DOWNLOAD 1
    LOG_UPDATE 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
    )                                                                            
    set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
                                                                                 
  add_subdirectory(                                                              
    ${nlohmann_json_sources_SOURCE_DIR}                                                  
    ${nlohmann_json_sources_BINARY_DIR}                                                  
    )
endif()

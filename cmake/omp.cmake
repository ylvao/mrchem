option(ENABLE_OPENMP "Enable OpenMP parallelization" OFF)

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(HAVE_OPENMP TRUE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

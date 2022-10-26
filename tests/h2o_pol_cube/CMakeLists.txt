if(ENABLE_MPI)
    set(_h2o_pol_cube_launcher "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1")
endif()

add_integration_test(
  NAME "H2O_polarizability_CUBE"
  LABELS "H2O_polarizability_CUBE;polarizability;mrchem;h2o_pol_cube"
  COST 100
  LAUNCH_AGENT ${_h2o_pol_cube_launcher}
  INITIAL_GUESS ${CMAKE_CURRENT_LIST_DIR}/initial_guess
  )

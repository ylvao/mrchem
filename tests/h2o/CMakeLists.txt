set(_h2o_launcher " ")

if(ENABLE_MPI)
    set(_h2o_launcher "mpirun -np 5")
endif()

add_integration_test(
  NAME h2o
  LABELS "mrchem;h2o"
  COST 100
  LAUNCH_AGENT ${_h2o_launcher}
  )

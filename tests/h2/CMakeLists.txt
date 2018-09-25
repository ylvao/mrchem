set(_h2_launcher " ")

if(ENABLE_MPI)
    set(_h2_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME h2
  LABELS "mrchem;h2"
  COST 100
  LAUNCH_AGENT ${_h2_launcher}
  )

set(_li_launcher " ")

if(ENABLE_MPI)
    set(_li_launcher "mpirun -np 2")
endif()

add_integration_test(
  NAME li
  LABELS "mrchem;li"
  COST 100
  LAUNCH_AGENT ${_li_launcher}
  )

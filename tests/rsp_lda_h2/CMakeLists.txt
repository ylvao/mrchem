set(_rsp_lda_h2_launcher " ")

if(ENABLE_MPI)
    set(_rsp_lda_h2_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME rsp_lda_h2
  LABELS "mrchem;rsp_lda_h2"
  COST 100
  LAUNCH_AGENT ${_rsp_lda_h2_launcher}
  )

set(_li_rsp_lda_launcher " ")

if(ENABLE_MPI)
    set(_li_rsp_lda_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME li_rsp_lda
  LABELS "mrchem;li_rsp_lda"
  COST 100
  LAUNCH_AGENT ${_li_rsp_lda_launcher}
  )

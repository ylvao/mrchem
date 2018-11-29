set(_li_fd_lda_launcher " ")

if(ENABLE_MPI)
    set(_li_fd_lda_launcher "mpirun -np 1")
endif()

add_integration_test(
  NAME li_fd_lda
  LABELS "mrchem;li_fd_lda"
  COST 100
  LAUNCH_AGENT ${_li_fd_lda_launcher}
  )

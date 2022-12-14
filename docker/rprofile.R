# This suppresses unwanted BLAS and 
RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

# Set the custom library path constantly
.libPaths(c("/r_user_lib/", .libPaths()))


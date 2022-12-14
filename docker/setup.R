#---------------------------------------------------------
# Install some basic packages during the container build
#---------------------------------------------------------

options(repos=c(CRAN="https://cloud.r-project.org"))

install.packages("pak")

pkgs <- c("renv", "RhpcBLASctl", "rmdformats", "bookdown", "markdown", "rmarkdown", "remotes")
pak::pkg_install(pkgs)
pak::cache_clean()

renv::snapshot(packages=rownames(installed.packages()), lockfile="/r_user_lib/lockfile_initial.json", prompt=FALSE)
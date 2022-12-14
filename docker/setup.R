#---------------------------------------------------------
# Install some basic packages during the container build
#---------------------------------------------------------

# pak and renv for package management and RhpcBLASctl to disable
# unwanted BLAS/OMP multithreading as this messes with CPU limits
pkgs <- c("pak", "renv", "RhpcBLASctl", "rmdformats")

install.packages(pkgs, repos=c(CRAN='https://cloud.r-project.org'))

# Core packages to run/render Rmarkdown inside RStudio and remotes
# for installing packages from GitHub
pkgs2 <- c("bookdown", "markdown", "rmarkdown", "remotes")

pak::pkg_install(pkgs2, upgrade=FALSE, ask=FALSE)

# Take a snapshot of this initial library to document versions
renv::snapshot(packages=rownames(installed.packages()),
               lockfile="/r_user_lib/lockfile_initial.json",
               prompt=FALSE)





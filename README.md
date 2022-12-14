
# Seasonal Project

This repository collects analysis code related to the circadian project, find it 
at [GitHub](https://github.com/ATpoint/seasonal_project).

## Data

All NGS data are stored on our WWU Cloud share (`e_expat/expatstorage/seasonal_project/`) which also contains all
documentation, command lines and software versions used for the preprocessing.

Processed data as required the downstream analysis code were copied to the expatspace in folder `dir_processedData/`.

## Docker image

We run all R-based analysis from inside a Docker container that provides a RStudio-Server
instance and comes with all dependencies to compile any Bioconductor packages.
Build it based on the official Bioconductor image:

```bash
cd dir_docker
docker build -t atpoint/seasonal:0.9a .
```

Then start up the container, change visibility and settings in RStudio to deactivate any auto-saves/auto-loads,
install the packages in `R/packages.R` and commit it to `atpoint/seasonal:x.x`. Then push to DockerHub to have
it globally available.

The current package versions are tracked in `R/lockfile_current.json`.

```bash
docker run -d -p 8787:8787 \
  -v "/home/atpoint/expatspace/projects/seasonal_project/":/projectdir \
  -e PASSWORD=sudo_drink_tequila \
  -e ROOT=TRUE \
  atpoint/seasonal:0.9
```

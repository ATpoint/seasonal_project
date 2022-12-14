#-----------------------------------------------------------------------------
# Install packages with pak
#-----------------------------------------------------------------------------

#/ Run this manually to get packages going
if(1 > 2){
  
  #/ base stock of packages
  to_install <- c("ashr",
                  "batchelor",
                  "BiocParallel",
                  "BiocStyle",
                  "biomaRt",
                  "Biostrings",
                  "circacompare",
                  "clusterProfiler",
                  "ComplexHeatmap",
                  "colorspace",
                  "cowplot",
                  "data.table",
                  "DESeq2",
                  "DropletUtils",
                  "edgeR",
                  "fgsea",
                  "fmsb",
                  "genomation",
                  "GenomicRanges",
                  "ggforce",
                  "ggplot2",
                  "ggpointdensity",
                  "ggpubr",
                  "ggrepel",
                  "ggridges",
                  "ggsignif",
                  "ggthemes",
                  "ggupset",
                  "glmGamPoi",
                  "goseq",
                  "gprofiler2",
                  "grid",
                  "GSVA",
                  "hexbin",
                  "limma",
                  "magick",
                  "magrittr",
                  "Matrix",
                  "matrixStats",
                  "openxlsx",
                  "org.Mm.eg.db",
                  "patchwork",
                  "PCAtools",
                  "ReactomeContentService4R",
                  "reshape2",
                  "RobustRankAggreg",
                  "rtracklayer",
                  "RUVSeq",
                  "S4Vectors",
                  "scater",
                  "scDblFinder",
                  "scran",
                  "scuttle",
                  "SingleCellExperiment",
                  "SingleR",
                  "shiny",
                  "slingshot",
                  "SummarizedExperiment",
                  "sva",
                  "tidyverse",
                  "tximport",
                  "carmonalab/UCell",
                  "UpSetR",
                  "uwot",
                  "viridis",
                  "xfun",
                  "writexl",
                  #
                  "atpoint/CreateGeneSignatures",
                  "atpoint/vizzy")
  
  options(repos=c(CRAN="https://cloud.r-project.org"))
  
  pak::pkg_install(to_install)
  pak::pak_cleanup()
  
  # lockfile with renv for reproduction outside the container
  renv::snapshot(packages=rownames(installed.packages()),
                 lockfile="/r_user_lib/lockfile_current.json",
                 prompt=FALSE)
  
  file.copy("/r_user_lib/lockfile_current.json", "/projectdir/seasonalProject/R/lockfile_current.json", 
            overwrite=TRUE)
  
}

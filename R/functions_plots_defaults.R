# Turn off BLAS multithreading
if("RhpcBLASctl" %in% rownames(installed.packages())){
  
  RhpcBLASctl::omp_set_num_threads(1)
  RhpcBLASctl::blas_set_num_threads(1)
  
}

if("data.table" %in% rownames(installed.packages()))
  invisible(data.table(setDTthreads(threads=1)))

#------------------------------------------
#
# Functions
#
#------------------------------------------

#' Oscillation detection via cosinor regression using limma
#' 
#' This function accepts a SummarizedExperiment with the first assay to be used for cosinor regression
#' via limma. By default limma-trend is used so the assay should be logcpm-like. Turning off the trend
#' option allows to use any prenormalized log2-ish data. Data are expected to be normalized, filtered and log2 scale.
#' Optionally, with use_weights to TRUE, the function uses arrayWeights to downweight outlier samples based on the 
#' linear model. The SummarizedExperiment must contain a numeric se$time column with the timepoint information.
#' 
#' @param a SummarizedExperiment with the first assay to be used in limma and a numeric time column
#' @param period numeric, the period of time to fit the cosisor to, default 24
#' @param trend logical, whether to use limma-trend
#' @param robust logical, whether to run eBayes as robust
#' @param use.weights logical, whether to use arrayWeights to downweight outlier samples
#' 
#' @details Reference for sample weights: https://academic.oup.com/nar/article/43/15/e97/2414275 and
#' https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-261
#' 
oscillation_call <- function(se, use.assay="logcpm", period=24, trend=TRUE, robust=FALSE, use.weights=FALSE, 
                             return.weights=FALSE){
  
  # Checks
  if(!is(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment")
  
  if(!use.assay %in% assayNames(se))
    stop("use.assay not in se")
  
  if(is.null(se$time) | !is.numeric(se$time)) stop("se$time must be present and numeric")
  
  if(length(unique(se$time)) < 4) warning("Fewer than four timepoints in se$time")
  
  if(!is.numeric(period)) stop("period must be numeric")
  
  # Setup the design
  colData(se) <- base::droplevels.data.frame(colData(se))
  
  exp_design <- base::cbind(data.frame(colData(se)),
                            cosphase = cos(2 * pi * se$time / period),
                            sinphase = sin(2 * pi * se$time / period))
    
  design <- model.matrix(~sinphase + cosphase, exp_design)
    
  # Use a EList so weights are used if not NULL automatically
  lcpm <- assay(se, use.assay)
  aw <- if(use.weights) limma::arrayWeights(lcpm, design) else NULL
  v <- new("EList", list(targets=data.frame(colData(se)), E=lcpm, weights=aw))
   
  # Get stats
  fit <- limma::lmFit(object=v, design=design)
  fit <- limma::eBayes(fit, robust=robust, trend=trend)
  tt  <- limma::topTable(fit=fit, coef=2:3, number=Inf, sort.by="none")
  tt  <- cbind(data.frame(gene=rownames(tt)), tt)
  
  # Peak-to-peak amplitude, acrophase in "time" units
  tt$amplitude <- 2*sqrt(base::rowSums(tt[,c("sinphase", "cosphase")]^2))
  #tt$acrophase <- (base::atan2(tt$sinphase, tt$cosphase)) %% (2*pi)
  tt$acrophase <- (base::atan2(tt$sinphase, tt$cosphase) / 2 / pi * period + period) %% period
  
  tt$intercept <- fit$coefficients[,"(Intercept)"]
  tt <- tt[,c("gene", "adj.P.Val", "amplitude", "acrophase", "AveExpr", "P.Value", "F", "intercept", "sinphase", "cosphase")]
  colnames(tt)[2] <- "fdr"
  colnames(tt)[6] <- "pvalue"
  
  if(!is.null(v$weights)) names(v$weights) <- colnames(v$E)
  
  slist <- SimpleList(results=tt, weights=v$weights)
  
  if(return.weights) return(slist) else slist$results
  
}

#' Classify oscillation between conditions using model selection
#' 
#' Inspired by the compareRhythms package this function classifies rhythmic genes into the five categories
#' gain, loss, same, change and arrhy contrasting two conditions. For this it fits five different models to the data
#' in a gene-wise fashion and then uses the BIC to select the best one, returning the category and a conditional
#' probability that the user can use for filtering purposes. The ref-level of levels(se$condition) is used as reference,
#' so category gain means more rhythmicity in that level compared to the other level.
oscillation_classify <- function(se, use.assay="logcpm", period=24, use.weights=FALSE, trend=TRUE, robust=FALSE){

  # Checks
  if(!is(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment")

  if(is.null(se$time) | !is.numeric(se$time)) stop("se$time must be present and numeric")
  
  if(length(unique(se$time)) < 4) warning("Fewer than four timepoints in se$time")
  
  if(!is.numeric(period)) stop("period must be numeric")

  if(!is(se$condition, "factor")) stop("se$condition must be a factor")

  lcond <- length(levels(se$condition))
  if(lcond !=2) stop("se$condition must have exactly two levels")
  
  # Setup the design
  colData(se) <- base::droplevels.data.frame(colData(se))
  lvls <- levels(se$condition)
  lcpm <- SummarizedExperiment::assay(se, use.assay)
  
  exp_design <- base::cbind(data.frame(colData(se)),
                            cosphase = cos(2 * pi * se$time / period),
                            sinphase = sin(2 * pi * se$time / period))
    
  # The full or "change" model so both conditions are rhythmic but parameters change.
  design_change <- stats::model.matrix(~condition + condition:sinphase + condition:cosphase,
                                       data=exp_design)
    
  # The "same" model, no difference, meaning no interaction terms
  design_same <- stats::model.matrix(~condition + sinphase + cosphase,
                                       data=exp_design)
    
  # The no rhythmicity (noR) model so without the rhythmic parameters
  design_noR <- stats::model.matrix(~condition, data=exp_design)
    
  # The loss and gain models so rhythmic parameters only for one group or the other
  design_loss <- design_change[,!grepl(paste0(lvls[2], ":"), colnames(design_change))]
  design_gain <- design_change[,!grepl(paste0(lvls[1], ":"), colnames(design_change))]
    
  # Combine into a list for limma::selectModel
  design_list <- list(arrhy=design_noR,
                      loss=design_loss,
                      gain=design_gain,
                      change=design_change,
                      same=design_same)
  
  # Loop through the designs, optionally using design-specific weights
  v <- new("EList", list(targets=data.frame(colData(se)), E=lcpm))
  
  sms <- sapply(names(design_list), function(ds){
    
    design_current <- design_list[ds]
    aw <- if(use.weights) limma::arrayWeights(object=lcpm, design=design_current[[1]])
    vx <- v
    vx$weights <- aw
    sm <- limma::selectModel(y=vx, designlist=design_current, criterion="bic")

  }, simplify=FALSE)
  
  # Get the preferred model
  models <- names(design_list)
  nmodels <- length(design_list)
  ic <- lapply(sms, function(x) x$IC) %>% do.call(cbind, .)
  pref <- factor(apply(ic, 1, which.min), levels=1:nmodels, labels=models)        

  # Calculate the conditional probability as in compareRhythms
  prob <- base::apply(ic, 1, function(x) base::max(exp(-0.5*x)/sum(exp(-0.5*x))))
  
  summarized <- data.frame(gene=rownames(v),
                           model_category=pref,
                           model_probability=prob)
    
  summarized$model_probability[is.nan(summarized$model_probability)] <- NA
  
  # Run a differential analysis based on the full (change) model with limma
  awx <- if(use.weights) limma::arrayWeights(v$E, design=design_change) else NULL
  fitx <- limma::lmFit(v$E, design=design_change, weights=awx)  
  fitx <- limma::eBayes(fitx, trend=trend, robust=robust)
  gx   <- grep(":sin|:cos", colnames(design_change))
  ttx  <- limma::topTable(fitx, coef=gx, number=Inf)
  ttx  <- ttx[,5:ncol(ttx)]
  ttx  <- data.frame(gene=rownames(ttx), ttx)
  colnames(ttx) <- c("gene", "AveExpr", "F", "pvalue", "fdr")
  
  r <- dplyr::left_join(x=summarized, y=ttx, by="gene")
  
  return(r)
    
}

#' Save plots
#' 
#' Save plots stored in a list as pdf and png
#' 
#' @param plotlist a named list with ggplot objects or compatible grobs
#' @param outdir output directory name
#' @param prefix a prefix to give the plot names
#' @param overwrite logical, overwrite if plot already exists
#' @param width width of plot
#' @param height height of plot
#' @param res resolution for the png
#' @param units units for the png 
#' @param verbose logical whether to tell which plot it stored to which location
#' @param create_outdir logical, guess what
#' 
save_plots <- function(plotlist, outdir, prefix="",
                       overwrite=TRUE, width=7, height=7, res=300, units="in",
                       verbose=TRUE, create_outdir=TRUE, skip_png=FALSE){
  
  if(length(plotlist)==0){
    message("[Info] save_plots did not do anything because the input list is empty!")
    return(NULL)
  }
  
  if(is.null(names(plotlist))) 
    stop("plotlist has no names() -- names(plotlist) will be the filenames!")
  
  if(!dir.exists(outdir))
    if(create_outdir){
      message("Creating outdir")
      dir.create(outdir, recursive=TRUE)
    } else stop("Outdir does not exist and create_outdir=FALSE")
  
  for(p in names(plotlist)){
    
    tmp.plot    <- plotlist[[p]]
    tmp.name    <- paste0(outdir, "/", prefix, p)
    tmp.message <- gsub("//|///", "/", paste("Saving", p, "to", paste0(tmp.name, ".png/.pdf")))
    
    if(verbose) tmp.message
    
    # pdf:
    save_pdf <- paste0(tmp.name, ".pdf")
    if(!file.exists(save_pdf) | (file.exists(save_pdf) & overwrite)){
      pdf(paste0(tmp.name, ".pdf"), width=width, height=height)
      print(tmp.plot); dev.off()
    } else message(paste(save_pdf, "exists and overwrite is set to FALSE"))
    
    # png:
    if(!skip_png){
    save_png <- paste0(tmp.name, ".png")
    if(!file.exists(save_png) | (file.exists(save_png) & overwrite)){
      png(paste0(tmp.name, ".png"), width=width, height=height, res=res, units=units)
      print(tmp.plot); dev.off()
    } else message(paste(save_png, "exists and overwrite is set to FALSE"))
    }
    
    rm(tmp.plot, tmp.name, tmp.message)      
    
  }
  
}

#' Function to save content of a list of table-like objects to disk as txt and xlsx
#'
#' @param l named list with elements to save
#' @param outdir path to output dir
#' @param create_outdir logical, guess what
#' @param overwrite logical, guess what
#' @param compress logical, whether to gzip-compress the text file and save as ".gz"
#' @param sep a valid field separator for the text file
#'  
save_lists <- function(l, outdir, create_outdir=TRUE, overwrite=TRUE, compress=FALSE, sep="\t"){
  
  suppressMessages({
    require(data.table)
    require(openxlsx)
  })
  
  if(!is(l, "list")) stop("l must be a list")
  
  if(length(l)==0){
    message("save_lists did not do anything because l was empty")
    return(NULL)
  }
  
  if(is.null(names(l))) stop("l must be named")
  
  if(!dir.exists(outdir))
    if(create_outdir){
      message("Creating outdir")
      dir.create(outdir, recursive=TRUE)
    } else stop("Outdir does not exist and create_outdir=FALSE")
  
  howtocompress <- if(compress) "gz" else "none"
  
  invisible(sapply(names(l), function(x){
    
    suffix <- if(compress) ".gz" else ""
    
    fl <- paste0(outdir, "/", x, ".txt", suffix)
    
    owr <- if(file.exists(fl) & overwrite) TRUE else FALSE
    owr <- if(!file.exists(fl)) TRUE else owr
    
    if(owr){
      
      data.table::fwrite(x=l[[x]], file=fl,
                         row.names=FALSE, col.names=TRUE, sep=sep, quote=FALSE,
                         compress=howtocompress, verbose=FALSE)
      
      openxlsx::write.xlsx(x=l[[x]], file=paste0(outdir, "/", x, ".xlsx"), overwrite=TRUE)
      
    } else message(paste(fl, "already exists -- skipping"))
    
    NULL
    
  }))
  
}

#' Function to define PCA outliers based on 3xMAD from the median based 
#' on the PC coordinates

#' @param x the content of the prcomp(logcounts)$x output
#' @param n_pcs number of PCs to use
#' @param threshold numeric, the threshold in median(x)+/-threshold(mad(x))
#' 
#' @details
#' Function will compute the median coordinate of each PC and then use threshold*MAD
#' to define outliers. Cell is an outlier of beyond the threshold in any of n_pcs
#' 
#' @examples
#' x <- data.frame(PC1=c(4,5,7,2,4,3,6,100, 2), PC2=c(3,4,5,3,1,2,3,4, 200), row.names=paste0("cell",1:9))
#' pca_outliers(x)
pca_outliers <- function(x, n_pcs=1:2, threshold=3){
  
  pc <- x[,n_pcs,drop=FALSE]
  
  median <- apply(pc, 2, stats::median)
  mad    <- apply(pc, 2, stats::mad)
  
  lapply(n_pcs, function(z){
    
    ol <- 
      x[,z,drop=TRUE] > (median[z]+threshold*mad[z]) |
      x[,z,drop=TRUE] < (median[z]-threshold*mad[z])  
    names(ol[ol])
    
  }) -> ol
  
  outliers <- unlist(ol)
  if(length(outliers)==0) return(NULL) else return(outliers)
  
}

#' Trim a matrix of df of counts by quantiles
#' 
#' @param counts a count matrix or data.frame
#' @param lower the lower quantile, e.g. .01
#' @param upper the upper quantile, e.g. .99
#' 
scale_by_quantile <- function (counts, lower = 0, upper = 1){
  
  if (class(counts)[1] != "matrix" & class(counts)[1] != "data.frame") {
    stop("counts must be a matrix or data.frame")
  }
  if (lower == 0 & upper == 1) 
    return(counts)
  if (class(counts)[1] == "data.frame") {
    cnames <- colnames(counts)
    counts <- as.matrix(counts)
    colnames(counts) <- cnames
  }
  if (upper < 1) {
    qt.upper <- as.numeric(quantile(counts, upper, na.rm = TRUE))
    counts[counts > qt.upper] <- qt.upper
  }
  if (lower > 0) {
    qt.lower <- as.numeric(quantile(counts, lower, na.rm = TRUE))
    counts[counts < qt.lower] <- qt.lower
  }
  return(counts)
  
}

#' Save components of a SummarizedExperiment or SingleCellExperiment to disk
#' 
#' Components are the assays, col- and rowData
#' 
#' @param x the SE
#' @param outdir outdir name
#' @param a basename to save the elements
#' 
save_se <- function(x, outdir, outname, create_outdir=TRUE){
  
  if(!is(x, "SummarizedExperiment") & !is(x, "SingleCellExperiment"))
    stop("x must be a SE or SCE")
  
  asynames <- assayNames(x)
  if(is.null(asynames)) stop("Nothing to be saved here -- no assays in x")
  
  delim <- if(outname=="") "" else "_"
  
  if(!dir.exists(outdir))
    if(create_outdir){
      message("Creating outdir")
      dir.create(outdir, recursive=TRUE)
    } else stop("Outdir does not exist and create_outdir=FALSE")
  
  if(nrow(colData(x))>0){
    data.table::fwrite(x=data.frame(colData(x)),
                       file=paste0(outdir, "/", outname, delim, "coldata.txt.gz"),
                       compress="gzip", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t",
                       verbose=FALSE)
  }
  
  if(nrow(rowData(x))>0){
    data.table::fwrite(x=data.frame(rowData(x)),
                       file=paste0(outdir, "/", outname, delim, "rowdata.txt.gz"),
                       compress="gzip", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t",
                       verbose=FALSE)
  }
  
  invisible(sapply(asynames, function(y){
    
    asy <- data.frame(assay(x, y))
    data.table::fwrite(x=asy,
                       file=paste0(outdir, "/", outname, delim, y, ".txt.gz"),
                       compress="gzip", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t",
                       verbose=FALSE)
    NULL
    
  }))
  
}

#' A function to run gprofiler2
#' 
#' @param input vector with genes to test, either Ensembl ID or symbol
#' @param background optional background genes to test against
#' @param species mmusculus or hsapiens
#' @param sources source databases, see gprofiler manual
#' @param FDR.threshold a cutoff to filter for significance
#' @param a multiple testing correction method as in \\?gost
#' @param ... further arguments passed to gost function
#'
run_gost <- function (input, background=NULL, species="mmusculus", 
                      sources=c("KEGG", "REAC"), FDR.threshold=0.05, correction_method="fdr", ...) {
  
  suppressMessages({
    require(gprofiler2)
    require(dplyr)
  })
  
  # set archive version to ensure reproducibility
  gprofiler2::set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e101_eg48_p14")
  gg <- gost(query=input, custom_bg=background, 
             organism=species, exclude_iea=TRUE, evcodes=TRUE, 
             user_threshold=FDR.threshold, sources=sources, correction_method=correction_method, ...)$result
  if (is.null(gg)) 
    return(NULL)
  gg <- data.frame(term =gg$term_name, pvalue=gg$p_value, 
                   isize=gg$intersection_size, tsize=gg$term_size, source=gg$source, 
                   genes=gg$intersection)
  gg$genes <- unlist(lapply(gg$genes, 
                            function(x) paste(sort(strsplit(x, split=",")[[1]]), collapse=",")))
  
  # remove the nonsensish root terms from the results
  gg <- gg %>% arrange(pvalue) %>% filter(!grepl("KEGG|REACTOME", term))
  return(gg)
}

#--------------------------------------
# Create standard folders if an "outdir"
# was set
#--------------------------------------

if("outdir" %in% ls()){

  suppressWarnings(dir.create(paste0(outdir, "/plots/"), recursive=TRUE))
  suppressWarnings(dir.create(paste0(outdir, "/lists/"), recursive=TRUE))
  suppressWarnings(dir.create(paste0(outdir, "/rdata/"), recursive=TRUE))  
  suppressWarnings(dir.create(paste0(outdir, "/rcache/"), recursive=TRUE))  
  
  cachedir <- paste0(outdir, "/rcache/")
  
  message(paste("The outdir is:", outdir, sep="\n"))
} else message("[Info] No outdir specified!")

#--------------------------------------
# save.image-like function 
# but with pigz
#--------------------------------------

save.pigz <- function(outname, use_cores=NULL){
  
  os <- Sys.info()[["sysname"]]
  
  if(os!="Windows"){
    
    if(is.null(use_cores)){
      use_cores <- if("mc_workers" %in% ls()) mc_workers else parallel::detectCores()
    }
    
    use_compressor <- if(Sys.which("pigz")!="") "pigz" else "gzip"
  
    con <- pipe(paste0(use_compressor, " -5 -p", use_cores, " > ", outname), "wb")
    save(list=ls(envir=.GlobalEnv), file=con)
    close(con)
  
  } else {
  
    save.image(file=outname)  
    
  }
  
}

# Find closest value of value in vector: https://www.oreilly.com/library/view/the-r-book/9780470510247/ch002-sec020.html
closest <- function(xv,sv) xv[which(abs(xv-sv)==min(abs(xv-sv)))]

#-------------------------------------
# Enrichment test for a set of genes
# relative to a background set for a 
# pathway collection using the hypergeometric
# distrubution
#-------------------------------------

#' @param foreground vector with genes names as foreground
#' @param background vector with background genes, should be all genes that were eligable
#' during the analysis. 
#' @param pathways data.frame with columns 'Term' and 'Genes'. Genes column is supposed to
#' be comma-separated and genes must match the input vectors
#' @param ignore_empty remove pathways that have no matches at all
#' @param min_size only test pathways that have at least that many genes overlapping with the total background set
#' @param all_genes an optional set of gene to consider the universe, e.g. when testing against a subset of a pathway database
#' one would add here all genes in that database -- a vector
enr_test <- function(foreground, background, pathways, ignore_empty=TRUE, min_size=5, all_genes=NULL){
  
  if(sum(c("Term", "Gene") %in% colnames(pathways)) < 2)
    stop("colnames of pathways must be Term and Gene")
  if(length(foreground)>length(background))
    stop("foreground must be equal or smaller than background")
  
  if(sum(foreground %in% background)!=length(foreground))
    stop("foreground must be a subset of background")
  
  # All genes in pathways
  all_genes_pathways <- if(is.null(all_genes)){
    sort(unique(unlist(strsplit(pathways$Genes, ","))))
  } else unique(all_genes)
  
  # all genes in fore- and background 
  all_genes_test <- unique(c(foreground, background))
  
  # the universe
  all_genes_universe <- intersect(all_genes_test, all_genes_pathways)
  
  # Run the test for every term
  phyper.test <- lapply(pathways$Term, function(term_name){
    
    term_gene <- strsplit(pathways[pathways$Term==term_name,]$Gene, ",")[[1]]
    
    # For phyper enrichment test
    # -- x: number of significant genes intersecting with the pathway
    # -- m: number of genes in the pathway and present in our test set
    # -- n: number of all tested genes with annotation anywhere in the pathway collection minus m
    # -- k: significant genes with any annotation in the database  
    .x <- sort(intersect(foreground, term_gene))
    .m <- sort(intersect(background, term_gene))
    .n <- setdiff(all_genes_universe, .m)
    .k <- intersect(foreground, all_genes_universe)
    
    .xl <- length(.x); .ml <- length(.m); .nl <- length(.n); .kl <- length(.k)
    
    pterm <- stats::phyper(.xl, .ml, .nl, .kl, lower.tail=FALSE)
    
    data.frame(term=term_name, pvalue=pterm, isize=.xl, bsize=.ml, tsize=length(term_gene), genes=paste(.x, collapse=","))
    
  }) %>% do.call(rbind, .)
  
  phyper.test <- phyper.test[phyper.test$bsize>min_size,]
  
  phyper.test <- if(ignore_empty) phyper.test[phyper.test$bsize>0,] else phyper.test
  
  # Correct for FDR
  phyper.test$fdr <- stats::p.adjust(phyper.test$pvalue, "BH")
  phyper.test <- phyper.test %>% dplyr::relocate(term, pvalue, fdr) %>% arrange(fdr)
  
  return(phyper.test)
  
}

#========================
# Function to make scatterplots based on sce expression values or colData columns
#========================

#' sce scatter plots
#' 
#' @param sce a SingleCellExperiment
#' @param gene a gene name to extract expr values for, assumed to be "genename" when rownames(sce) are geneid_genename,
#' @param gene_col the color for points with high expression values, low expr is hardcoded as grey
#' @param by use this column from colData for coloring
#' @param text_add logical, whether to add the by label to the plot, e.g. the name of the cluster
#' @param text_col color for text_add
#' @param text_use_label logical, whether to add text as label box rather than just letters (geom_label instead of geom_text)
#' @param legend_title guess what
#' @param embed_onto name of the reducedDim to use as manifold
#' @param embed_dims vector with two integers, the column idx of embed_onto, e.g. c(1,2) when using the first and second UMAP dimensions
#' @param return_data logical, whether to return a list with the plot and the data.frame with the underlying data
#' @param x/y_shift shift labels by these values
#' @param ncol ncols for legend
#' @param nrow nrows for legend
plot_umap <- function(sce, 
                      gene=NULL, 
                      gene_col="darkblue",
                      by=NULL, text_add=TRUE, text_col="black", text_use_label=FALSE,
                      legend_title=NULL,
                      embed_onto="UMAP", embed_dims=c(1,2), # embed_dims are the columns or reducedDims, e.g. UMAP1+2
                      return_data=FALSE, 
                      point_alpha=list.ggplot$umap_alpha,
                      point_size=list.ggplot$pointsize,
                      text_size=list.ggplot$textsize,
                      x_shift=0, y_shift=0, 
                      ncol=NULL, nrow=NULL){
  
  suppressMessages({
    
    require(ggrepel)
    require(SingleCellExperiment)
    require(scater)
    require(tidyverse)
    require(viridis)
    
  })
  
  if(!is(sce, "SingleCellExperiment"))
    stop("sce must be a SingleCellExperiment")
  
  #/ check that only by or gene is used:
  if(!is.null(gene) & !is.null(by))
    stop("Use either by or gene, not both!")
  
  #/ the dim reduction to use as backbone for the plot:
  if(!embed_onto %in% names(reducedDims(sce)))
    stop("embed_onto is not a reducedDim in the sce!")
  
  plotobj <- as.data.frame(reducedDim(sce, embed_onto))
  # otherwise this will be part of the ggplot environment!
  
  if(!all(embed_dims %in% 1:ncol(plotobj)))
    stop("Check embed_dims -- some are out of bounds")
  plotobj <- plotobj[,embed_dims]
  colnames(plotobj) <- paste0("dim", embed_dims)
  
  #/ plot gene expression:
  if(!is.null(gene)){
    
    grepped <- grep(paste0("_", gene, "$"), rownames(sce), value=TRUE)
    if(length(grepped)==0)
      stop("gene not found!")
    
    plotobj$values <- as.numeric(retrieveCellInfo(sce, grepped)$value)
    
    if(is.null(legend_title)) legend_title <- gene
    
    ggobj <- 
      plotobj %>%
      arrange(values) %>%
      ggplot(aes(x=dim1, y=dim2, color=values)) +
      geom_point(alpha=point_alpha, size=point_size) +
      scale_color_gradientn(name=gene, colours=c("grey90", gene_col))
    
  }
  
  #/ plot colData elements:
  if(!is.null(by)){
    
    retrieved_by <- tryCatch(expr=scater::retrieveCellInfo(sce, by),
                             error=function(e) stop("by not found!"))
    
    if(!class(retrieved_by$value) %in% c("integer", "numeric")){
      
      plotobj$group <- factor(retrieved_by$value)  
      
      text_df <- 
        sapply(levels(plotobj$group), function(x){
          data.frame(group=x,
                     x=median(plotobj[plotobj$group==x,"dim1"]) + x_shift,
                     y=median(plotobj[plotobj$group==x,"dim2"]) + y_shift)
        }, simplify=FALSE) 
      
      text_df <- do.call(rbind, text_df)
      
    } else {
      text_add <- FALSE
      plotobj$group <- retrieved_by$value
    }
    
    if(is.null(legend_title)) legend_title <- by
    ggobj <- 
      ggplot(plotobj, aes(x=dim1, y=dim2, color=group)) +
      geom_point(alpha=point_alpha, size=point_size)
    
    if(class(retrieved_by$value) %in% c("numeric", "integer")){
      
      ggobj <- 
        ggobj + viridis::scale_color_viridis(name=legend_title)
      
    } else {
      
      ggobj <- 
        ggobj + 
        scale_color_manual(values=list.ggplot$colorblind_cols) +
        guides(colour=guide_legend(title=legend_title, 
                                   ncol=ncol, nrow=nrow,
                                   override.aes=list(alpha=1, size=text_size)))
      
    }
    
    
    #/ add name of the by into the center of the by points:
    if(text_add) {
      
      if(!text_use_label){
        
        ggobj <- ggobj + geom_text(data=text_df, mapping=aes(x=x, y=y, label=group),
                                   size=text_size, color=text_col)
        
      } else {
        
        ggobj <- ggobj + geom_label_repel(data=text_df, mapping=aes(x=x, y=y, label=group),
                                          size=text_size, color=text_col,
                                          min.segment.length=0, max.overlaps=Inf)
        
      }
    }}
  
  ggobj <- 
    ggobj +
    xlab(paste0(embed_onto, 1)) +
    ylab(paste0(embed_onto, 2))
  
  if(class(plotobj$values) != "numeric"){
    ggobj <- ggobj + theme(legend.position="top", legend.justification="center", legend.direction="horizontal")
  } 
  
  #/ either return plot+data or plot only:
  rm(sce)
  if(return_data){
    return(list(plot=ggobj, data=plotobj))
  } else {
    rm(plotobj)
    return(ggobj)
  }
  
}

#========================
# Normalize by deconvolution
#========================

#/ Wrapper for scran normalization by deconvolution:
norm_scran <- function(sce, seed=1, graph.fun=igraph::cluster_louvain){
  
  suppressMessages({
    require(BiocParallel)
    require(scran)
    require(scuttle)
  })
  
  #/ quick clustering:
  set.seed(seed)
  quicky <- scran::quickCluster(x=sce, assay.type="counts", graph.fun=graph.fun)
  
  #/ get sum factors:
  set.seed(seed)
  sce <- scran::computeSumFactors(sce, clusters=quicky)
  
  set.seed(seed)
  sce <- scuttle::logNormCounts(sce, assay.type="counts", size_factors=sce$sizeFactor)
  return(sce)
  
}

#========================
# Convenience functions related to the single-cell workflow:
#========================

#/ graph-based clustering wrapper scran/iGraph
GBclust <- function(x, k=10, type="jaccard", 
                    rseed=1, clusterfun=igraph::cluster_louvain){
  
  suppressMessages({
    require(bluster)
    require(BiocParallel)
    require(igraph)
  })
  
  set.seed(1)
  g <- bluster::makeSNNGraph(x, type=type, k=k, BPPARAM=SerialParam())
  
  set.seed(rseed)
  cl <- clusterfun(g)
  
  return(cl$membership)
  
}

#========================
# UCell to vector
#========================

convert_ucell <- function(ucell){
  
  lapply(1:nrow(ucell), function(x) {
    
    u <- ucell[x,]
    umax <- which(u==max(u))
    if(length(umax)>1) return(NA)
    return(gsub("_UCell", "", names(umax)))
    
  }) -> uc
  
  unlist(uc)
  
}

#------------------------------------
# Customized version of fgsea::plotEnrichment
#------------------------------------

plot_fgsea <- function(pathway, stats, gseaParam=1, ticksSize=0.3,
                       lwd=0.5, color_line="darkorchid", color_segment="black",
                       xlab="ranked genes", ylab="enrichment score", add_cutoff=FALSE,
                       add_caption=FALSE){
  
  suppressMessages({
    require(fgsea)
    require(tidyverse)
  })
  
  l_pathway <- length(pathway)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats=pathway, 
                          returnAllExtremes=TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  diff <- (max(tops) - min(bottoms))/8
  x=y=NULL
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  minpos <- min(stats[stats>0])
  maxneg <- max(stats[stats<0])
  if(length(minpos) > 0) use_cutoff <- minpos else use_cutoff <- maxneg
  
  g <- 
    ggplot(toPlot, aes(x=x, y=y)) + 
    geom_hline(yintercept=0, color="black", lwd=lwd) +
    geom_line(color=color_line) + 
    geom_segment(data=data.frame(x=pathway), 
                 mapping=aes(x=x, y=-diff/2, xend=x, yend=diff/2), 
                 size=ticksSize, color=color_segment) + 
    labs(x=xlab, y=ylab)
  
  if(add_cutoff){
    g <- 
      g +
      geom_segment(data=data.frame(x=as.numeric(which(stats==use_cutoff)[1])), 
                   mapping=aes(x=x, y=-diff*1.1, xend=x, yend=diff*1.1), 
                   size=ticksSize*2)
  }
  g + 
    labs(caption=paste(l_pathway, "genes in pathway --", length(pathway), "found in ranking"))
  if(add_caption) {
    return(g + theme(plot.caption=element_text(hjust=0)))
  } else return(g)
}

# UpSetR::fromList but preserving rownames
fromlist <- function(input){
  
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  rownames(data) <- elements
  return(data)
  
}

# https://stackoverflow.com/a/6463946/5074060
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

# Given a matrix-like object with genes being rows and samples being columns,
# and a group information for the columns, calculate the percentage of columns per group
# that have counts > threshold per gene. Efficiently done via base R.
get_pexpr <- function(data, group, threshold=0, digits=2){
  
  if(ncol(data)!=length(group)) stop("ncol(data) != length(group)")
  if(!is.numeric(threshold) | threshold < 0) stop("threshold must be numeric and > 0")
  
  datar <- (data>threshold) * 1
  a <- base::rowsum(x=t(datar), group=group)
  b <- as.numeric(table(group)[rownames(a)])
  f <- base::round(100*t(apply(a, 2, function(x) x/b)), digits=digits)
  f
  
}

#--------------------------------------
# ggplot defaults
#--------------------------------------

list.ggplot <- list()

# basic palette of 11 colors -- 10 are actual colors plus grey and black
list.ggplot$colorblind_cols <-  c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                                  "#0072B2", "#D55E00", "#CC79A7", "#661100", 
                                  "#332288", "gray40", "black")
# make brighter and darker
list.ggplot$colorblind_cols <- 
  c(list.ggplot$colorblind_cols,
    colorspace::lighten(list.ggplot$colorblind_cols, .5),
    colorspace::lighten(list.ggplot$colorblind_cols, .9))

# some global options for ggplot, to achieve some standard appearances:
list.ggplot$theme <- ggplot2::theme_bw                  # default theme
list.ggplot$themesize        <- 12                      # base_size option
list.ggplot$pointsize        <- 0.5                     # point size for UMAPs
list.ggplot$pointsize_larger <- 3                       # for dotplots with few points
list.ggplot$textsize         <- 3                       #
list.ggplot$ablinesize       <- .5                      #
list.ggplot$linesize         <- 1                       #
list.ggplot$umap_alpha       <- 0.6                     # for UMAPs: alpha of points
list.ggplot$legendsize       <- 3

# set the global theme:
theme_set(list.ggplot$theme(base_size=list.ggplot$themesize))
theme_update(axis.text=element_text(size=list.ggplot$themesize),
             axis.title=element_text(size=list.ggplot$themesize),
             legend.text=element_text(size=list.ggplot$themesize),
             legend.title=element_text(size=list.ggplot$themesize))

base::options(ggplot2.discrete.colour=list.ggplot$colorblind_cols, 
              ggplot2.discrete.fill=list.ggplot$colorblind_cols)


#!/usr/bin/env Rscript

# Author: Sean Maden
#
# Preprocess available sample metadata
#

require(data.table); require(rjson)
library(recountmethylationManuscriptSupplement)

library(recountmethylationManuscriptSupplement)
library(HDF5Array)
library(recountmethylation)
library(wateRmelon) # for age predictions
library(minfi) # for cell type and sex predictions

library(recountmethylationManuscriptSupplement)
library(recountmethylation)
library(HDF5Array)
library(minfi)
library(ewastools)


#--------------
# dnam metadata
#--------------
# model-based predictoins

#' Make model-based metadata predictions 
#'
#' Make model-based metadata predictions using DNAm assays. Predictions for 
#' age, sex, and blood cell fractions are produced.
#' 
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param rgset.fname Name of the compilation file containing red signals 
#' extracted from IDATs (character).
#' @param grset.fname Name of the compilation file containing green signals 
#' extracted from IDATs (character).
#' @param mdmod.fname Name of the table of model-based predictions produced 
#' ("mdmod_dnam-predictions").
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @param md.dname Name of directory, in files.dname, containing the instance
#' metadata files ("metadata).
#' @param nsamp.block Number of samples per data block processed (integer, 50).
#' @param rgset.path Path 
#' @param compilations.dname Name of directory, in files.dname, containing the
#' compilation files with red and green signals for prediction calculations 
#' ("compilations")
#' @return NULL, produces table of model-based metadata predictions.
#' @seealso get_qcmetrics(); md_preprocess(); md_postprocess()
#' @export
md_predictions <- function(ts, rgset.fname, grset.fname,
                           mdmod.fname = "mdmod_dnam-predictions", 
                           files.dname = "recount-methylation-files",
                           md.dname = "metadata", nsamp.block = 50,
                           compilations.dname = "compilations"){
  md.dpath <- file.path(files.dname, md.dname)
  rgset.path <- file.path(files.dname, compilations.dname, rgset.fname)
  grset.path <- file.path(files.dname, compilations.dname, grset.fname)
  if(!file.exists(rgset.path)){stop("Couldn't locate red signals compilation ",
                                    "file at:\n",rgset.path)}
  if(!file.exists(grset.path)){stop("Couldn't locate red signals compilation ",
                                    "file at:\n",grset.path)}
  message("Loading h5se datasets...")
  rgset <- loadHDF5SummarizedExperiment(rgset.path)
  grset <- loadHDF5SummarizedExperiment(grset.path)
  mdmod <- matrix(nrow = 0, ncol = 8)
  message("Getting model predictions for blocks of ",nsamp.block," samples...")
  if(!ncol(rgset) == ncol(grset)){
    message("Warning, ncol not identical for rgset and grset. Choosing the ",
            "lowest column count for sample indexing...")
    col.tot <- ifelse(ncol(rgset) > ncol(grset), ncol(grset), ncol(rgset))
  } else{ncol.tot <- ncol(rgset)};blocks <- getblocks(ncol.tot, nsamp.block)
  for(b in blocks){
    rgf <- rgset[, unlist(b)]
    rgf.matrix <- RGChannelSet(Green = as.matrix(getGreen(rgf)),
                               Red = as.matrix(getRed(rgf)),
                               annotation = annotation(rgf))
    celltypepred <- estimateCellCounts(rgf.matrix)
    msf <- mapToGenome(preprocessRaw(rgf));sexpred <- getSex(msf)
    grf <- grset[,colnames(rgf)];predage <- agep(getBeta(grf))
    mdf <- cbind(predage, cbind(sexpred[,3], celltypepred))
    mdmod <- rbind(mdmod, mdf);message("Completed index ", max(unlist(b)))}
  message("Formatting mdmod...")
  mdmod <- as.data.frame(mdmod, stringsAsFactors = FALSE)
  colnames(mdmod) <- c("predage", "predsex", 
                       paste0("predcell.", colnames(mdmod)[3:8]))
  mdmod$gsm <- gsub("\\..*", "", rownames(mdmod))
  mdmod <- mdmod[!duplicated(mdmod$gsm),]
  mdmod.fpath <- file.path(md.dpath, paste0(mdmod.fname,"_", ts, ".rda"))
  message("Saving mdmod to ", mdmod.fpath, "...")
  save(mdmod, file = paste0(table.fn, ".rda")); return(NULL)
}

# get qc metrics from red/grn signals

#' Get quality metrics from DNAm assays
#'
#' Calculates quality metrics from DNAm assays contained in the compilation 
#' files. These include signals for 17 BeadArray controls, 2 controls for the 
#' methylated and unmethylated signals, and predicted genotypes by sample
#' using methods in the ewastools package.
#' 
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param mdqc.fname Name of the quality metrics table output ("mdqc")
#' @param athresh Similarity threshold (percent similarity) for the ewastools 
#' predicted genotype (decimal, 0.1).
#' @param nsamp.block Samples per data block processed (integer, 50).
#' @param md.dname Name of directory, in files.dname, containing the instance
#' metadata files ("metadata).
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @return NULL, produces a table of quality metrics.
#' @seealso md_predictions(); md_preprocess(); md_postprocess()
#' @export
get_qcmetrics <- function(ts, mdqc.fname = "mdqc", athresh = 0.1,
                          nsamp.block = 50, md.dname = "metadata",
                          files.dname = "recount-methylation-files"){
  mdqc.fpath <- file.path(md.dpath, paste0(mdqc.fname, "_", ts, ".rda"))
  message("Working on file: ", mdqc.fpath, "...")
  md.dpath <- file.path(files.dname, md.dname)
  blocks <- getblocks(slength = ncol(rgset), bsize = nsamp.block)
  ms <- matrix(nrow = 0, ncol = 20)
  cdf <- as.data.frame(getProbeInfo(rgset, type = "Control"))
  message("Getting quality metrics from blocks of ",nsamp.block," samples...")
  for(bi in seq(length(blocks))){
    b <- blocks[[bi]];rgf <- rgset[, b]; 
    colnames(rgf) <- gsub("\\..*", "", colnames(rgf))
    redsignal <- getRed(rgf); greensignal <- getGreen(rgf)
    basignals <- bactrl(rs = t(redsignal), gs = t(greensignal), cdf = cdf)
    mset <- preprocessRaw(rgf);ms <- getMeth(mset); us <- getUnmeth(mset);
    meth.l2med <- apply(ms, 2, function(x){log2(median(x))})
    unmeth.l2med <- apply(us, 2, function(x){log2(median(x))})
    mi <- cbind(basignals, data.frame(meth.l2med = meth.l2med, 
                                      unmeth.l2med = unmeth.l2med,
                                      stringsAsFactors = FALSE))
    ms <- rbind(ms, mi); message("Finished block ", bi)}
  message("Detecting replicates from genotypes using ",
          "similarity threshold ",athresh,"...")
  snp1.info <- getProbeInfo(rgset, type = "SnpI")
  snp2.info <- getProbeInfo(rgset, type = "SnpII")
  snp.addr <- c(snp1.info$AddressA, snp1.info$AddressB, snp2.info$AddressA)
  beta.snp <- getSnpBeta(rgset)
  colnames(beta.snp) <- gsub("\\..*", "", colnames(beta.snp))
  ama = matrix(nrow = 0, ncol = 3) # main identity matrix
  colnames(ama) <- c("gsm", "gseid", "cgsnp.gsm.geno")
  for(id in unique(md$gseid)){
    gsmv <- md[md$gseid == id,]$gsm
    sbi <- beta.snp[, colnames(beta.snp) %in% gsmv, drop = FALSE]
    sbi <- as.matrix(sbi); class(sbi) <- "numeric"
    sgenoi <- call_genotypes(sbi, learn = FALSE)
    sagree <- check_snp_agreement(sgenoi, donor_ids = colnames(sbi), 
                                  sample_ids = colnames(sbi))
    for(sai in sagree){
      am <- data.frame(gsm = unique(c(sai$sample1, sai$sample2)), 
                       stringsAsFactors = FALSE)
      am$gseid <- id; am$gsm.geno <- "NA"
      for(gsm in am$gsm){
        sai.cond <- (sai$sample1==gsm | sai$sample2==gsm) & 
          sai$agreement > athresh
        sai.gsm <- sai[sai.cond,]
        sai.id <- unique(c(sai.gsm$sample1, sai.gsm$sample2))
        am[am$gsm == gsm,]$gsm.geno <- paste(sai.id, collapse = ";")
      };ama <- rbind(ama, am)
    };message("finished study ", id)
  }
  ama$num.shared <- unlist(lapply(ama$gsm.geno, function(x){
    length(unique(unlist(strsplit(x, ";"))))}))
  message("Appending replicates info...")
  mdqc <- ms; ama <- ama[ama$gsm %in% mdqc$gsm,]
  ama <- ama[order(match(ama$gsm, mdqc$gsm)),]
  if(identical(ama$gsm, mdqc$gsm)){
    mdqc$cgsnp.gsm.geno <- ama$gsm.geno;mdqc$cgsnp.nshared <- ama$num.shared
  } else{stop("Couldn't match GSM IDs for mdqc and genotypes.")}
  save(mdqc, file = mdqc.fpath);return(NULL)
}

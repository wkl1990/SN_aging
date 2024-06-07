library(tidyverse)
library(rlang)
library(rlang)
library(data.table)
library(reticulate)
library(furrr)
plan(multicore, workers = 5)
options(future.globals.maxSize= 15e9)
library(future.apply)

library(gamlss.dist)
library(gamlss)
library(fitdistrplus)

#rdir <- file.path(here::here(), "package/R")
#import::from(.from = "cembav2env.R", .directory = rdir,
#  cembav2env, Sa2PeakCalling)

usePython <- "/bin/python"
reticulate::use_python(usePython)
pd <- reticulate::import("pandas", convert = FALSE)

# * meta
#atacMeta <- readRDS(cembav2env$sa2metaFile)
projdir <- here::here()
atacmetafnm <- file.path(projdir, "rds/RNA/cluster/integration/meta/SN.allen.integration.meta.final.barcode.csv")
atacMeta <- read.csv(atacmetafnm, row.names=1)
atacMeta$barcode2 <- paste(atacMeta$sampleID, atacMeta$barcode, sep=".")
rownames(atacMeta) <- atacMeta$barcode2

workdir <- file.path(projdir, "peak_calling")
default_cutoff <- 0.001
pep <- 0.001
rnd_upbound <- 0.1
qval <- 0.01
outdir <- file.path(workdir, "after_integra/scfilter/")

fitbgdir <- file.path(workdir, "after_integra/scfilter/",
  "fitfrac_bg")
fitbgsuffix <- ".fitPeakModel.para.csv"

loadfitbg <- function(l4) {
  data.table::fread(
    file = file.path(fitbgdir, paste0(l4, fitbgsuffix)),
    sep = ",",
    header = TRUE,
    data.table = FALSE)
}


# * read union peak
peakfrac_union <- pd$read_pickle(
  file.path(workdir, "after_integra/scfilter/", "peakfrac_union.pkl"))
peakfrac <- reticulate::py_to_r(x = peakfrac_union)

# * check all the models are converged
# all of them converged
subclasses <- rownames(peakfrac)
#subclass_names <- gsub(" ", "_", subclasses)
checkConverge <- future_map_lgl(subclasses,
  \(subclass) {
    m <- loadfitbg(gsub(" ", "_", subclass))
    if(m$converged > 0) {
      TRUE
    } else {
      FALSE
    }
  }
)

# * fit model per pL4
pvalues <- future_map(.x = subclasses, .f = \(subclass) {
  m <- loadfitbg(gsub(" ", "_", subclass))
  p <- 1 - gamlss.dist::pBEZI(
    unlist(peakfrac[subclass, ]), mu = m$mu, sigma = m$sigma,
    nu = m$nu, lower.tail = TRUE, log.p = FALSE)
  return(p)
})
names(pvalues) <- subclasses

# check pvalues
## (subclass <- sample(subclasses, 1))
## p1 <- pvalues[[subclass]]
## m <- loadfitbg(subclass)
## t1 <- 1-gamlss.dist::pBEZI(
##   unlist(peakfrac[subclass, ]), mu = m$mu, sigma = m$sigma,
##   nu = m$nu, lower.tail = TRUE, log.p = FALSE)
## sum(p1 == t1) == length(p1)

## save pvalues as matrix
pvaluesMat <- do.call(rbind, pvalues)
rownames(pvaluesMat) <- subclasses
colnames(pvaluesMat) <- colnames(peakfrac)
# > 6G
saveRDS(pvaluesMat, file = file.path(outdir, "subclass.unionpeak.pvalue.mat.rds"))

# * qvalue and filter
fullpeaksetFile = file.path(projdir,
      "peak_calling/after_integra/macs2/parsePeak/",
      "SN_integra.macs2.filteredNfixed.union.peakSet")
readUnionPeakSetFile <- function(sep1 = ":", sep2 = "-") {
      r <- data.table::fread(file = fullpeaksetFile,
        sep = "\t", header = TRUE, data.table = FALSE)
      rownames(r) <- with(r, paste(seqnames,
        paste(start, end, sep = sep2), sep = sep1))
      # make sure bed files are sorted
      r <- r[order(r$seqnames, as.numeric(r$start)), ]
      return(r)
    }

unionpeaks <- readUnionPeakSetFile()
spm <- 5
filterUnionPeaks <- rownames(unionpeaks)[unionpeaks$spm >= spm]
filterPvalMat <- pvaluesMat[, filterUnionPeaks]
#identical(filterPvalMat,pvaluesMat) #TRUE

qvaluesUnion <- vapply(subclasses, FUN = \(subclass) {
  message(subclass)
  p.adjust(filterPvalMat[subclass, ], method = "BH")
}, FUN.VALUE = filterPvalMat[1, ])

qvaluesUnion <- t(qvaluesUnion)

min_qval <- furrr::future_map_dbl(.x = seq_len(ncol(qvaluesUnion)),
  .f = \(i) {
    min(qvaluesUnion[, i])
  }
)
names(min_qval) <- colnames(qvaluesUnion)

scfilteredPeaks <- names(min_qval)[min_qval <= qval]

saveRDS(scfilteredPeaks, file=file.path(outdir, "scfilteredPeaks.rds"))
## not all the peaks pass the sc-filtering

# check all union peak
qvaluesAll <- future.apply::future_vapply(
  subclasses, FUN = \(subclass) {
    message(subclass)
    p.adjust(pvaluesMat[subclass, ], method = "BH")
  }, FUN.VALUE = pvaluesMat[1, ])
colnames(qvaluesAll) <- subclasses

## all the peaks pass the sc-filtering
min_qval_all <- furrr::future_map_dbl(.x = seq(nrow(qvaluesAll)),
  .f = \(i) {min(qvaluesAll[i, ])})
sum(min_qval_all <= qval)

# * test
## i <- 1
## pL4 <- rownames(peakfrac)[i]
## m <- loadfitbg(pL4)
## p <- 1 - gamlss.dist::pBEZI(
##   unlist(peakfrac[i, ]), mu = m$mu, sigma = m$sigma,
##   nu = m$nu, lower.tail = TRUE, log.p = FALSE)
## q <- p.adjust(p, method = "BH")


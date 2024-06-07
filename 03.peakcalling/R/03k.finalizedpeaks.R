library(tidyverse)
library(rlang)
Sys.setenv("_R_USE_PIPEBIND_" = TRUE)
library(future)
library(data.table)

projdir <- here::here()
#rdir <- file.path(projdir, "package/R")
#import::from(.from = "cembav2env.R", .directory = rdir,
#  cembav2env, Sa2PeakCalling)

workdir <- file.path(projdir, "peak_calling")
peaksumdir <- file.path(workdir, "after_integra/final")

# * prepare finalized peaks (SPM >= 5)
#unionpeaks <- Sa2PeakCalling$readUnionPeakSetFile()
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

#finalpeaks <- subset(unionpeaks, spm >= 5)
scfilteredPeaks <- readRDS(file=file.path(workdir, "after_integra/scfilter/", "scfilteredPeaks.rds"))
finalpeaks <- unionpeaks[which(rownames(unionpeaks) %in% scfilteredPeaks), ]

fpknms <- rownames(finalpeaks)
fpkbed <- data.frame(
  chrom = finalpeaks$seqnames,
  start = finalpeaks$start,
  end = finalpeaks$end,
  name = fpknms
)
rownames(fpkbed) <- fpkbed$name

write.table(fpkbed, file = file.path(peaksumdir,
  "SN_integra.final.peak.srt.bed"), quote = FALSE,
  sep = "\t", col.names = FALSE, row.names = FALSE)

# * load finalized peak
#fpkbed <- Sa2PeakCalling$loadFinalPeakBed()
finalpeakBedFile = file.path(peaksumdir,
      "SN_integra.final.peak.srt.bed")
readBed4File = function(fnm,
                           sep1 = ":", sep2 = "-", header = FALSE,
                           sep = "\t") {
      beds <- data.table::fread(fnm, header = header, sep = "\t",
        data.table = FALSE) |>
        setNames(object = _, nm = c("chr", "start", "end", "name"))
      rownames(beds) <- beds$name
      return(beds)
    }

fpkbed <- readBed4File(finalpeakBedFile)

fpknms <- fpkbed$name

# * preprare finalized peaks for each subclass (binary pmat)
#subclassunionpeaks <- Sa2PeakCalling$loadpL4UnionPeakBeds()
readBedFile <- function(fnm,
                           sep1 = ":", sep2 = "-", header = FALSE,
                           sep = "\t") {
      beds <- data.table::fread(fnm, header = header, sep = "\t",
        data.table = FALSE) |>
        setNames(object = _, nm = c("chr", "start", "end"))
      rownames(beds) <- with(beds,
        paste(chr, paste(start, end, sep = sep2), sep = sep1))
      return(beds)
    }
subclassUnionPeakBedDir = file.path(here::here(), "peak_calling", "after_integra/macs2/celltype/")
loadpL4UnionPeakBeds <- function() {
      files <- list.files(path = subclassUnionPeakBedDir,
        pattern = "SN_integra.macs2.celltype.*.bed")
      cls <- basename(files) |>
        gsub(".bed", "", x = _) |>
        gsub("SN_integra.macs2.celltype.", "", x = _)
      bedList <- lapply(cls, \(x) {
        fnm <- file.path(subclassUnionPeakBedDir,
          paste0("SN_integra.macs2.celltype.", x, ".bed"))
        readBedFile(fnm)
      })
      names(bedList) <- cls
      return(bedList)
    }
subclassunionpeaks <- loadpL4UnionPeakBeds()

groupSubclassnms <- names(subclassunionpeaks)
cleanpSubclassnms <- gsub("^.*\\.", "", groupSubclassnms)
subclassfinalpeaks <- lapply(seq_along(subclassunionpeaks), \(i) {
  subclass <- names(subclassunionpeaks[i])
  rawp <- rownames(subclassunionpeaks[[i]])
  r <- fpknms[fpknms %in% rawp]
  message(length(r), " final peaks in ", subclass)
  return(r)
})
names(subclassfinalpeaks) <- groupSubclassnms

# check we have all the peaks
# all_peaks <- unique(unlist(subclassfinalpeaks))

outsubclassfinalpeak_dir <- file.path(peaksumdir, "subclass.final.peak.srt")
if (!dir.exists(outsubclassfinalpeak_dir)) {
  dir.create(outsubclassfinalpeak_dir)
}

subclassfp_df <- lapply(seq_along(subclassfinalpeaks), \(i) {
  rawnm <- names(subclassfinalpeaks)[i]
  nm <- cleanpSubclassnms[i]
  fnm <- file.path(outsubclassfinalpeak_dir, paste0(nm, ".bed"))
  message(nm, " from ", rawnm, " to outfile: ", fnm)
  rdf <- fpkbed[subclassfinalpeaks[[i]], ]
  write.table(rdf, file = fnm, quote = FALSE, sep = "\t",
    col.names = FALSE, row.names = FALSE)
  return(rdf)
})

# add subclass binary pmat
# copy from sa2.Fig2.R in figures
## construct subclass binary pmat.
atacmetafnm <- file.path(projdir, "rds/RNA/cluster/integration/meta/SN.allen.integration.meta.final.barcode.csv")
atacMeta <- read.csv(atacmetafnm, row.names=1)
atacMeta$barcode2 <- paste(atacMeta$sampleID, atacMeta$barcode, sep=".")
rownames(atacMeta) <- atacMeta$barcode2

subclasses <- unique(atacMeta$subclass_label_id)
subclasspeakdir <- file.path(projdir,
  "peak_calling/after_integra/final/subclass.final.peak.srt")
subclassPeaks <- lapply(subclasses, \(subclass) {
  fnm <- file.path(subclasspeakdir, paste0(gsub(" ", "_", subclass), ".bed"))
  readBed4File(fnm)
})
names(subclassPeaks) <- subclasses

peakBed <- fpkbed

bpmat <- vapply(subclasses, \(x) {
  r <- rep(0, nrow(peakBed))
  r[match(subclassPeaks[[x]]$name , peakBed$name)] <- 1
  return(r)
}, rep(0, nrow(peakBed))) |> t()
rownames(bpmat) <- subclasses
colnames(bpmat) <- peakBed$name
saveRDS(bpmat, file.path(peaksumdir, "SN_integra.binary.pmat.subclass.rds"))



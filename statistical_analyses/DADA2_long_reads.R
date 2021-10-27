if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

### set environment
library(dada2);packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")

path <- "permeability/PACB/DADA2" 
setwd(path)
path.out <- "Figures/"
path.rds <- "RDS/"
fns <- list.files(path, pattern="fq", full.names=TRUE)
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())
sample.names <- sapply(strsplit(basename(fns), ".fq"), `[`, 1)

### remove primers and orient reads

nops <- file.path(path, "noprimers", basename(fns))
for(i in seq_along(fns)) {
  fn <- fns[[i]]; nop <- nops[[i]]
  dada2:::removePrimers(fn, nop, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE, verbose=TRUE)
}
# generate histogram
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
postscript("hist.eps",hor=T)
hist(lens, 100)
dev.off()

### quality filtering using length and EE=6
filtsLen <- file.path(path, "noprimers", "filteredLen", basename(fns))
filtsMinQ <- file.path(path, "noprimers", "filteredMinQ", basename(fns))
trackLen <- filterAndTrim(nops, filtsLen, minQ=3, minLen=1430, maxLen=1600, maxN=0, rm.phix=FALSE)
trackMinQ <- filterAndTrim(nops, filtsMinQ, minQ=3, maxN=0, rm.phix=FALSE)
filts <- file.path(path, "noprimers", "filtered", basename(fns))
filts_track <- filterAndTrim(nops, filts, minQ=3, minLen=1430, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)

write.csv(file="filter_stat.csv", filts_track)
write.csv(file="filter_statEE2.csv", trackEE2)
write.csv(file="filter_statMinQ.csv", trackMinQ)

### run DADA2

# 1. dereplicate
drp <- derepFastq(filts, verbose=TRUE)
err <- learnErrors(drp, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
postscript("err.eps",hor=T)
plotErrors(err)
dev.off()
saveRDS(err, file.path(path.rds, "err.rds"))
err=readRDS("RDS/err.rds")

# 2. denoise
dd <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE) 
saveRDS(dd, file.path(path.rds, "dd.rds"))
dd=readRDS("RDS/dd.rds")
st <- makeSequenceTable(dd); dim(st)
saveRDS(st, file.path(path.rds, "st.rds"))
st=readRDS("RDS/st.rds")

# 3. remove chimera
# Higher MFPOA to avoid flagging intra-genomic variants
bim <- isBimeraDenovo(st, minFoldParentOverAbundance=4.5, multithread=TRUE)
table(bim) 
rbi=removeBimeraDenovo(st, minFoldParentOverAbundance=4.5, multithread=TRUE)
saveRDS(rbi, file.path(path.rds, "rbi.rds"))
dim(rbi) #[1]   187 1376
rbi=readRDS("RDS/rbi.rds")

#4. get stats
getN <- function(x) sum(getUniques(x))
track <- cbind(filts_track, sapply(dd, getN), rowSums(rbi))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(file="total_stat.csv", track)


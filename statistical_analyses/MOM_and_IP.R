##
## IP vs total valume of mother's milk
##

library(rstan)

setwd("leaky_gut/data")

st <- read.csv("subject.csv", row.names=1)

save(st, file="subject_table.rda")
load("subject_table.rda")

mm.dir <- "../pics/mm_dir"
if (!file.exists(mm.dir))
    dir.create(mm.dir)

## IP vs mother's milk

ip <- st$IP
mm <- st$doseMBM
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]


file <- paste0(mm.dir,"/ip_vs_mothers_milk_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
par(op)
dev.off()
system(paste0("open ",file))


## totalBM - total mom's breast milk + donor

tBM <- st$totalBM
ip <- st$IP
idx <- is.finite(ip) & is.finite(tBM)
ip <- ip[idx]
tBM <- tBM[idx]
idx <- tBM<600
ip <- ip[idx]
tBM <- tBM[idx]


file <- paste0(mm.dir,"/ip_vs_mother_plus_donor_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(tBM,ip, xlab="Cumulative mother's milk and donor's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
par(op)
dev.off()
system(paste0("open ",file))



## IP vs mother's milk with percentage of formula

ip <- st$IP
mm <- st$doseMBM
fp <- as.numeric(gsub("\\%","",st$FormulaPerc))
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
fp <- fp[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
fp <- fp[idx]

myHist(fp)

fp20 <- ifelse(fp>20,1,0)+1
fp40 <- ifelse(fp>40,1,0)+1
fp60 <- ifelse(fp>60,1,0)+1
fp80 <- ifelse(fp>80,1,0)+1

fp.col.tbl <- c('black', 'red')
fp.pch.tbl <- c(1, 1.5)

file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_pFormula20_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- fp20==2
points(mm[idx],ip[idx], pch=19, col='red', cex=1.25)
par(op)
dev.off()
system(paste0("open ",file))


file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_pFormula80_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- fp80==2
points(mm[idx],ip[idx], pch=19, col='red', cex=1.25)
par(op)
dev.off()
system(paste0("open ",file))


## dose of formula

myHist(log10(st$doseFormula+1))

log10.f <- log10(st$doseFormula+1)
dbm.bin <- ifelse(st$doseDBM>0,1,0)
ip <- st$IP
mm <- st$doseMBM
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
log10.f <- log10.f[idx]
dbm.bin <- dbm.bin[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
log10.f <- log10.f[idx]
dbm.bin <- dbm.bin[idx]



file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_doseFormula_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- log10.f>0
points(mm[idx],ip[idx], pch=19, col='red', cex=log10.f[idx])
par(op)
dev.off()
system(paste0("open ",file))


##  with percentage of formula and percentage of donor's milk

ip <- st$IP
mm <- st$doseMBM
fp <- as.numeric(gsub("\\%","",st$FormulaPerc))
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
fp <- fp[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
fp <- fp[idx]

myHist(fp)

myHist(log10(fp+1))


ip <- st$IP
mm <- st$doseMBM

describe(st$feedingGroup)

log10.fp <- log10(as.numeric(gsub("\\%","",st$FormulaPerc))+1)
log10.dbmp <- log10(as.numeric(gsub("\\%","",st$DBMperc))+1)
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
log10.fp <- log10.fp[idx]
log10.dbmp <- log10.dbmp[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
log10.fp <- log10.fp[idx]
log10.dbmp <- log10.dbmp[idx]

file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_pFormula_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- log10.fp>0
points(mm[idx],ip[idx], pch=19, col='red', cex=log10.fp[idx]/2)
points(mm[!idx],ip[!idx], pch=19, col='black',cex=0.5)
par(op)
dev.off()
system(paste0("open ",file))


file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_pFormula_marking_mothers_with_donor_milk_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- log10.fp>0
points(mm[idx],ip[idx], pch=19, col='red', cex=log10.fp[idx]/2)
points(mm[!idx],ip[!idx], pch=19, col='black',cex=0.5)
idx <- log10.dbmp>0
points(mm[idx],ip[idx], pch=19, col='green', cex=log10.dbmp[idx]/2)
legend('topright',legend=c("only mother's milk",'donor','formula'), pch=19, col=c('black','green','red'), inset=0.05)
par(op)
dev.off()
system(paste0("open ",file))

ip <- st$IP
mm <- st$doseMBM
fp <- as.numeric(gsub("\\%","",st$FormulaPerc))
fGroup <- st$feedingGroup
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
fp <- fp[idx]
fGroup <- fGroup[idx]
## idx <- mm<600
## ip <- ip[idx]
## mm <- mm[idx]
## fp <- fp[idx]

describe(st$feedingGroup)
## Value          D     F   MBM   MDM   MFM  None
## Frequency      2     7    63    12    31     1
## Proportion 0.017 0.060 0.543 0.103 0.267 0.009


ip <- st$IP
mm <- st$doseMBM
fGroup <- st$feedingGroup
log10.fp <- log10(as.numeric(gsub("\\%","",st$FormulaPerc))+1)
log10.dbmp <- log10(as.numeric(gsub("\\%","",st$DBMperc))+1)
idx <- is.finite(ip) & is.finite(mm)
fGroup <- fGroup[idx]
ip <- ip[idx]
mm <- mm[idx]
log10.fp <- log10.fp[idx]
log10.dbmp <- log10.dbmp[idx]
## idx <- mm<600
## ip <- ip[idx]
## mm <- mm[idx]
## log10.fp <- log10.fp[idx]
## log10.dbmp <- log10.dbmp[idx]




file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_pFormula_marking_mothers_with_donor_milk_loess_plot_v2.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP", pt.col='white')
abline(v=140,col='gray')
abline(h=0.05,col='gray')
idx <- log10.fp>0
points(mm[idx],ip[idx], pch=20, col='red', cex=0.3)
idx <- log10.fp>0 & mm==0
points(mm[idx],ip[idx], pch=19, col='red', cex=1)
idx <- log10.fp>0 & mm>0
points(mm[idx],ip[idx], pch=1, col='red', cex=1.15*log10.fp[idx])
## idx <- !is.na(fGroup) & fGroup!='None' & !(log10.fp>0)
idx <- fGroup=="MBM" & mm>0
points(mm[idx],ip[idx], pch=19, col='black',cex=1)
idx <- log10.dbmp>0
points(mm[idx],ip[idx], pch=20, col='green', cex=0.3)
idx <- log10.dbmp>0 & mm==0
points(mm[idx],ip[idx], pch=19, col='green', cex=1)
idx <- log10.dbmp>0 & mm>0
points(mm[idx],ip[idx], pch=1, col='green', cex=1.15*log10.dbmp[idx])
idx <- !is.na(fGroup) & fGroup=='None'
##points(mm[idx],ip[idx], pch=1, col='white',cex=2)
points(mm[idx],ip[idx], pch=1, col='black',cex=1)
legend('topright',legend=c("only mother's milk",'donor','formula','combined donor', 'combined formula'),
       pch=c(19,19,19,1,1), pt.cex=c(1,1,1,2,2),
       col=c('black','green','red', 'green', 'red'), inset=0.05)
par(op)
dev.off()
system(paste0("open ",file))


table(log10.fp>0, log10.dbmp)


##  with doseage of formula and dosage of donor's milk

ip <- st$IP
mm <- st$doseMBM
log10.fd <- log10(st$doseFormula+1)
log10.dbmd <- log10(st$doseDBM+1)
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
log10.fd <- log10.fd[idx]
log10.dbmd <- log10.dbmd[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
log10.fd <- log10.fd[idx]
log10.dbmd <- log10.dbmd[idx]

file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_dFormula_dDonor_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- log10.fd>0
points(mm[idx],ip[idx], pch=19, col='red', cex=log10.fd[idx]/2)
points(mm[!idx],ip[!idx], pch=19, col='black',cex=0.5)
idx <- log10.dbmd>0
points(mm[idx],ip[idx], pch=19, col='green', cex=log10.dbmd[idx]/2)
legend('topright',legend=c("only mother's milk",'donor','formula'), pch=19, col=c('black','green','red'), inset=0.05)
par(op)
dev.off()
system(paste0("open ",file))



##
##   with percentage of formula and percentage of donor's milk
##


o <- order(x)
x <- x[o]
y <- y[o]

m.dat <- list(y=y, xvar1=x, J=length(y))
m.dat2 <- spmrf.get.data(m.dat)

nItr=2000; thin=1; nChains=3; nCores=3

m <- sampling(spmrf.normal.o2.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores) # control=list(adapt_delta=0.96, max_treedepth=12)

alpha <- 0.05
plow <- alpha/2
phigh <- 1 - alpha/2

theta <- rstan::extract(m, "theta")[[1]]
y.med <- apply(theta, 2, median)
y.mad <- apply(theta, 2, mad)
y.l <- apply(theta, 2, quantile, probs = plow)
y.u <- apply(theta, 2, quantile, probs = phigh)

x.uq <- unique(m.dat2$xvar1)
ylim <- range(c(y.l, y.u))

file <- paste0(formula.dir,"/IP_spmrf_normal_o2_model_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(3.75, 4.0, 0.5, 0.5), mgp=c(2.75,0.6,0),tcl = -0.3)
plot(x.uq, y.med, type = "n", las=1, ylim=ylim, xlab="Formula [%]", ylab="IP")
polygon(c(x.uq, rev(x.uq)), c(y.l, rev(y.u)),  border = NA, col='gray90')
lines(x.uq, y.med, col="blue", lwd=1)
par(op)
dev.off()
system(paste0("open ",file))


spmrf.norm.o2.fn <- function(x, y, nItr=2000, thin=1, nChains=3, nCores=3, alpha=0.05)
{
    o <- order(x)
    x <- x[o]
    y <- y[o]

    m.dat <- list(y=y, xvar1=x, J=length(y))
    m.dat2 <- spmrf.get.data(m.dat)

    r <- sampling(spmrf.normal.o2.model, data=m.dat2, iter=nItr, thin=thin, chains=nChains, cores=nCores)

    theta <- rstan::extract(r, "theta")[[1]]
    plow <- alpha/2
    phigh <- 1 - alpha/2
    y.med <- apply(theta, 2, median)
    y.mad <- apply(theta, 2, mad)
    y.l <- apply(theta, 2, quantile, probs = plow)
    y.u <- apply(theta, 2, quantile, probs = phigh)

    idx <- seq(y.med)
    r2 <- gEff.pval(x[idx], y.med, y.mad)

    x.uq <- unique(m.dat2$xvar1)
    ylim <- range(c(y.l, y.u))

    list(x=x,
         x.uq=x.uq,
         ylim=ylim,
         y=y,
         r=r,
         r2=r2,
         y.med=y.med,
         y.mad=y.mad,
         y.l=y.l,
         y.u=y.u)
}

ip <- st$IP
mm <- st$doseMBM
log10.fp <- log10(as.numeric(gsub("\\%","",st$FormulaPerc))+1)
log10.dbmp <- log10(as.numeric(gsub("\\%","",st$DBMperc))+1)
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
log10.fp <- log10.fp[idx]
log10.dbmp <- log10.dbmp[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
log10.fp <- log10.fp[idx]
log10.dbmp <- log10.dbmp[idx]



r <- spmrf.norm.o2.fn(mm,ip)

file <- paste0(formula.dir,"/IP_spmrf_normal_o2_model_plot.pdf")
pdf(file, width=6, height=6)

op <- par(mar=c(3.75, 4.0, 0.5, 0.5), mgp=c(2.75,0.6,0),tcl = -0.3)
plot(r$x.uq, r$y.med, type = "n", las=1, ylim=r$ylim, xlab="Formula [%]", ylab="IP")
polygon(c(r$x.uq, rev(r$x.uq)), c(r$y.l, rev(r$y.u)),  border = NA, col='gray90')
lines(r$x.uq, r$y.med, col="blue", lwd=1)
par(op)


dev.off()
system(paste0("open ",file))



file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_pFormula_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
idx <- log10.fp>0
points(mm[idx],ip[idx], pch=19, col='red', cex=log10.fp[idx]/2)
points(mm[!idx],ip[!idx], pch=19, col='black',cex=0.5)
par(op)
dev.off()
system(paste0("open ",file))



## PMAatDosing

ip <- st$IP
pma <- st$PMAatDosing
idx <- is.finite(ip) & is.finite(pma)
ip <- ip[idx]
pma <- pma[idx]

file <- paste0(mm.dir,"/ip_vs_PMA_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.25,0.6,0),tcl = -0.3)
plot.loess(pma,ip, xlab="Postmenstrual Age [weeks]", ylab="IP")
par(op)
dev.off()
system(paste0("open ",file))



ip <- st$IP
mm <- st$doseMBM
pma <- st$PMAatDosing
idx <- is.finite(ip) & is.finite(mm)
ip <- ip[idx]
mm <- mm[idx]
pma <- pma[idx]
idx <- mm<600
ip <- ip[idx]
mm <- mm[idx]
pma <- pma[idx]

myHist(pma)

pma1 <- pma-min(pma)
pma.scaled <- pma1/max(pma1)

myHist(pma.scaled)


file <- paste0(mm.dir,"/ip_vs_mothers_milk_with_PMA_loess_plot.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(4, 4, 0.5, 0.5), mgp=c(2.5,0.6,0),tcl = -0.3)
plot.loess(mm,ip, xlab="Cumulative mother's milk amount [g]", ylab="IP")
abline(v=130,col='gray')
points(mm,ip, pch=1, cex=(1-pma.scaled)+1)
par(op)
dev.off()
system(paste0("open ",file))

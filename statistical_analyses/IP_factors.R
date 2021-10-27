##
## What factors are associated with IP?
##

##
## to look into it without any assumptions on the nature of association,
## Hilbert Schmidt independence criterion (HSIC) statistic will be used
##

## save(subj.mt, file="subj_mt.rda")
load("subject_table.rda")


library(dHSIC)

## example
x <- matrix(rbinom(100,1,0.5),ncol=1)
y <- matrix(2*x[,1], ncol=1)
X <- list(x,y)
dhsic.test(X)

colnames(subj.mt)[4]

ip.dhsic.stat <- c()
ip.dhsic.pval <- c()
for ( i in seq(ncol(subj.mt)) )
{
    idx <- is.finite(subj.mt[,i])
    if ( is.numeric(subj.mt[idx,i]) && !(colnames(subj.mt)[i] %in% c("IP", "Ippattern", "No", "antsteroids")) )
    {
        idx <- is.finite(subj.mt$IP) & is.finite(subj.mt[,i]) & subj.mt[,i]>0
        if ( sum(idx) > 10 )
        {
            var.name <- colnames(subj.mt)[i]
            x <- matrix(subj.mt$IP[idx], ncol=1)
            y <- matrix(subj.mt[idx,i], ncol=1)
            r <- dhsic.test(x,y)
            ip.dhsic.stat[var.name] <- r$statistic
            ip.dhsic.pval[var.name] <- r$p.value
        }
    }
}


file <- "../pics/IP_vs_num_factors_dHSIC.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
plot(ip.dhsic.stat, log10(ip.dhsic.pval), ylab="log10(HSIC p-val)", xlab="HSIC statistic", las=1)
par(op)
dev.off()
file
system(paste0("open ",file))

names(ip.dhsic.stat)[ip.dhsic.stat>1]


ip.dhsic.qval <- p.adjust(ip.dhsic.pval, method="fdr")

file <- "../pics/IP_vs_num_factors_dHSIC_qvals.pdf"
pdf(file, width=6, height=6)
op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
plot(ip.dhsic.stat, log10(ip.dhsic.qval), ylab="log10(HSIC q-val)", xlab="HSIC statistic", las=1)
abline(h=log10(0.05), col='gray')
par(op)
dev.off()
file
system(paste0("open ",file))


thld <- 1
idx <- ip.dhsic.qval<thld
ns <- names(ip.dhsic.stat)[idx]
X <- cbind(ip.dhsic.stat[idx], ip.dhsic.pval[idx], ip.dhsic.qval[idx])
rownames(X) <- ns
colnames(X) <- c("HSIC statistic", "p-val","q-val")
o <- order(X[,1], decreasing = T)
X <- X[o,]
X

write.csv(X, file="IP_vs_num_factors_dHSIC.csv", quote=F)


myHist(log10(ip.dhsic.pval), xlab="log10(HSIC p-val)")


thld <- 0.05
(sel.factors <- names(ip.dhsic.stat)[ip.dhsic.qval<thld])


out.dir <- "../pics/dhsic_dir/"
dir.create(out.dir)

for ( sel in sel.factors )
{
    print(sel)
    idx <- is.finite(subj.mt$IP) & is.finite(subj.mt[,sel]) & subj.mt[,sel]>0
    y <- subj.mt$IP[idx]
    x <- subj.mt[idx,sel]
    file <- paste0(out.dir,"IP_vs_",sel,".pdf")
    pdf(file, width=6, height=6)
    op <- par(mar=c(3.5, 3.5, 0.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
    plot.loess(x, y, las=1, span=0.95,
               xlab=sel, ylab="IP")
    par(op)
    dev.off()
}



library(corrgram)

PMA.factors <- c("PMAatDosing", "PMAenrollment", "GA", "BW", "weightDosing")

file <- paste0(out.dir,"PMA_factors_pearson_correlogram.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(1.5, 1.5, 0.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
corrgram(subj.mt[,PMA.factors], lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)
par(op)
dev.off()
system(paste0("open ",file))


PMA.factors <- c("PMAatDosing", "PMAenrollment", "GA", "BW", "weightDosing")

file <- paste0(out.dir,"PMA_factors_pearson_correlogram.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(1.5, 1.5, 0.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
corrgram(subj.mt[,PMA.factors], lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)
par(op)
dev.off()
system(paste0("open ",file))



(milk.factors <- setdiff(sel.factors, PMA.factors))
## [1] "doseMBM"   "totalBM"   "totalMilk" "ratioBMF"

file <- paste0(out.dir,"milk_factors_pearson_correlogram.pdf")
pdf(file, width=6, height=6)
op <- par(mar=c(1.5, 1.5, 0.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
corrgram(subj.mt[,milk.factors], lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)
par(op)
dev.off()
system(paste0("open ",file))



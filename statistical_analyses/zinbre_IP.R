## ----------------------------------------------
## zero-inflated negative-binomial models
## of relative abundance vs log(LR-ration)
## ----------------------------------------------

## selecting phylotypes with absolute value of correlation greater or equal to
## 0.5 and the number of samples where it was detected at least 3.

pacman::p_load(MASS,FNN,mgcv,rms,lme4,nlme) # https://cran.r-project.org/web/packages/FNN/FNN.pdf # mgcv for GAMs
setwd("ZINBRE/")

library(rjags)
source("jags_utils.R") #load ci.mcmc(), plot.ci()
library(parallel)

jagsFn <- function(i) {
  ph <- selPh[i]
  j.dat <- list(y=as.vector(ct[,ph]),
                log.total=log.total,
                x=log(LR.rat),
                subjID=as.integer(factor(subjID)),
                nSubjIDs=length(unique(subjID)))
  str(j.dat)
  j.dat
  
  j.m <- jags.model("zinbin_withOffset_contExplVar_ri_hCauchy.bug", j.dat, n.chains=2)
  nIter <- 10000
  update(j.m, nIter)
  j.pars <- c ("a","b","c","d","e","sigma.eta")
  j.out <- coda.samples(j.m, j.pars, thin=10, n.iter=nIter)
  gelman.diag(j.out)
  (j.ci <- ci.mcmc(j.out))
  
  i.b <- 2
  i.e <- 5
  
  ## computing p-values
  j.mat <- as.matrix(j.out)
  j.b <- j.mat[,"b"]
  pval <- NA
  if ( median(j.b) < 0 ){
    pval <- 1 - pnorm(0, mean=median(j.b), sd=sd(j.b))
  } else {
    pval <- pnorm(0, mean=median(j.b), sd=sd(j.b))
  }
  
  list(c(j.ci[i.b,c(1,3:4)], pval),
       expit(j.ci[i.e,c(1,3:4)]))
}

#~~~~~ASV~~~~~~~~~
load("ASV.RData")
LR.rat <- mt$PMA
idx <- !is.na(LR.rat)
LR.rat <- LR.rat[idx]
ptC <- pt[,idx]

expit <- function(x) 1/(1 + exp(-x))
nProc <- 8
selPh=selPh.10
log.total <- as.vector(log(rowSums(ct)))
length(log.total) # 59
nrow(ct)          # 59

subjID <- mt[rownames(mt), "subjID"]
names(subjID) <- rownames(mt)
table(subjID)

res <- mclapply( seq(selPh), mc.cores = nProc, jagsFn)
#save(res, file="Rdata/perm_resIP_ASV.RData")

bCoef <- matrix(nrow=length(selPh), ncol=4)
colnames(bCoef) <- c("mean","2.5%","9.75%","pval")
rownames(bCoef) <- selPh
eCoef <- matrix(nrow=length(selPh), ncol=3)
colnames(eCoef) <- c("mean","2.5%","9.75%")
rownames(eCoef) <- selPh

for ( i in seq(selPh) ) {
  print(i)
  #x <- jagsFn(i)
  x <- res[[i]]
  bCoef[i,] <- x[[1]]
  eCoef[i,] <- x[[2]]
}

## looking for significant resluts
(i <- which( bCoef[,2] * bCoef[,3]>0 ))

bCoef[i,]
qvals <- p.adjust(bCoef[,4], method='fdr')
length(qvals[qvals<0.1]) # 0
min(qvals) # [1] 0.1592211
which.min(qvals)
bCoef <- cbind(bCoef, qvals)
o <- order(qvals)
bCoef <- bCoef[o,]

write.csv(file="result/w2t_PMA_ASV.csv", bCoef)


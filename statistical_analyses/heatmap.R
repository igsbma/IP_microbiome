##heatmap

mytotal=read.table("ASV.txt",header=T,row.names=1,check.names=F,sep="\t")

dim(mytotal)
head(mytotal)
mytotal[1:5,1:5]
mydata=mytotal
mydata[1:5,1:5]

IPf=factor(mytotal$category)
levels(IPf)
IPTbl=c()
IPTbl[levels(IPf)[1]]="red"
IPTbl[levels(IPf)[2]]="green"

sitef=factor(mytotal$PMA)
levels(sitef)
siteTbl=c()
siteTbl[levels(sitef)[1]]="navy"
siteTbl[levels(sitef)[2]]="orange"
siteTbl[levels(sitef)[3]]="deeppink"

cs=colSums(mydata)
o=order(cs,decreasing=T)
ct=mydata[,o]
ct[1:5,1:5]
pt=t(apply(ct,1,function(x) 100*x/sum(x)))
pt[1:5,1:5]
nCols=37
dim(pt)
hc.ward=hclust(dist(pt),method="ward.D")

nClrs=6
memb.ward=cutree(hc.ward,k=nClrs)
table(memb.ward)
write.csv(memb.ward, file = "ward_sample_cluster.csv")

lx=0.86
sideBars.ward=cbind(colorTbl[memb.ward],IPTbl[IPf],siteTbl[sitef])
colnames(sideBars.ward)=c("Cluster","IP","PMA")
rownames(sideBars.ward)=memb.ward

postscript("heatmap_ward.eps",hor=T)

heatmap2(as.matrix(pt)[,1:nCols],         # show only nCols first columns of the table tbl
         col=rainbow(50,start=1/6,end=0),
         #col=greenred(100), # set heatmap color palette
         #col=colorRampPalette(c("darkblue","black","yellow"))(100),
         Colv=NA,                         # comment this out if you want to have columns clustered as well
         #Rowv = NA,
         Rowv=as.dendrogram(hc.ward),  
         #Rowv=as.dendrogram(hc.ward),    # use our h. clustering
         RowSideColors=sideBars.ward,     # RowSideColors=sideBars,          # this sets three color columns: clustering, pH and Nugent score
         RowSideTitleCex=0.5,             # scaling factor for the titles of the color bars
         RowSideTitleLine=0.5,            # the distance of the titles for the color bars
         margins=c(18,4),                 # =c(bottom margin, right margin)
         #labRow=rownames(mydata),         # add row labels
         labRow=NA,                       # suppress row labels
         cexCol=0.8,                      # magnification factor for column labels
         xlas=2,                          # column labels are to be perpendicular to the x-axis
         main="")                         # main title

legend(lx,0.97,legend=c("high","low"),bg="white",fill=IPTbl,title="La/Rh ratio",cex=0.7)
legend(lx,0.87,legend=c("early","late","very early"),bg="white",fill=siteTbl,title="PMA",cex=0.7)

dev.off()

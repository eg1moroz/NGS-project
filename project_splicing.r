# 3.4 AS ######

library(SAJR)
as = readRDS('Rdata/as.Rdata')
length(as)
as[1:2,]
table(as$seg$sites,as$seg$position)
hist(log10(1+as$i+as$e))
as$ir[as$i+as$e<10] = NA
as = as[as$seg$sites=='ad' & as$seg$type=='ALT',] # focus on cassette exons
table(apply(is.na(as$ir),1,sum))
as = as[apply(is.na(as$ir),1,sum) < 2 ,]
length(as)

mds = cmdscale(1-cor(as$ir,u='p'),k=2)
plot(mds, pch=19,col=meta$col,cex=meta$nage*3+0.6,xlab='Dim 1',ylab='Dim 2',main='Spearman'   ,bty='n')

meta
?fitSAGLM

as.pv = fitSAGLM(as,x ~ region,meta,return.pv = TRUE)
as.pv[1:10,-1]

hist(as.pv[,2])
apply(as.pv<0.05,2,sum,na.rm=T)

as.qv = apply(as.pv,2,p.adjust,m='BH')
as.qv[1:10,]

apply(as.qv<0.05,2,sum,na.rm=T)

# lets check most significant exon
as[order(as.qv[,'region'])[1:10],]
meta
jpeg(fil = "plotASit.jpg")
plot(as$ir['ENSG00000077380.s11',],pch=19,col=meta$col)
dev.off()
# 3.5 Visualization #####
BiocManager::install('GenomicAlignments')

# ������������ ���������� - ��� BAM ������

# 3.6 GO genes with dif AS


library(topGO)
BiocManager::install("org.Hs.eg.db")
str(as.qv)

rownames(as.qv) = substr(rownames(as.qv),1,15)
as.integer(apply(as.qv<0.05,1,sum)>0)
length(apply(as.qv<0.05,1,sum))
s = factor(as.integer(apply(as.qv<0.05,1,sum)>0))
names(s) = rownames(as.qv)
tgo1 = new("topGOdata", ontology = "BP",
           allGenes = s,
           nodeSize = 6,
           annotationFun = annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl') 

### tgo1 - NA significant genes (?)
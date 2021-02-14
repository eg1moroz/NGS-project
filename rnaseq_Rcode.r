
# 2.0 count reads in R ###

library(Rsubread)

counts = featureCounts(list.files('project/bam1','*.bam$',full.names = T),annot.ext = 'project/index/Homo_sapiens.GRCh38.102.chr.gtf', isGTFAnnotationFile=TRUE, isPairedEnd = TRUE, nthreads = 8)

###avarage asigned aligments - 45% (??)

saveRDS(counts,'project/Rdata/counts.Rdata')
counts = readRDS('project/Rdata/counts.Rdata')

counts = counts$counts

colnames(counts) = sub('.bam','',colnames(counts))
for (i in 1:4) {
  colnames(counts) = sub(colnames(counts)[i],paste('Cwhite',as.character(i), sep ='_'),colnames(counts))
}  # change sample name to CWhiteX - Cerebellar White Matter
for (i in 5:8) {
  colnames(counts) = sub(colnames(counts)[i],paste('Prefro',as.character(i), sep ='_'),colnames(counts))
} 
counts[1:2,] 

meta = data.frame(region=substr(colnames(counts),1,6))
rownames(meta) = colnames(counts)
meta$col = c(Cwhite='blue',Prefro='red')[meta$tissue]

# 2.1 filter by coverage

table(apply(counts,1,sum) == 0)
table(apply(counts,1,mean) >= 10)
counts = counts[apply(counts,1,mean) >= 10,]

# 2.2 quality control

hist(log(counts[,1]))
jpeg(fil = "project/Rdata/corplot.jpg")
plot(counts[,1],counts[,2],log='xy')
dev.off()


cors = cor(counts,m='sp')
jpeg(fil = "project/Rdata/heatmapspir.jpg")
heatmap(1-cors,symm = T,distfun = function(x){as.dist(x)},ColSideColors = meta$col)
dev.off()

# 3.0 Dif Exp ###

library(edgeR)

groups = c(1,1,1,1,2,2,2,2)
edger = DGEList(counts, group = groups)
edger = calcNormFactors(edger,method='RLE')
edger$samples
design = model.matrix(~ region,data = meta)
edger = estimateDisp(edger,design)

# find DE genes, by two different ways

glm = glmFit(edger,design)
exat = exactTest(edger)

pv = cbind(tissue= glmLRT(glm,2)$table$PValue)
pv2 = cbind(exat$table$PValue)

rownames(pv) = rownames(counts)
rownames(pv2) = rownames(counts)

hist(pv[,1])
hist(pv2[,1])

qv = apply(pv,2,p.adjust,m='BH')
qv2 = apply(pv2, 2, p.adjust,m='BH')

hist(qv[,1])
hist(qv2[,1])
apply(qv < 0.05,2,sum)
apply(qv2 < 0.05,2,sum)


# 3.3 GO ######

library(topGO)
BiocManager::install("org.Hs.eg.db")

s = factor(as.integer(apply(qv<0.05,1,sum)>0))
names(s) = rownames(qv)
tgo1 = new("topGOdata", ontology = "BP",
           allGenes = s,
           nodeSize = 8,
           annotationFun = annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl') ### nodesize? Если категория представлена меньше 8-ми генов - не считать её. Но это же уменьшает точноть?

resultFish = runTest(tgo1, algorithm = "classic", statistic = "fisher")
hist(score(resultFish))
qvGO = p.adjust(score(resultFish),m='BH')
hist(qvGO)
sort(qvGO)[1:10]
sort(score(resultFish))[1:10]
GenTable(tgo1,resultFish,topNodes=10)
dev.off()
showSigOfNodes(tgo1, score(resultFish), firstSigNodes = 10, useInfo ='all')

#3.4 Clustering ####

cpm = cpm(edger)
cpm.s = cpm[apply(qv<0.05,1,sum)>0,]
dim(cpm.s)
cpm.s = t(scale(t(cpm.s)))

hcl = hclust(as.dist(1-cor(t(cpm.s))))
plot(hcl)
cl = cutree(hcl,2)
table(cl)

renameClustsBySize = function(x){
  t = table(x)
  n = names(t)
  o = order(t,decreasing=T)
  r = x
  for(i in 1:length(o))
    r[x == n[o[i]]] = i
  r
}
cl = renameClustsBySize(cl)
table(cl)

      
means1 = apply(cpm.s[cl==1,,drop=F],2,mean) # среднее значение каждого образца в кластере 1
means2 = apply(cpm.s[cl==2,,drop=F],2,mean) 

plot(means1, pch=19, col=meta$col,main=paste0('c',1,' (',sum(cl==1),')'), cex=3) #график уровня экспрессиии между образцами в одном кластере
plot(means2, pch=19, col=meta$col,main=paste0('c',2,' (',sum(cl==2),')'), cex=3)

#GO of clusters


library(topGO)
BiocManager::install("org.Hs.eg.db")

s = factor(as.integer(apply(matrix(cl),1,sum)>1)) #сравнение GO категорий группы. Значимой группой является вторая - понижение экспрессии
names(s) = rownames(qv1)
tgo1 = new("topGOdata", ontology = "BP",
           allGenes = s,
           nodeSize = 8,
           annotationFun = annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl') 

resultFish = runTest(tgo1, algorithm= "classic", statistic = "fisher")
hist(score(resultFish))
dev.off()
qvGO = p.adjust(score(resultFish),m='BH')
hist(qvGO)
sum(qvGO<0.05)
sort(qvGO)[1:10]
resultFish
sort(score(resultFish))[1:10]

allRes <- GenTable(tgo1,resultFish,topNodes=10)
write.xlsx(GenTable(tgo1,resultFish,topNodes=20), file = 'ClDownGO.xls')

printGraph(tgo1, resultFish, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

s = factor(as.integer(apply(matrix(cl),1,sum)<2)) #сравнение GO категорий группы. Значимой группой является первая - повыщение
names(s) = rownames(qv1)
tgo1 = new("topGOdata", ontology = "BP",
           allGenes = s,
           nodeSize = 8,
           annotationFun = annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl') 

resultFish = runTest(tgo1, algorithm= "classic", statistic = "fisher")
hist(score(resultFish))
dev.off()
qvGO = p.adjust(score(resultFish),m='BH')
hist(qvGO)
sum(qvGO<0.05)
sort(qvGO)[1:10]

sort(score(resultFish))[1:10]

allRes <- GenTable(tgo1,resultFish,topNodes=10)
write.xlsx(GenTable(tgo1,resultFish,topNodes=20), file = 'ClDownGO.xls')
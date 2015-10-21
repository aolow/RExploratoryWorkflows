# Data Exploration Oncogene Project
# 31st July 2015
# A. Olow

setwd("/Users/aolow/Box\ Sync/Documents/VantveerLab/Moasser_project/Analysis")

# Load data (M1E samples swapped - corrected version here)

load("/Users/aolow/Box\ Sync/Documents/VantveerLab/Moasser_project/Analysis/Cdata.RData")
load("/Users/aolow/Box\ Sync/Documents/VantveerLab/Moasser_project/Analysis/Mdata.RData")

# combine the two datasets

Cdat <- Cdata[,11:33]
row.names(Cdat) <- Cdata[,2]

Mdat <- Mdata[,10:33]
row.names(Mdat) <- Mdata[,1]

# keep only the common probes
cmb <- intersect(row.names(Cdat), row.names(Mdat))

Cdat2 <- subset(Cdat, row.names(Cdat) %in% cmb)
Mdat2 <- subset(Mdat, row.names(Mdat) %in% cmb)

# merge datasets
Cdat2$name <- row.names(Cdat2)
Mdat2$name <- row.names(Mdat2)
names(Mdat2)[19:24] <- sapply(names(Mdat2)[19:24], function(x) paste0("M_", x))

Data <- merge(Cdat2, Mdat2, by="name")
rownames(Data) <- Data$name
Data <- Data[,-1]

### Make matrix

DataM <- data.matrix(Data)
#save(DataM, file="./DataM.RData")
#load(file="./DataM.RData")

library(limma)

nDataM <- normalizeQuantiles(DataM, ties=TRUE)
#save(nDataM, file="./nDataM.RData")
#load(file="./nDataM.RData")

# Visualize pre-after normalization
par(mfrow = c(1,2))
par(oma=c(8,0,0,0))

boxplot(DataM, las=2, cex = 0.1)
boxplot(nDataM, las=2, cex = 0.1)

dev.off()

# make into an expression dataset
library(Biobase)

samples <- as.character(c(rep("10A-NeuT",5), rep("C-1", 6), rep("C-2",6),
                          rep("C-5",6), rep("M1-A",6), rep("M1-B", 6), rep("M1-E",6),
                          rep("M_MCF10A-NeuT",6)))

culture <- as.character(c(rep("P", 5), rep("C", 18), rep("M", 18), rep("P", 6)))
                        
                        
treatment <- as.character(c(rep("DMSO",3), "LAP", "LAP", 
                            rep("DMSO",3), rep("LAP", 3),
                            rep("DMSO",3), rep("LAP", 3),
                            rep("DMSO",3), rep("LAP", 3),
                            rep("DMSO",3), rep("LAP", 3),
                            rep("DMSO",3), rep("LAP", 3),
                            rep("DMSO",3), rep("LAP", 3),
                            rep("DMSO",3), rep("LAP", 3)))

addiction <- as.character(c(rep("NO",5), rep("YES",6), rep("NO",12),
                          rep("YES",18), rep("NO",6)))

pData <- cbind(samples,treatment,addiction, culture)
rownames(pData) <- colnames(nDataM)
pData <- as.data.frame(pData)

metadata <- data.frame(labelDescription=
                         c("Sample origin",
                           "Treatment",
                           "Oncogene Addiction Status Yes/No",
                           "Cell Culture Type In vitro only (C) / harvested from mice (M)"),
                       row.names=c("samples", "treatment", "addiction", "culture"))

phenoData <- new("AnnotatedDataFrame",
                 data=pData, varMetadata=metadata)

experimentData <- new("MIAME",
                      name = "M. Moasser",
                      lab = "Moasser Lab",
                      title = "Oncogene Addiction Project")

eset <- ExpressionSet(assayData=nDataM,
                       phenoData=phenoData,
                       experimentData=experimentData, 
                      annotation="lumiHumanAll.db")

#save(eset, file="./eset.RData")
load(file="./eset.RData")

### exploring expression dataset

table(pData(eset)$treatment, pData(eset)$addiction)

# exploring the data
par(mfrow=c(1,1))
par(oma=c(8,0,0,0))
boxplot(exprs(eset), col=as.numeric(pData(eset)$addiction), las=2, cex=0.1)
dev.off()

# MA plot
library(data.table)

Index <- as.numeric(pData(eset)$addiction)
d <- rowMeans(exprs(eset)[,Index==2]) - rowMeans(exprs(eset)[,Index==1])
a <- rowMeans(exprs(eset))

smoothScatter(a,d, main="MA plot", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")

################
# PCA

# not scaling as may not make sense to scale expression data
tst <- prcomp(t(nDataM))

palette(c("blue", "turquoise"))
plot(tst$x[,1], tst$x[,2], pch = 20 , cex = 2, col=as.double(as.factor(treatment)),
     xlab = "PC1", ylab="PC2", main="PCA", ylim=c(-18,40))
palette(c("black", "red"))
text(tst$x[,1]+0.1, tst$x[,2]+0.6, samples, cex=.8, col=(as.double(as.factor(addiction))))
legend(-44,40,
       c("Non-Addicted","Addicted"),
       cex=0.8,
       text.col=c("black","red")) 
legend(-32,40,
       c("DMSO","LAP"),
       cex=0.8,
       pch=20,
       col=c("blue","turquoise")) 
dev.off()

#MDS plot (PCA)
y <- exprs(eset)
plotMDS(y)

####
######
#########

data <- nDataM[grep("LAP", colnames(nDataM)),]
labs <- samples[grep("LAP", colnames(nDataM))]

table(labs)

# PCA

pr.out <- prcomp(data, scale=T)

#assign distinct color to each element of numeric vector (each of the 64 cell lines, based on cancer type)

Cols <- function(vec){
  cols <- rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

# plot principal component score vectors
par(mfrow=c(1,2))
plot(pr.out$x[,1:2],
     col=Cols(labs), 
     pch=19,
     xlab="Z1",
     ylab="Z2")

plot(pr.out$x[,c(1,3)],
     col=Cols(labs), 
     pch=19,
     xlab="Z1",
     ylab="Z3")

# proportion of variance explained

summary(pr.out)
plot(pr.out)

# pve for each principal component (scree plots)

pve <- 100* pr.out$sdev^2/sum(pr.out$sdev^2)
par(mfrow=c(1,2))

plot(pve, type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="brown3")

# also you can find pve:
plot(summary(pr.out)$importance[2,], type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(summary(pr.out)$importance[3,],  type="o", ylab="Cumulative PVE", xlab="Principal Component", col="brown3")

#### Clustering observations of the DMSO data

sd.data <- scale(data)
par(mfrow=c(1,3))

data.dist <- dist(sd.data)

plot(hclust(data.dist), labels = labs,
     main = "Complete Linkage",
     xlab = "", ylab = "", sub = "")

plot(hclust(data.dist, method = "average"), labels = labs,
     main = "Average Linkage",
     xlab = "", ylab = "", sub = "")

plot(hclust(data.dist, method = "single"), labels = labs,
     main = "Single Linkage",
     xlab = "", ylab = "", sub = "")

# cut into clusters
hc.out <- hclust(dist(sd.data))
hc.clusters <- cutree(hc.out,4)
table(hc.clusters, labs)

# plot the cut dendrogram
par(mfrow=c(1,1))
plot(hc.out, labels = labs)
abline(h=10, col="red")

hc.out

set.seed(2)
km.out <- kmeans(sd.data, 4, nstart=20)
km.clusters <- km.out$cluster
table(km.clusters, hc.clusters)

# hierarchical clustering on the first few principal component score vectors

hc.out <- hclust(dist(pr.out$x[,1:5]))

plot(hc.out, labels = labs, main = "Hier. Clust. DMSO on First Five Score Vectors")
table(cutree(hc.out, 4), labs)






# within cell line origin correlation

# samples from same culture are not independent; estimating within-culture correlation
ct <- factor(eset$culture)
design <- model.matrix(~0+ct)
colnames(design) <- levels(ct)
dupcor <- duplicateCorrelation(y, design, block=eset$samples)
dupcor$consensus.correlation


# DE genes
library(limma)
y <- exprs(eset)

design <- model.matrix(factor(eset$treatment)~factor(eset$addiction))
fit <- lmFit(y, design
             #, 
             #block=eset$samples, 
             #correlation=dupcor$consensus.correlation
             )

ebayes <- eBayes(fit)

tab <- topTable(ebayes, coef=2, adjust="fdr", n=150)
tab

### heatmap

par(mar=c(5,4,4,2)+0.1, 
    oma = c(7, 1,1,1),
    cex.lab=0.01)


# clustered just per sample variance
library(gplots)
evar <- apply(exprs(eset),1,var)
evar_ordr <- order(evar, decreasing = T)
inds <- head(evar_ordr, length(evar)*.1)
heatmap(exprs(eset)[inds, ], col = redgreen(75))


evar <- apply(nDataM,1,var)
evar_ordr <- order(evar, decreasing = T)
inds <- head(evar_ordr, length(evar)*.1)
heatmap(nDataM[inds, ], col = redgreen(75))

library(gplots)

# add color bar for addiction status
color.map <- function(addiction) { if (addiction=="YES") "#FF0000" else "#0000FF" }
addict_colors <- unlist(lapply(eset$addiction, color.map))

heatmap(nDataM[inds, ], col = redgreen(75), 
        ColSideColors=addict_colors)

heatmap.2(nDataM[inds,], col = redgreen(75), scale="row")

############
## explore clustering 
inds <- head(evar_ordr, length(evar)*.1)
d <- dist(x = t(nDataM[inds, ]), method = 'manhattan')

par(mfrow = c(2, 2))
plot(hclust(d =  d, method = 'ward.D2'))
plot(hclust(d =  d, method = 'average'))
plot(hclust(d =  d, method = 'single'))
plot(hclust(d =  d, method = 'complete'))
# 
# d <- dist(x = t(nDataM[inds, ]), method = 'euclidean')
# 
# plot(hclust(d =  d, method = 'ward'))
# plot(hclust(d =  d, method = 'average'))
# plot(hclust(d =  d, method = 'single'))
# plot(hclust(d =  d, method = 'complete'))
# 
# d <- dist(x = t(nDataM[inds, ]), method = 'maximum')
# 
# plot(hclust(d =  d, method = 'ward'))
# plot(hclust(d =  d, method = 'average'))
# plot(hclust(d =  d, method = 'single'))
# plot(hclust(d =  d, method = 'complete'))


cor_mat <- cor(nDataM[inds, ], method = 'pearson')
d <- as.dist(1 - cor_mat^2)

par(mfrow = c(2, 2))
plot(hclust(d =  d, method = 'ward.D'))
plot(hclust(d =  d, method = 'average'))
plot(hclust(d =  d, method = 'single'))
plot(hclust(d =  d, method = 'complete'))
dev.off()

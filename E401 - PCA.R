### PCA for E401 -- Andrew R Gross -- 11JUN21
### INPUT: Expression data; 
### OUTPUT: PCA Plots

############################################################################################
### Header
library(ggplot2)
library(ggbiplot)

############################################################################################
### Functions

### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iEC cell pellets and CM for Proteomics core/Proteomics results/')
list.files()
iec.data <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/2021_64_DataReport_Optra_PanHuman_mapDIA.csv')
iec.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/E401-metadata.csv')

############################################################################################
### Format

### Rename columns

### Filter low expression genes
summary(iec.data)
iec.data.max <- apply(iec.data, 1, max)
rows.to.keep <- which(iec.data.max >10)
length(rows.to.keep)
### Reorder columns
results.iec <- iec.data[rows.to.keep,]

### Replace row names
#results.als <- convert.ids(results.als)

### Convert to matrix
results.iec <- as.matrix(results.iec)

############################################################################################
### Calculate Principle Components

### Calculate the actual components
pca.iec <- prcomp(t(results.iec), scale = TRUE)

### Calculate the percent variation accounted for by each component
pca.data.var <- pca.iec$sdev^2
pca.data.var.per <- round(pca.data.var/sum(pca.data.var)*100, 1)

### Identify the genes with the largest influence
# PC1
l.score.pc1 <- pca.iec$rotation[,1]
l.score.pc1.ranked <- sort(abs(l.score.pc1), decreasing = TRUE)
l.score.pc1[names(l.score.pc1.ranked)][1:10]

# PC2
l.score.pc2 <- pca.iec$rotation[,2]
l.score.pc2.ranked <- sort(abs(l.score.pc2), decreasing = TRUE)
l.score.pc2[names(l.score.pc2.ranked)][1:10]

############################################################################################
### Plot
subtitle = ''

plot(pca.iec$x[,1], pca.iec$x[,2])
plot(pca.iec$x[,2], pca.iec$x[,3])

barplot(pca.data.var.per, main = 'Scree Plot', xlab = 'Principle Component', ylab = 'Percent Variation')

### GGPlot

pca.data.to.plot <- data.frame(Sample = rownames(pca.iec$x), 
                              PC1 = pca.iec$x[,1],
                              PC2 = pca.iec$x[,2],
                              PC3 = pca.iec$x[,3],
                              PC4 = pca.iec$x[,4])

(pca.data.to.plot <- cbind(pca.data.to.plot, iec.metadata))

pca.plot <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC2, label = Group, color = Group)) +
  geom_point() + 
  geom_text() + 
  xlab(paste('PC1 - ', pca.data.var.per[1], '%', sep = '')) +
  ylab(paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
  xlim(c(-50,70)) +
  theme_bw() +
  ggtitle('PCA of E401')

pca.plot








ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC2, label = ShortName)) +
  geom_point(size = 4, aes(color = Group)) +
  geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_color_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
  labs(title="PCA of RNA seq Expression", 
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) 

### PC1 v PC3
ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC3, label = ShortName)) +
  geom_point(size = 4, aes(color = Group)) +
  geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_color_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
  labs(title="PCA of RNA seq Expression", 
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC3 - ', pca.data.var.per[3], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) 




### Alternative coloration


ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC2, label = Group)) +
  geom_point(size = 4, aes(color = Group)) +
  geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_color_manual(values = c("#fab9b6", "#a6acf7", "#fe1c1c", "#000dc4")) +
  labs(title="PCA of RNA seq Expression", 
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12), 
        legend.position="top", legend.title = element_blank(), legend.key=element_blank(),
        legend.spacing.x = unit(5, "mm")) 



### PC2 v PC3
ggplot(data = pca.data.to.plot, aes(x = PC2, y = PC3, label = Sample)) +
  geom_text() + 
  theme_bw() +
  labs(title="PCA of E283", 
       subtitle = subtitle,
       x = paste('PC2 - ', pca.data.var.per[2], '%', sep = ''), 
       y = paste('PC3 - ', pca.data.var.per[3], '%', sep = ''))


### PC1 v PC3
ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC3, label = Sample)) +
  geom_text() + 
  theme_bw() +
  labs(title="PCA of E283", 
       subtitle = subtitle,
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC3 - ', pca.data.var.per[3], '%', sep = ''))


###############################################################################################
### Find the loadings that MATTER
### Convert all loadings into vectors of the first two components.  Find the axis of interest.
### Then, identify the loadings with the greatest magnitudes along this axis.

##############################################
### Find the slope and magnitude of all loadings
loadings <- data.frame(l.score.pc1, l.score.pc2)
  loadings <- loadings*100
loadings$slope = loadings[,2]/loadings[,1]
loadings$degree = atan(loadings$slope)*180/pi
loadings$magnitude = abs((loadings[,1]^2+loadings[,2]^2)^0.5)
loadings$mag2 = abs(loadings[,1])+abs(loadings[,2])

### Check the distribution of loadings
plot(dist(loadings$magnitude[sample(1:nrow(loadings),100,replace = FALSE)]))
plot(dist(loadings$mag2[sample(1:nrow(loadings),100,replace = FALSE)]))

subsample <- loadings[sample(1:nrow(loadings),1000,replace = FALSE),]
plot(subsample$l.score.pc1, subsample$magnitude)
plot(subsample$l.score.pc2, subsample$magnitude)

plot(subsample$l.score.pc1, subsample$mag2)
plot(subsample$l.score.pc2, subsample$mag2)
hist(loadings$magnitude)
hist(loadings$mag2)

plot(subsample$degree, subsample$magnitude)
plot(subsample$degree, subsample$mag2)

subsample$mag3 <- (subsample$magnitude + subsample$mag2*0.5)
plot(subsample$degree, subsample$mag3)

##############################################
### Find the axis of interest
### Find center point of control samples:
mean.x.of.ctr <- mean(pca.iec.to.plot$X[1:4])
mean.y.of.ctr <- mean(pca.iec.to.plot$Y[1:4])
mean.ctr <- c('CTR mean', mean.x.of.ctr, mean.y.of.ctr)

mean.x.of.als <- mean(pca.iec.to.plot$X[5:8])
mean.y.of.als <- mean(pca.iec.to.plot$Y[5:8])
mean.als <- c('ALS mean', mean.x.of.als, mean.y.of.als)

#means.x <- c(mean.x.of.ctr, mean.x.of.als)
#means.y <- c(mean.y.of.ctr, mean.y.of.als)
#centers.of.sample.groups <- data.frame(means.x, means.y)

pca.iec.to.plot <- data.frame(Sample = c(rownames(pca.iec$x),'CTR mean','ALS mean'), 
                              X = c(pca.iec$x[,1],mean.x.of.ctr, mean.x.of.als),
                              Y = c(pca.iec$x[,2],mean.y.of.ctr, mean.y.of.als),
                              color = c('Pink','Pink','Pink','Pink','Blue','Blue','Blue','Blue','Red','Dark_blue'))

### Plot with means of each group
pca.plot <- ggplot(data = pca.iec.to.plot, aes(x = X, y = Y, label = Sample)) +
    geom_text(aes(color = color)) + 
    #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","Red")) +
    scale_color_manual(values = c('Blue','Dark Green','Red','Orange')) +
    xlab(paste('PC1 - ', pca.data.var.per[1], '%', sep = '')) +
    ylab(paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
    theme_bw() + 
    ggtitle('PCA of E099') 

pca.plot


### Identify slope of interest
(slope.of.interest <- mean.y.of.als/mean.x.of.als)

### Visualize range of interest
range.1 = 1.2
range.2 = 2
range.3 = 10
pca.plot.w.slopes <- pca.plot + 
  geom_abline(slope=slope.of.interest*range.1, color = 'Gray20')+ geom_abline(slope=slope.of.interest/range.1, color = 'Gray20')+
  geom_abline(slope=slope.of.interest*range.2, color = 'Gray50')+ geom_abline(slope=slope.of.interest/range.2, color = 'Gray50')+
  geom_abline(slope=slope.of.interest*range.3, color = 'Gray70')+ geom_abline(slope=slope.of.interest/range.3, color = 'Gray70')
# geom_abline(slope = slope.of.interest)
pca.plot.w.slopes

##############################################
### Downselect loadings based on angle
rows.selected <- c()
### SET 1: Range 3
rows.of.interest <- intersect(which(loadings$slope < slope.of.interest/range.3),
                              which(loadings$slope > slope.of.interest*range.3))
#length(rows.of.interest)/nrow(loadings)
length(rows.of.interest)
loadings.of.interest <- loadings[rows.of.interest,]

subsample <- loadings.of.interest[sample(1:nrow(loadings),1000,replace = FALSE),]
plot(subsample$degree, subsample$magnitude)

### Find the nth percentile of magnitude
#(percentile.98 <- quantile(loadings.of.interest$magnitude, 0.9))
#rows.to.keep <- row.names(loadings.of.interest[which(loadings.of.interest$magnitude>percentile.98),])
#length(rows.to.keep)
#rows.selected <- c(rows.selected, rows.to.keep)
#length(rows.selected)
#loadings.of.interest <- loadings.of.interest[rows.to.keep,]

### SET 2: Range 2
rows.of.interest <- intersect(which(loadings$slope < slope.of.interest/range.2),
                              which(loadings$slope > slope.of.interest*range.2))
length(rows.of.interest)
loadings.of.interest <- loadings[rows.of.interest,]

### Find the nth percentile of magnitude
#(percentile.98 <- quantile(loadings.of.interest$magnitude, 0.99))
#rows.to.keep <- which(loadings.of.interest$magnitude>percentile.98)
#length(rows.to.keep)
#rows.selected <- c(rows.selected, rows.to.keep)
#length(rows.selected)

### SET 3: Range 1
rows.of.interest <- intersect(which(loadings$slope < slope.of.interest/range.1),
                              which(loadings$slope > slope.of.interest*range.1))
length(rows.of.interest)
loadings.of.interest <- loadings[rows.of.interest,]
### Find the nth percentile of magnitude
#(percentile.98 <- quantile(loadings.of.interest$magnitude, 0.98))
#rows.to.keep <- which(loadings.of.interest$magnitude>percentile.98)
#length(rows.to.keep)
#rows.selected <- c(rows.selected, rows.to.keep)
#length(rows.selected)

##############################################################################
#genes.to.plot <- loadings[1:2][unique(rows.selected),]*30
genes.to.plot <- loadings.of.interest
genes.to.plot[1:2] <- genes.to.plot[1:2]*30
genes.to.plot$x0 = 0
genes.to.plot$y0 = 0
genes.to.plot$Sample = 1
### Filter by level
genes.to.plot <- genes.to.plot[which(genes.to.plot$magnitude>1.4),]

#ggplot(data = genes.to.plot, aes(x = x0, y = y0, xend = l.score.pc1, yend = l.score.pc2, color = mag2)) +
#  geom_segment() +
#  scale_color_gradient(low = "white", high = "black")

pca.plot.w.slopes + geom_segment(data = genes.to.plot, aes(x = x0, y = y0, xend = l.score.pc1, yend = l.score.pc2))


##############################################################################
### Format row names
### Convert row names into gene names for easy reading or Entrez ids for GAGE analysis

loadings.null <- loadings[sample(1:nrow(loadings.of.interest),replace = FALSE),]

loadings.of.interest <- convert.ids(loadings.of.interest)
tpm.of.interest <- tpm.als[row.names(loadings.of.interest),]
tpm.of.interest <- convert.ids(tpm.of.interest)

tpm.o.i.ez <- convert.to.entrez(tpm.of.interest)
str(tpm.o.i.ez)
tpm.o.i.ez[4]
tpm.o.i.ez2 <- tpm.o.i.ez[[1]]
grep(8623, tpm.o.i.ez2$join)
tpm.o.i.ez2 <- tpm.o.i.ez2[-c(2202,2203),]
tpm.o.i.ez2 <- tpm.o.i.ez2[-581,]

### Assign new IDs to row names
row.names(tpm.o.i.ez2) <- tpm.o.i.ez2$join



loadings.null <- convert.ids(loadings.null)
tpm.null <- tpm.als[row.names(loadings.null),]
tpm.null <-convert.ids(tpm.null)
tpm.null <- convert.to.entrez(tpm.null)
tpm.null2 <- tpm.null[[1]]
tpm.null2 <- tpm.null2[-c(),]











### Filter TPM list based on loadings of interest

for.gage.loadings <- tpm.als[]

pway.loadings <- gage(for.gage.loadings, gsets = kegg.gs, ref = ctrl.index, samp = ko.index, compare = "unpaired")





### Annotate!

add.description(loadings.of.interest, 'external_gene_name')



load.scores.als <- pca.iec$rotation[,1]
gene.scores <- abs(load.scores.als)
gene.scores.ranked <- sort(gene.scores, decreasing = TRUE)
gene.scores.ranked <- names(gene.scores.ranked[1:10])
pca.iec$rotation[gene.scores.ranked,1]


setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/DEG analyses/')
row.num.pos = 10
gene.name <- results.als$Gene[row.num.pos]
data.als <- results.als[row.num.pos,][-c(1,2,3,12,13)]
data.ko <- results.ko[row.num.pos,][-c(1,2,3,10,11)]
expression <- t(data.als)
expression.ko <- t(data.ko)

expression <- rbind(expression,expression.ko)
Disease <- c("CTR", "CTR", "CTR", "CTR", "ALS", "ALS", "ALS", "ALS", 'WT', 'WT', 'WT', 'KO', 'KO', 'KO')

expression <- data.frame(expression, Disease)
names(expression) <- c('tpm', 'dis')

expression$dis <- factor(expression$dis, c('CTR','ALS','WT','KO'))

boxplot(tpm~dis, data=expression, main=gene.name, xlab="Genome type", ylab="Expression [TPM]")

png(paste0(gene.name,'.png'))
boxplot(tpm~dis, data=expression, main=gene.name, xlab="Genome type", ylab="Expression [TPM]")
dev.off()


p <- ggplot(test, aes(Dis, TPM))
p + geom_boxplot()

c('Normalized Expression (TPM)', 'Genenome type')
test <- t(row.of.interest[-c(1,2,3,12,13)])
Dis  <- c("CTR", "CTR", "CTR", "CTR", "ALS", "ALS", "ALS", "ALS")
test <- cbind(test,Dis)
test <- data.frame(test)
names(test) <- c('TPM', 'Dis')

test[1] <- as.numeric(levels(test[,1]))[test[,1]]
boxplot(TPM~Dis, data=test, main=row.of.interest$gene.names, xlab="Expression", ylab="Diseases state")

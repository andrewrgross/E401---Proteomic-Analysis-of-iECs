### PCA for E401 -- Andrew R Gross -- 11JUN21
### INPUT: Expression data; 
### OUTPUT: PCA Plots

############################################################################################
### Header
library(ggplot2)
library(ggbiplot)

############################################################################################
### Functions
reassign.protein.desc <- function(data.frame) {
  protein.list.reordered <- protein.list[names(data.frame),]
  data.frame <- cbind(data.frame, protein.list.reordered)
  return(data.frame)
}

############################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/')
list.files()
iec.data <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/2021_64_DataReport_Optra_PanHuman_mapDIA.csv', fileEncoding="UTF-8-BOM")
iec.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/E401-metadata.csv', fileEncoding="UTF-8-BOM")

############################################################################################
### Format
### Select data to plot
iec.metadata <- iec.metadata[c(1,2,3,11,12,13,14,15,16,17,18,19),]

### Rename rows
row.names(iec.data) <- iec.data$Protein
### Generate Protein name and description table
protein.list <- iec.data[c(2,3)]
### Remove extraneous columns
iec.data.f <- iec.data[-c(1,2,3,23,24,25,26,27)]
iec.data.f <- iec.data[iec.metadata$Sample]
summary(iec.data.f)
### Replace NAs with 0
iec.data.f[is.na(iec.data.f)] <- 1
iec.data.max <- apply(iec.data.f, 1, max)
rows.to.keep <- iec.data.max > 5000
summary(rows.to.keep)
### Reorder columns
results.iec <- iec.data.f[rows.to.keep,]

### Convert to matrix
results.iec <- as.matrix(results.iec)
summary(results.iec)
############################################################################################
### Calculate Principle Components

### Calculate the actual components
pca.iec <- prcomp(t(results.iec), scale = TRUE)
#pca.iec <- prcomp(t(results.iec), scale = FALSE)

### Calculate the percent variation accounted for by each component
pca.data.var <- pca.iec$sdev^2
pca.data.var.per <- round(pca.data.var/sum(pca.data.var)*100, 1)

############################################################################################
### Plot Data
title = 'PCA analysis of Endothelial Proteomes vs iPSCs'
subtitle = ''

#plot(pca.iec$x[,1], pca.iec$x[,2])
#plot(pca.iec$x[,2], pca.iec$x[,3])

barplot(pca.data.var.per, main = 'Scree Plot', xlab = 'Principle Component', ylab = 'Percent Variation')

### Define data frame for ggplot
pca.data.to.plot <- data.frame(Sample = rownames(pca.iec$x), 
                              PC1 = pca.iec$x[,1],
                              PC2 = pca.iec$x[,2],
                              PC3 = pca.iec$x[,3],
                              PC4 = pca.iec$x[,4],
                              PC5 = pca.iec$x[,5])

(pca.data.to.plot <- cbind(pca.data.to.plot, iec.metadata[-1]))
#pca.data.to.plot$Group <- pca.data.to.plot$CellLine

### Basic plot of PC1 v PC2
(pca.plot <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC2, label = Group, color = Group)) +
  geom_point() + 
  geom_text() + 
  xlab(paste('PC1 - ', pca.data.var.per[1], '%', sep = '')) +
  ylab(paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
  xlim(c(-50,70)) +
  theme_bw() +
  ggtitle('PCA of E401'))

############################################################################################
### Generate formatted pca plots
### PC1 v PC2
subtitle = 'PC1 v PC2'
(pca.1.v.2 <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC2, label = ShortName)) +
  geom_point(size = 5, aes(fill = Group), color = 'black', pch = 21) +
  #geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_fill_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
  labs(title = title, 
       subtitle = subtitle,
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC2 - ', pca.data.var.per[2], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) )


### PC1 v PC3
subtitle = 'PC1 v PC3'
(pca.1.v.3 <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC3, label = ShortName)) +
  geom_point(size = 5, aes(fill = Group), color = 'black', pch = 21) +
  #geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_fill_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
  labs(title = title, 
       subtitle = subtitle,
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC3 - ', pca.data.var.per[3], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) )

### PC2 v PC3
subtitle = 'PC2 v PC3'
(pca.2.v.3 <- ggplot(data = pca.data.to.plot, aes(x = PC2, y = PC3, label = ShortName)) +
    geom_point(size = 5, aes(fill = Group), color = 'black', pch = 21) +
    #geom_text(hjust=-0.2,vjust=0.5) + 
    xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
    scale_fill_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
    labs(title = title, 
         subtitle = subtitle,
         x = paste('PC2 - ', pca.data.var.per[2], '%', sep = ''), 
         y = paste('PC3 - ', pca.data.var.per[3], '%', sep = '')) +
    theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
          axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
          axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
          panel.background = element_rect(fill = 'white', color = 'black'),
          plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) )

### PC1 v PC4
subtitle = 'PC1 v PC4'
(pca.1.v.4 <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC4, label = ShortName)) +
  geom_point(size = 5, aes(fill = Group), color = 'black', pch = 21) +
  #geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_fill_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
  labs(title = title, 
       subtitle = subtitle,
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC4 - ', pca.data.var.per[4], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) )

### PC1 v PC5
subtitle = 'PC1 v PC5'
(pca.1.v.5 <- ggplot(data = pca.data.to.plot, aes(x = PC1, y = PC5, label = ShortName)) +
  geom_point(size = 5, aes(fill = Group), color = 'black', pch = 21) +
  #geom_text(hjust=-0.2,vjust=0.5) + 
  xlim(min(pca.data.to.plot$PC1)-5,max(pca.data.to.plot$PC1+20)) +
  scale_fill_manual(values = c("#fe1c1c", "#fab9b6", "#a6acf7", "#000dc4")) +
  labs(title = title, 
       subtitle = subtitle,
       x = paste('PC1 - ', pca.data.var.per[1], '%', sep = ''), 
       y = paste('PC5 - ', pca.data.var.per[5], '%', sep = '')) +
  theme(plot.title = element_text(color="black", face="bold", size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=14,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=14,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black'),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 12)) )



############################################################################################
### Write to folder
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/PCA/')


### Save plot
tiff(filename= paste0('PCA 1v2.tiff'), width = 800, height = 800, units = "px", pointsize = 12)
pca.1.v.2
dev.off()

tiff(filename= paste0('PCA 1v3.tiff'), width = 800, height = 800, units = "px", pointsize = 12)
pca.1.v.3
dev.off()

tiff(filename= paste0('PCA 1v4.tiff'), width = 800, height = 800, units = "px", pointsize = 12)
pca.1.v.4
dev.off()

tiff(filename= paste0('PCA 1v5.tiff'), width = 800, height = 800, units = "px", pointsize = 12)
pca.1.v.5
dev.off()


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


### Identify the genes with the largest influence
# PC1
l.score.pc1 <- pca.iec$rotation[,1]
l.pc1.up <- sort(l.score.pc1, decreasing = TRUE)
l.pc1.down <- sort(l.score.pc1, decreasing = FALSE)
l.pc1.abs <- l.score.pc1[order(abs(l.score.pc1), decreasing = TRUE)]
quantile(l.pc1.abs[,1], c(0.1, 0.2, 0.8, 0.9, 0.95))

### Reassign protein descriptors
l.pc1.up <- reassign.protein.desc(l.pc1.up)
l.pc1.down <- reassign.protein.desc(l.pc1.down)
l.pc1.abs <- reassign.protein.desc(l.pc1.abs)

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/PCA/PCA-loadings-tables/')
write.csv(l.pc1.up, 'PrComp1-loadings-UP.csv')
write.csv(l.pc1.down, 'PrComp1-loadings-DOWN.csv')
write.csv(l.pc1.abs, 'PrComp1-loadings-ABS_VALUES.csv')

# PC2
l.score.pc2 <- pca.iec$rotation[,2]
l.score.pc2.ranked <- sort(abs(l.score.pc2), decreasing = TRUE)
l.score.pc2[names(l.score.pc2.ranked)][1:10]

# PC3
l.score.pc3 <- pca.iec$rotation[,3]
l.score.pc3.ranked <- sort(abs(l.score.pc3), decreasing = TRUE)
l.score.pc3[names(l.score.pc3.ranked)][1:10]

# PC4
l.score.pc4 <- pca.iec$rotation[,4]
l.score.pc4.ranked <- sort(abs(l.score.pc4), decreasing = TRUE)
l.score.pc4[names(l.score.pc4.ranked)][1:10]

# PC5
l.score.pc5 <- pca.iec$rotation[,5]
l.pc5.up <- sort(l.score.pc5, decreasing = TRUE)
l.pc5.down <- sort(l.score.pc5, decreasing = FALSE)
l.pc5.abs <- l.score.pc5[order(abs(l.score.pc5), decreasing = TRUE)]

### Reassign protein descriptors
l.pc5.up <- reassign.protein.desc(l.pc5.up)
l.pc5.down <- reassign.protein.desc(l.pc5.down)
l.pc5.abs <- reassign.protein.desc(l.pc5.abs)

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/PCA/PCA-loadings-tables/')
write.csv(l.pc5.up, 'PrComp5-loadings-UP.csv')
write.csv(l.pc5.down, 'PrComp5-loadings-DOWN.csv')
write.csv(l.pc5.abs, 'PrComp5-loadings-ABS_VALUES.csv')



##############################################
############################################################################################
### sCRATCHWORK
############################################################################################
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

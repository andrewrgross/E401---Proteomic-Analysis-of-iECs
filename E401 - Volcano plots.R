### DESeq2 Pipeline - Andrew R Gross - 25MAY21
### A Volcano plot generation script  
### INPUT: This script requires a counts table.  
### OUTPUT: This script generates a table normalized expression with fold change and p-values between sample groups.

####################################################################################################################################################
### Header
#library("pasilla")
#library("DESeq2")
#library("biomaRt")
library("heatmap2")
#library("VennDiagram")
library (ggplot2)
library(gplots)

####################################################################################################################################################
### Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/Differential Expression Tables/DE csv/')

### Input DE results
iec.de <- read.csv('DE_iEC & HUVECs vs. iPSCs.csv', header = TRUE, row.names = 1) ; title = 'iEC & HUVECs vs. iPSCs'
iec.de <- read.csv('DE_iEC vs. HUVECs.csv', header = TRUE, row.names = 1) ; title = 'iEC vs. HUVECs'
iec.de <- read.csv('DE_iEC vs. iPSCs.csv', header = TRUE, row.names = 1) ; title = 'iEC vs. iPSCs'

iec.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/E401-metadata.csv', fileEncoding="UTF-8-BOM")

####################################################################################################################################################
### Format
##########################################################################
### Reassign names
#names(counts.cov) <- sample.names
geneNames <- iec.de$gene.table.row.names.iec.de....

geneNames = c()
dataframe <- iec.de
for (row.num in 1:(nrow(dataframe))){
  current.row <- dataframe[row.num,]
  firstGene <- strsplit(current.row[,1], split = ' ')[[1]][1]
  geneNames <- c(geneNames, firstGene)
}
geneNames <- make.unique(geneNames)
nrow(iec.de) == length(unique(geneNames))

row.names(iec.de) <- geneNames

### Filter p-value expression genes
summary(iec.de)
iec.de <- iec.de[iec.de$padj < 0.05,]
iec.de <- iec.de[order(iec.de$log2FoldChange),]

### Reorder columns and convert to matrix
expression.df <- iec.de[8:ncol(iec.de)]

newNames             <- c('PosCtrl-HUVEC 1',     'PosCtrl-HUVEC 2',     'PosCtrl-HUVEC 3',
                          'NegCtrl - EDI028 2',  'NegCtrl - EDI028 3',  
                          'NegCtrl - EDI42A 1',  'NegCtrl - EDI42A 3',  
                          'NegCon_03n14 1',      'NegCon_03n14 2', 'NegCon_03n14 3',
                          'iEC(d21) - EDI028 1', 'iEC(d21) - EDI028 2', 'iEC(d21) - EDI028 3',
                          'iEC(d21) - EDI42A 1', 'iEC(d21) - EDI42A 2', 'iEC(d21) - EDI42A 3',
                          'iEC(d21) - 03n14 1',  'iEC(d21) - 03n14 2',  'iEC(d21) - 03n14 3')
data.frame(names = names(expression.df), newNames)
names(expression.df) <- newNames


expression.df <- expression.df[c(1,2,3, 11,12,13,14,15,16,17,18,19, 4,5,6,7,8,9,10)]

names(expression.df)

expression.m <- as.matrix(expression.df)

####################################################################################################################################################
### Format Results


##########################################################################
### Heatmaps

pal <- colorRampPalette(c('yellow','black','blue'))(100)

#heatmap.2(expression.m, scale = 'row', col = pal, trace = 'none', dendrogram="none", 
#          Rowv=FALSE, symm=TRUE, density.info='none', labRow=NA,
#          lmat=rbind(c(4, 2), c(1, 3)), lhei=c(2, 8), lwid=c(4, 1), margins = c(9,0))

heatmap(expression.m, col = pal, Rowv = NA, Colv = NA, scale = 'row', margins = c(8,1))
##########################################################################
### Volcano Plots

### Format #############################################################################################
volcano.data.cov <- iec.de
p.less.than = volcano.data.cov$padj <0.5
p.less.than[is.na(p.less.than)] = FALSE
volcano.data.cov <- volcano.data.cov[p.less.than,]
volcano.data.cov$Gene = volcano.data.cov[,1]
volcano.data.cov$Gene = rownames(volcano.data.cov)
genes.ds <- volcano.data.cov$Gene[1:20]    # Genes of interest

### Converting FC from log2 to log10
volcano.data.cov$log10FC <- log10(2^volcano.data.cov$log2FoldChange)

volcano.data.cov.sig <- volcano.data.cov[volcano.data.cov$padj<0.005,]
volcano.data.cov.sig.up <- volcano.data.cov.sig[volcano.data.cov.sig$log2FoldChange>0.5,]
volcano.data.cov.sig.down <- volcano.data.cov.sig[-volcano.data.cov.sig$log2FoldChange>0.5,]

volcano.data.cov.sig <- volcano.data.cov.sig[order(volcano.data.cov.sig$baseMean),]


### Selecting genes to list
volcano.data.cov.sig$Edge.up <- -log(volcano.data.cov.sig$padj) * volcano.data.cov.sig$log2FoldChange
volcano.data.cov.sig$Edge.down <- -log(volcano.data.cov.sig$padj) * -volcano.data.cov.sig$log2FoldChange

volcano.text.up <- volcano.data.cov.sig[volcano.data.cov.sig$Edge.up > 9,]; nrow(volcano.text.up)
volcano.text.down <- volcano.data.cov.sig[volcano.data.cov.sig$Edge.down > 9,]; nrow(volcano.text.down)

### Genes of interest
volcano.genes.of.interest <- volcano.data.cov.sig[volcano.data.cov.sig$Gene %in% genes.ds,]
volcano.genes.of.interest.null <- volcano.data.cov[volcano.data.cov$Gene %in% genes.ds,]

### Manually filtering
#genes.to.omit <- c('ANKRD36', 'AC022400.7', 'LINC00342' )
#volcano.text.up <- volcano.text.up[!volcano.text.up$Gene %in% genes.to.omit,]
#volcano.text.down <- volcano.text.down[!volcano.text.down$Gene %in% genes.to.omit,]
text.size = 2


volcano <- ggplot( ) +
  geom_point(data = volcano.data.cov, aes(x=log2FoldChange, y = log10(padj) ),color = 'grey', size = 0.9) +
  geom_point(data = volcano.data.cov.sig.up, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'red') +
  geom_point(data = volcano.data.cov.sig.down, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'blue') +
#  geom_text(data = volcano.text.up, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), 
#            hjust = 0, vjust = 0.1, size = 3, check_overlap = TRUE) +
#  geom_text(data = volcano.text.down, aes(x=log2FoldChange - 0.2, y = log10(padj) - 0, label = Gene), 
#            hjust = 1, vjust = 0.1, size = 3, check_overlap = TRUE) +
  geom_point(data = volcano.genes.of.interest, aes(x=log2FoldChange, y = log10(padj), size = log10(baseMean)), color = 'yellow') +
  geom_text(data = volcano.genes.of.interest, aes(x=log2FoldChange + 0.2, y = log10(padj) - 0, label = Gene), hjust = 0, vjust = 0.1, size = 5, check_overlap = TRUE) +  
  scale_color_gradient(low="pink", high="red") +
  scale_size('Log10 Expression', range = c(0.5,4)) +
  ylim(c(0, min(log10(volcano.data.cov$padj))-2)) +
  xlim(c(min(volcano.data.cov$log2FoldChange-0.1), max(volcano.data.cov$log2FoldChange)+2)) +
  labs(title= title,
       x = 'Log2 Fold Change', 
       y = 'Log10 P-value (Adjusted)') +
  theme(plot.title = element_text(color="black", face="bold", hjust = 0.5, size=22, margin=margin(10,0,20,0)),
        axis.title.x = element_text(face="bold", size=18,margin =margin(20,0,10,0)),
        axis.title.y = element_text(face="bold", size=18,margin =margin(0,20,0,10)),
        panel.background = element_rect(fill = 'white', color = 'black', size = 3),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"), axis.text = element_text(size = 16),
        legend.position = 'none') 


volcano

############################################################################################
### Write to folder
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/Volcano Plots/')


### Save plot
tiff(filename= paste0('Volcano - ',title, '.tiff'), width = 2000, height = 1600, units = "px", pointsize = 20, res = 250)
volcano
dev.off()



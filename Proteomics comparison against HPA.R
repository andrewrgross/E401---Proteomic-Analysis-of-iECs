### Proteomics comparison to HPA -- Andrew R Gross -- 13MAY21
### A comparison of iec data against known protein lists from the Human Protein Atlas

#########################################################################################################
### 1. Header
library("DESeq2")
library(zoo)
library(pathfindR)
library("biomaRt")

ensembl = useMart(biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

#########################################################################################################
### 2. Functions
rm.na <- function(dataframe){
  sum.column <- apply(dataframe, 1, sum)
  na.rows <- is.na(sum.column)
  return(dataframe[!na.rows,])
}

convert.proteins.to.genes <- function(dataframe){
  dataframe.genes.expanded <- dataframe[0,]
  for (row.num in 1:(nrow(dataframe))){
    current.row <- dataframe[row.num,]
    gene.text <- current.row[,1]
    gene.text <- gsub(',', '', gene.text)
    gene.text <- gsub('/', ' ', gene.text)
    genelist <- strsplit(gene.text, split = ' ')[[1]]
    gene.num = length(genelist)
    for (num in 1:gene.num){
      current.row[1] <- genelist[num]
      dataframe.genes.expanded <- rbind(dataframe.genes.expanded, current.row)
    }
  }
  return(dataframe.genes.expanded)
}

assign.go.terms <- function(dataframe) {                                                       # Assign go terms to all rows in a dataframe with protein ID row names
  output.df <- data.frame(uniprotswissprot = '', bio.process = '', molecular.fun = '', cell.comp = '')                        # Define an empty data frame
  go.output <- getBM(attributes=c('uniprotswissprot','name_1006','namespace_1003'), filters='uniprotswissprot', values=row.names(dataframe), mart=ensembl)    # perform a biomart search
  
  unique.proteins <- unique(go.output$uniprotswissprot)                                        # For each protein, isolate the GO terms reported
  
  for(current.protein in unique.proteins) {
    current.df <- go.output[which(go.output$uniprotswissprot == current.protein),]
    bio.process <- current.df$name_1006[which(current.df$namespace_1003=='biological_process')]
    bio.process <- paste(bio.process, collapse = '; ')
    molecular.fun <- current.df$name_1006[which(current.df$namespace_1003 == 'molecular_function')]
    molecular.fun <- paste(molecular.fun, collapse = '; ')
    cell.comp <- current.df$name_1006[which(current.df$namespace_1003=='cellular_component')]
    cell.comp <- paste(cell.comp, collapse = '; ')
    
    new.row <- data.frame(uniprotswissprot = current.protein, bio.process, molecular.fun, cell.comp)
    output.df <- rbind(output.df, new.row)
  }
  return(output.df)
}

#########################################################################################################
### 3. Input
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iEC cell pellets and CM for Proteomics core/Proteomics results/')
list.files()
iec.data <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/2021_64_DataReport_Optra_PanHuman_mapDIA.csv')
iec.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/E401-metadata.csv')

hpa.reference.data <- read.table('C:/Users/grossar/Box/Sareen Lab Shared/Data/Reference Data/proteinatlas.tsv', sep = '\t', header = TRUE)
#hpa.reference.data <- read.table('C:/Users/grossar/Desktop/proteinatlas.tsv', sep = '\t', header = TRUE)

#########################################################################################################
### 4. Compose lists of rows with proteins detected in all samples of a type  
detected.in.all.iecs <- which(iec.data$COUNT_iEC >= 8)          ;       length(detected.in.all.iecs)
detected.in.all.huvecs <- which(iec.data$COUNT_PC >= 3)         ;       length(detected.in.all.huvecs)
detected.in.all.ipscs <- which(iec.data$COUNT_NC >= 6)          ;       length(detected.in.all.ipscs)
detected.in.all.ec <- intersect(detected.in.all.iecs, detected.in.all.huvecs)    ;    length(detected.in.all.ec)
detected.in.all.total <- intersect(detected.in.all.ec, detected.in.all.ipscs)    ;    length(detected.in.all.total)

detected.only.in.iecs <- setdiff(detected.in.all.iecs, union(detected.in.all.huvecs, detected.in.all.ipscs))  ;  length(detected.only.in.iecs)

detected.only.in.huvecs <- setdiff(detected.in.all.huvecs, union(detected.in.all.iecs, detected.in.all.ipscs))  ;  length(detected.only.in.huvecs)
detected.only.in.ec <- setdiff(detected.in.all.ec,detected.in.all.ipscs)          ;  length(detected.only.in.ec)

#########################################################################################################
### 5. Formatting
## - Remove unwanted columns
gene.table <- iec.data[2]
rownames(gene.table) <- iec.data$ï..Protein

iec.data <- iec.data[-c(1,2,3,23,23,24,25,26,27)]
rownames(iec.data) <- rownames(gene.table)

hpa.cell.type <- hpa.reference.data[-seq(17,241)]
hpa.ec <- hpa.cell.type[c(1,2,3,4,5,8,9,10,11,33)]

#########################################################################################################
### 6: Analysis
### - Filter based on expression level

sig.expression <- hpa.ec$Single.Cell.Type.RNA...Endothelial.cells..NX.> quantile(hpa.ec$Single.Cell.Type.RNA...Endothelial.cells..NX, na.rm = TRUE)[4][[1]]

hpa.ec <- hpa.ec[sig.expression,]
hpa.ec <- hpa.ec[order(hpa.ec$Single.Cell.Type.RNA...Endothelial.cells..NX.,decreasing = TRUE),]

iec.prot <- iec.data[grep('iEC',iec.metadata$Group)]
iec.prot$mean <- apply(iec.prot[seq(1,9)],1,mean)
iec.prot$median <- apply(iec.prot[seq(1,9)],1,median)
iec.prot$sd <- apply(iec.prot[seq(1,9)],1,sd)
iec.prot <- iec.prot[order(iec.prot$mean, decreasing = TRUE),]

### - Identify intersection of lists
shared <- intersect(rownames(iec.prot), hpa.ec$Uniprot)
shared.ec.prot <- iec.prot[shared,]
shared.annotations <- hpa.ec[c(4,6,7,8)][match(shared,hpa.ec$Uniprot),]
rownames(shared.annotations) <- shared
shared.ec.prot <- cbind(shared.ec.prot,shared.annotations)

######################################################################################
### 7: Differential Expression

### 7.1: Select data
# Counts of iEC vs HUVECs
counts.data <- iec.data[c(1,2,3,11,12,13,14,15,16,17,18,19)] ; title = 'iEC vs. HUVECs'; columnData <- data.frame(rownames = names(counts.data), celltype = c("HUVEC","HUVEC","HUVEC","iEC","iEC","iEC","iEC","iEC","iEC","iEC","iEC","iEC"))

# Counts of iEC & HUVECs vs iPSC
counts.data <- iec.data ; title = 'iEC & HUVECs vs. iPSCs'; columnData <- data.frame(rownames= names(counts.data), celltype = c('EC','EC','EC','iPSC','iPSC','iPSC','iPSC','iPSC','iPSC','iPSC','EC','EC','EC','EC','EC','EC','EC','EC','EC'))

# Counts of iEC vs iPSC
counts.data <- iec.data[4:19] ; title = 'iEC vs. iPSCs'; columnData <- data.frame(rownames= names(counts.data), celltype = c('iPSC','iPSC','iPSC','iPSC','iPSC','iPSC','iPSC','EC','EC','EC','EC','EC','EC','EC','EC','EC'))

# Remove NAs
counts.data <- rm.na(counts.data)

# Convert counts into a matrix of integers
counts.data.m <- round(as.matrix(counts.data), 0)

####################################################################################################################################################
### 8: Differential Expression
### Make our DESeq data sets
dds.iec <- DESeqDataSetFromMatrix(countData = counts.data.m, colData = columnData, design = ~ celltype)
### Run DESeq
dds.iec <- DESeq(dds.iec)
iec.de <- as.data.frame(results(dds.iec))
### Join intensity to data
genes <- as.data.frame(gene.table[row.names(iec.de),])
row.names(genes) <- row.names(iec.de)
iec.de <- cbind(genes,iec.de, counts.data)

####################################################################################################################################################
### 9: Format Results
### Order by padj
iec.de <- iec.de[order(iec.de$padj, decreasing = FALSE),]

### Output distribution
hist(log10(iec.de$padj[iec.de$padj<0.5]), breaks = seq(round(log10(min(iec.de$padj[1]))), 0))

length(which(iec.de$padj<0.5))
length(which(iec.de$padj<0.005))

### Assign GO terms
go.results <- assign.go.terms(iec.de)
row.names(go.results) <- go.results$uniprotswissprot
### Join
missing.rows <- data.frame(uniprotswissprot = setdiff(row.names(iec.de), row.names(go.results)), bio.process = '', molecular.fun = '', cell.comp = '')
row.names(missing.rows) <- missing.rows$uniprotswissprot

go.results <- rbind(go.results, missing.rows)
go.results <- go.results[row.names(iec.de),]

iec.de <- cbind(iec.de, go.results)

####################################################################################################################################################
### Output DE results
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iEC cell pellets and CM for Proteomics core/Proteomics results/E401 - Analysis of proteomics data/')
write.csv(iec.de, paste0('DE_',title,'-with GO.csv'))

### Re-Input results
iec.de <- read.csv('DE_iEC & HUVECs vs. iPSCs.csv', header = TRUE, row.names = 1) ; title = 'iEC & HUVECs vs. iPSCs'
iec.de <- read.csv('DE_iEC vs. HUVECs.csv', header = TRUE, row.names = 1) ; title = 'iEC vs. HUVECs'
iec.de <- read.csv('DE_iEC vs. iPSCs.csv', header = TRUE, row.names = 1) ; title = 'iEC vs. iPSCs'


####################################################################################################################################################
### Expand Gene lists
### Separate upregulated and downregulated top 100
print(title)
nrow(iec.de)
iec.de.up <- iec.de[which(iec.de$log2FoldChange>0),]
iec.de.down <- iec.de[which(iec.de$log2FoldChange<0),]

iec.de.up <- iec.de.up[1:100,]
iec.de.down <- iec.de.down[1:100,]

iec.de.up.expanded <- convert.proteins.to.genes(iec.de.up)
nrow(iec.de.up.expanded)
write.csv(iec.de.expanded$gene.table.row.names.iec.de...., paste0('genes-for-pathway-analysis_',title,'-UP.csv'), row.names = FALSE)


iec.de.down.expanded <- convert.proteins.to.genes(iec.de.down)
nrow(iec.de.down.expanded)
write.csv(iec.de.down.expanded$gene.table.row.names.iec.de...., paste0('genes-for-pathway-analysis_',title,'-DOWN.csv'), row.names = FALSE)



hpa.filtered <- hpa.ec[1:100,]
hpa.filtered$Gene <- paste(hpa.filtered$Gene, hpa.filtered$Gene.synonym)
dataframe <- hpa.filtered[1:10,]



hpa.filt.expanded <- convert.proteins.to.genes(hpa.filtered)
nrow(hpa.filt.expanded)
write.csv(hpa.filt.expanded, paste0('genes-for-pathway-analysis_HPA100.csv'), row.names = FALSE)

top <- head(iec.de, 30)
test <- iec.de
input_df <- head(iec.de,20)


####################################################################################################################################################
### Pathway identifier using pathfindR
### testing
output_df <- run_pathfindR(input_df)



### Output
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/E401 - Analysis of proteomics data')
write.csv(shared.ec.prot, 'Proteins shared by iecs and HPA ECs.csv')
write.csv(top, 'Proteins upregulated in iEC vs HUVECs.csv')

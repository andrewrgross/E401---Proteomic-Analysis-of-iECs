### Proteomics comparison to HPA -- Andrew R Gross -- 13MAY21
### A comparison of iec data against known protein lists from the Human Protein Atlas

#########################################################################################################
### 1. Header
#########################################################################################################
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
#########################################################################################################
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
#########################################################################################################
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/')
list.files()
#iec.data <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/Results/2021/i-ECs_2021/iEC cell pellets and CM for Proteomics core/Proteomics results/2021_64_DataReport_Optra_PanHuman_mapDIA.csv', fileEncoding="UTF-8-BOM")
iec.data <-read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Proteomics Datasets/DIANN_iPSC_EC_TotalProteome_ProteotypicIdentifications.csv', fileEncoding = 'UTF-8-BOM')

iec.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/E401-metadata.csv', fileEncoding = 'UTF-8-BOM')

hpa.reference.data <- read.table('C:/Users/grossar/Box/Sareen Lab Shared/Data/Reference Data/proteinatlas.tsv', sep = '\t', header = TRUE)

#########################################################################################################
### 4. Compose lists of rows with proteins detected in all samples of a type  
#########################################################################################################
### 4.1: Format expression data
gene.table <- iec.data[2]
rownames(gene.table) <- iec.data$Protein
rownames(iec.data) <- rownames(gene.table)


### 4.2 - Formate HPA data
hpa.cell.type <- hpa.reference.data[-seq(17,241)]
hpa.ec <- hpa.cell.type[c(1,2,3,4,5,8,9,10,11,33)]

### Isolating relevant proteins:
## Find the quartiles:
all_expression = iec.data#[4:22]

all_values = c()
for (columnno in seq(1,ncol(all_expression))){
  column = all_expression[,columnno]
  all_values = c(all_values, unlist(column))
}
length(all_expression)
summary(all_expression)

# Remove all NAs and all proteins with median expression below the first quartile
test <- which(apply(iec.data[11:19],1,min) >= 7593)



### I need 10 sets. I'll make a data frame to hold them.
vennDdata = data.frame(counts = rep(0,13))
rownames(vennDdata) = c('iec', 'ipsc', 'huvec', 'in_iec_hu', 'in_iec_ipsc', 'in_ipsc_hu', 'in_all', 'only_iec_hu', 'only_iec_ipsc', 'only_ipsc_hu', 'only_iec', 'only_hu', 'only_ipsc')

detected_in_iec <- which(iec.data$COUNT_iEC >= 5)          ;       vennDdata['iec',] = length(detected_in_iec)
detected_in_huvec <- which(iec.data$COUNT_PC >= 2)         ;       vennDdata['huvec',] = length(detected_in_huvec)
detected_in_ipsc <- which(iec.data$COUNT_NC >= 4)          ;       vennDdata['ipsc',] = length(detected_in_ipsc)

intersection_iec_huvec <- intersect(detected_in_iec, detected_in_huvec)    ;  vennDdata['in_iec_hu',]   = length(intersection_iec_huvec)
intersection_iec_ipsc  <- intersect(detected_in_iec, detected_in_ipsc)     ;  vennDdata['in_iec_ipsc',] = length(intersection_iec_ipsc)
intersection_ipsc_hu   <- intersect(detected_in_ipsc, detected_in_huvec)   ;  vennDdata['in_ipsc_hu',]  = length(intersection_ipsc_hu)

intersection_all       <- intersect(intersection_iec_huvec, intersection_ipsc_hu) ;  vennDdata['in_all',] = length(intersection_all)

only_iec_hu   <- setdiff(intersection_iec_huvec, intersection_all)  ;  vennDdata['only_iec_hu',]   = length(only_iec_hu)
only_iec_ipsc <- setdiff(intersection_iec_ipsc,  intersection_all)  ;  vennDdata['only_iec_ipsc',] = length(only_iec_ipsc)
only_ipsc_hu  <- setdiff(intersection_ipsc_hu,   intersection_all)  ;  vennDdata['only_ipsc_hu',]  = length(only_ipsc_hu)

only_iec      <- setdiff(detected_in_iec,union(detected_in_ipsc, detected_in_huvec)) ; vennDdata['only_iec',] = length(only_iec)
only_ipsc     <- setdiff(detected_in_ipsc,union(detected_in_iec, detected_in_huvec)) ; vennDdata['only_ipsc',] = length(only_ipsc)
only_hu       <- setdiff(detected_in_huvec,union(detected_in_ipsc, detected_in_iec)) ; vennDdata['only_hu',] = length(only_hu)

vennDdata

#########################################################################################################
### 5 - Generate tables
#########################################################################################################

data.iec <- iec.data[only_iec,]
data.hu  <- iec.data[only_hu,]
data.ipsc <- iec.data[only_ipsc,]
data.iec_hu <- iec.data[only_iec_hu,]
data.iec_ipsc <- iec.data[only_iec_ipsc,]
data.ipsc_hu <- iec.data[only_ipsc_hu,]
data.all  <- iec.data[intersection_all,]


#########################################################################################################
### 6 - Export tables
#########################################################################################################
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/Venn diagram tables/')

vennDdata
write.csv(data.iec, paste0('Venn_iec.csv'))
write.csv(data.hu, paste0('Venn_hu.csv'))
write.csv(data.ipsc, paste0('Venn_ipsc.csv'))
write.csv(data.all, paste0('Venn_detected_in_all.csv'))
write.csv(data.iec_hu, paste0('Venn_iec_hu_only.csv'))
write.csv(data.iec_ipsc, paste0('Venn_iec_ipsc_only.csv'))
write.csv(data.ipsc_hu, paste0('Venn_ipsc_hu_only.csv'))


#########################################################################################################
### 7: Differential Expression
#########################################################################################################
iec.data = iec.data[-c(1:3,24:29)]
names(iec.data) = iec.metadata$ShortName

### 7.1: Select data
# Counts of iEC vs HUVECs
counts.data <- iec.data[iec.metadata$CellType == 'EC'] ; title = 'iEC vs. HUVECs'; columnData <- data.frame(rownames = names(counts.data), celltype = c("HUVEC","HUVEC","HUVEC","iEC","iEC","iEC","iEC","iEC","iEC","iEC","iEC","iEC"))

# Counts of iEC & HUVECs vs iPSC
counts.data <- iec.data ; title = 'iEC & HUVECs vs. iPSCs'; columnData <- data.frame(rownames= names(counts.data), celltype =  iec.metadata$CellType)

# Counts of iEC vs iPSC
counts.data <- iec.data[iec.metadata$Group != 'HUVEC'] ; title = 'iEC vs. iPSCs'; columnData <- data.frame(rownames= names(counts.data), celltype = iec.metadata$Group[iec.metadata$Group != 'HUVEC'])

# Remove NAs
counts.data <- rm.na(counts.data)

# Convert counts into a matrix of integers
counts.data.m <- round(as.matrix(counts.data/100), 0)

#########################################################################################################
### 8: Differential Expression
#########################################################################################################
### Make our DESeq data sets
dds.iec <- DESeqDataSetFromMatrix(countData = counts.data.m, colData = columnData, design = ~ celltype)
### Run DESeq
dds.iec <- DESeq(dds.iec)
iec.de <- as.data.frame(results(dds.iec))
### Join intensity to data
genes <- as.data.frame(gene.table[row.names(iec.de),])
row.names(genes) <- row.names(iec.de)
iec.de <- cbind(genes,iec.de, counts.data)

#########################################################################################################
### 9: Format Results
#########################################################################################################
### Order by padj
iec.de <- iec.de[order(iec.de$padj, decreasing = FALSE),]
### Flip the sign of the FC if necessary to make sure the results are normalized to iPSCs
iec.de$log2FoldChange = -iec.de$log2FoldChange

### Output distribution
hist(log10(iec.de$padj[iec.de$padj<0.5]), breaks = seq(round(log10(min(iec.de$padj[1])))-1, 0))

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

#########################################################################################################
### 10: Output DE results
#########################################################################################################
print(title)
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E401 - Analysis of proteomics data/Differential Expression Tables/')
write.csv(iec.de, paste0('DE_',title,'-with GO.csv'))

### Re-Input results
iec.de <- read.csv('DE_iEC & HUVECs vs. iPSCs.csv', header = TRUE, row.names = 1) ; title = 'iEC & HUVECs vs. iPSCs'
iec.de <- read.csv('DE_iEC vs. HUVECs.csv', header = TRUE, row.names = 1) ; title = 'iEC vs. HUVECs'
iec.de <- read.csv('DE_iEC vs. iPSCs.csv', header = TRUE, row.names = 1) ; title = 'iEC vs. iPSCs'




#########################################################################################################
### X: Depreciated steps
#########################################################################################################


all_expression[]
detected_in_all.iecs <- which(iec.data$COUNT_iEC >= 8)          ;       length(detected_in_all.iecs)


#########################################################################################################
### 6: Analysis
### - Filter based on expression level

sig.expression <- hpa.ec$Single.Cell.Type.RNA...Endothelial.cells..NX.> quantile(hpa.ec$Single.Cell.Type.RNA...Endothelial.cells..NX, na.rm = TRUE)[4][[1]]

hpa.ec <- hpa.ec[sig.expression,]
hpa.ec <- hpa.ec[order(hpa.ec$Single.Cell.Type.RNA...Endothelial.cells..NX.,decreasing = TRUE),]

cutoff = 30000

### Compose list of proteins in iEC
iec.prot <- iec.data[grep('iEC',iec.metadata$Group)+3][detected_in_all.iecs,]
iec.prot$mean <- apply(iec.prot[seq(1,9)],1,mean)
iec.prot$median <- apply(iec.prot[seq(1,9)],1,median)
iec.prot$sd <- apply(iec.prot[seq(1,9)],1,sd)
iec.prot$min <- apply(iec.prot[seq(1,9)],1,min)
iec.prot <- iec.prot[iec.prot$min > cutoff,]
iec.prot <- iec.prot[order(iec.prot$median, decreasing = TRUE),] ; nrow(iec.prot)   # 3207

### Compose list of proteins in HUVEC
huvec.prot <- iec.data[grep('HUVEC', iec.metadata$Group) + 3][detected_in_all.huvecs,]
huvec.prot$median <- apply(huvec.prot[seq(1,3)], 1, median)
huvec.prot$min <- apply(huvec.prot[seq(1,3)], 1, min)
huvec.prot <- huvec.prot[huvec.prot$min > cutoff,]
huvec.prot <- huvec.prot[order(huvec.prot$median, decreasing = TRUE),] ; nrow(huvec.prot) # 3429

### Compose list of proteins in iPSC
ipsc.prot <- iec.data[grep('iPSC', iec.metadata$Group) + 3][detected_in_all.huvecs,]
ipsc.prot$median <- apply(ipsc.prot[seq(1,3)], 1, median)
ipsc.prot$min <- apply(ipsc.prot[seq(1,3)], 1, min)
ipsc.prot <- ipsc.prot[ipsc.prot$min > cutoff,]
ipsc.prot <- ipsc.prot[order(ipsc.prot$median, decreasing = TRUE),] ; nrow(ipsc.prot)   # 3666

### - Identify intersection of lists
shared.iec.huvec = intersect(rownames(iec.prot), rownames(huvec.prot))   # 2583
shared.iec.ipsc = intersect(rownames(iec.prot), rownames(ipsc.prot))   # 2583
shared.ipsc.huvec = intersect(rownames(ipsc.prot), rownames(huvec.prot)) #2333
shared.all = intersect(shared.iec.huvec, shared.ipsc.huvec)

### - Separate out proteins only in each group pair
only.iec.huvec = setdiff(shared.iec.huvec, shared.all)
only.iec.ipsc = setdiff(shared.iec.ipsc, shared.all)
only.ipsc.huvec = setdiff(shared.ipsc.huvec, shared.all)

### - Separate out proteins only in each single group
only.iec = setdiff(rownames(iec.prot), union(shared.iec.huvec, shared.iec.ipsc))        ; length(only.iec)
only.huvec = setdiff(rownames(huvec.prot), union(shared.iec.huvec, shared.ipsc.huvec))  ; length(only.huvec)
only.ipsc = setdiff(rownames(ipsc.prot), union(shared.iec.ipsc, shared.ipsc.huvec))     ; length(only.ipsc)
  

### Find intersection with HPA
shared <- intersect(rownames(iec.prot), hpa.ec$Uniprot)
shared.ec.prot <- iec.prot[shared,]
shared.annotations <- hpa.ec[c(4,6,7,8)][match(shared,hpa.ec$Uniprot),]
rownames(shared.annotations) <- shared
shared.ec.prot <- cbind(shared.ec.prot,shared.annotations)

#########################################################################################################
### 6.1 - Plot Eulerr/Venn diagram
### - Filter based on expression level



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

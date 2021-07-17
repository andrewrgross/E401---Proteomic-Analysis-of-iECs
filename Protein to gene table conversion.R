### Protein to gene table conversion -- Andrew R Gross -- 13MAY21
### A pipeline to convert a table of protein expression values into tables of individuals genes which makeup the proteins

### Header


### Functions

### Input

iec.data <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/2021_64_DataReport_Optra_PanHuman_mapDIA.csv')
iec.metadata <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/iec cell pellets and CM for Proteomics core/Proteomics results/E401-metadata.csv')

iec.de.data <- read.csv('C:/Users/grossar/Box/Sareen Lab Shared/Data/Roberta/2021/Results and reports/i-ECs/...')
protein.to.gene.lookup.table <- read.csv('C:/User/gros')

### Formatting
## - Remove unwanted columns
gene.table <- iec.data[2]
rownames(gene.table) <- iec.data$ï..Protein

### Expand all protein rows
## For each protein row, separate out the associated genes
genes.temp <- strsplit(geneNames, split = ' ')[[1]]

genes <- as.data.frame(gene.table[row.names(iec.de),])
row.names(genes) <- row.names(iec.de)

dataframe <- head(iec.de)

iec.de.genes.expanded <- iec.de[0,]

for (row.num in 1:(nrow(dataframe))){
  current.row <- dataframe[row.num,]
  genelist <- strsplit(current.row[,1], split = ' ')[[1]]
  gene.num = length(genelist)
  for (num in 1:gene.num){
    current.row[1] <- genelist[num]
    iec.de.genes.expanded <- rbind(iec.de.genes.expanded, current.row)
  }
}

input.test <- iec.de.genes.expanded[c(1,3,7)]
input.test <- head(input.test,100)

output_df <- run_pathfindR(input.test)

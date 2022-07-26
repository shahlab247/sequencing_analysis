
### Base analytical pipeline for analysis of RNA-seq read alignment
### Performed for evaluation of effects of immune stimuli on PAF1 knockout cells
### Author: Matthew W. Kenaston
### Date: January 2022

### ***plotting functions are NOT included in this summary of code.

### Folllowing packages and their dependants are used:
### biomaRt, DESeq2, fgsea, readr, dplyr


#### Differential gene expression analysis using DESEQ2 -----------------

#importing STAR alignment data from featureCounts matrix
source <- "[FILENAME]"
#create table of raw values
raw<-read.csv(source, header=TRUE, row.names = 'id') 
#remove low read counts
raw <- raw [which(rowMeans(raw) > 1), ]

#create ENSG to common gene name list for labeling use
id <- read.csv(source, header=TRUE)
id <- id[,c(1,2)]
head(id)

### clean up data frame of raw data
df <- raw
df <- as.matrix(df)
df <- na.omit(df) 
df <- df[!duplicated(df), ] 
head(df)

### Create factors to describe data
# NOTE: these factors are specific to the experiment this pipeline is designed around.
stim <- factor(c(rep("Mock", 3),rep("polyIC", 3), rep("polydAdT", 3), rep("IFNB", 3), 
                 rep("LPS", 3), rep("TNFA", 3), rep("Mock", 3),rep("polyIC", 3), rep("polydAdT", 3), rep("IFNB", 3), 
                 rep("LPS", 3), rep("TNFA", 3)))
cell <- factor(c(rep("A549", 18), rep("PAF1KO",18)))
day <- factor(rep(c(1,2,3),12))

### Establish a correlation between experiment and mock trials for DESeq2 to use for analysis
### Create DESeq Data Set from matrix and condition objects
coldata <- data.frame(row.names=colnames(df), cell, stim, day)
dds <- DESeqDataSetFromMatrix(countData=df, colData=coldata)
colnames(dds) <- paste(rownames(colData(dds)), as.character( dds$stim), 
                       as.character( dds$cell),sep="_"  )
colnames(dds) = gsub("sample","S",colnames(dds))
### specify design of the experiment
dds$group <- factor(paste0(dds$cell, dds$stim))
design(dds) <- ~ group + day
### make ncgRNA mock the base level
dds$group <- relevel(dds$group, ref = "A549Mock")

### DESEQ2 pipeline implemented
dds <- DESeq(dds, full=design(dds),minReplicatesForReplace = 3)

### Plot dispersions to qualitatively check data is appropriate for further analysis:
plotDispEsts(dds, main="Dispersion plot")
### Create DESEQ2 results object
res <- results(dds)

### PCA Plot -------------------------------
vsd <- vst(dds, blind=FALSE)
### remove batch effects as needed
new <- ComBat(assay(vsd), vsd$day)
assay(vsd) <- new
counts_batch_corrected <- assay(vsd)

### plot PCA
PCA <- plotPCA(vsd, intgroup=c('cell'), ntop=1000, returnData = T)

percentVar <- round(100 * attr(PCA, "percentVar")) 
ggplot(PCA, aes(x = PC1, y = PC2, color = stim, shape =cell)) + 
  geom_point(size =4, aes(fill=stim), alpha=0.7) + 
  geom_point(size =4, alpha = 0.7) + 
  xlab(paste("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste("PC2: ", percentVar[2], "% variance")) + 
  theme(text=element_text(size=14,  family="arial")) + theme_few() 

### Contrast  function ----------------------
### outputs are available in Table S1 of corresponding publication.
### applied for every contrast/comparison made

### form results based on a specific comparison
### NOTE: numbers adjusted by contrast, showing no contrast currently
res.XXX <- results(dds, contrast=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)) 
### create corresponding data frame
res.XXX_data <- merge(as.data.frame(res.XXX), as.data.frame(counts(dds, normalized=TRUE)), 
                      by="row.names", sort=FALSE)
### add common gene labels
sig_XXX <- data.frame(id[match(res.XXX_data$Row.names, 
                               id$id),], res.XXX_data[,2:42])


### GSEA function -------------------------
### use outputs created from previous contrast function
gseaDat <- read.csv('[FILENAME]')
### rank genes by "stat" variable, weighted statistic for log2FC and p.adjusted
ranks <- gseaDat$stat
names(ranks) <- gseaDat$gene
path <- ('GENESET')
### run fgsea on ranks against chosen geneset
fgseaRes <- fgsea(path, sort(ranks, decreasing = T), 
                  minSize=10, maxSize = 1000, nPermSimple = 5000)

### BioMart selection of immune genes for filtering of downstream analyses --------------
### Ensembl database
bm <- useMart("ensembl") 
### human genes and related terms
bm <- useDataset("hsapiens_gene_ensembl", mart=bm) 
### extract genes and Gene Ontology (GO) attributes
go <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id', 'entrezgene_id'),filter ="go_parent_term",
            values = c("GO:0002376"), useCache = T)




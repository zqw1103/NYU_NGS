#Step 1: Construct a data frame with sample information

sampleFiles <- paste(c('PDAC253','PDAC282','PDAC286','PDAC316','PDAC266','PDAC273','PDAC306','PDAC318'),'.htseq_count.txt',sep="")
sampleCondition <- c(rep('highSucrose',4),rep('lowSucrose',4))
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
View(sampleTable)

#Step 2: Read htseq-count files
library(DESeq2)
sampleTable$condition <- factor(sampleTable$condition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9",
                                  design = ~ condition)
dds
assay(dds)
#Q1.1. Report the output indicating what type of object the dds variable is, how many genes, and how many samples are stored in the object? [ 1 point ]
class(dds) # Determine the type of object
dds

#Task 2: Pre-filter low count genes, normalize counts and run DESeq
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
class(dds) # Determine the type of object
rawcounts.matrix <- counts(dds,normalized=F)
normalizedcounts.matrix <- counts(dds,normalized=T)
class(rawcounts.matrix)
View(normalizedcounts.matrix)

#Q2.1 Enter dds at the console to summarize your object. Report the output for your answer and how many genes were retained after removing the low count genes [ 1 point ].
dds
assays(dds)

#Q3.1 Include the hierarchical clustering result, your plotPCA command, and the PCA plot for you answer [ 1 point ].

rld <- rlog(dds)
assay(rld)
head(t(assay(rld)))

dists <- dist(t(assay(rld)))
hclust(dists)

png("week10_task2_hcluster.png")
plot(hclust(dists))
dev.off()


png("week10_task2_PCA.png")
plotPCA(rld, intgroup=c("condition"))
dev.off()


res <- results(dds, contrast = c('condition','lowSucrose','highSucrose') )
res05 <- results(dds, contrast = c('condition','lowSucrose','highSucrose'), alpha=0.05 )
class(res) 
res
summary(res) #0.1 as threshold
summary(res05) #0.05

metadata(res)$filterThreshold

resOrdered <- res[order(res$pvalue),] 
res05Ordered <- res05[order(res05$pvalue),] 
head(resOrdered,10) # View the 10 most significant genes
head(res05Ordered,10)
sum(is.na(resOrdered$pvalue))
resOrdered[is.na(resOrdered$padj),]
View(resOrdered[is.na(resOrdered$pvalue),])
resOrdered[is.na(resOrdered$log2FoldChange),]

install.packages("gridExtra")
library(gridExtra)
resOrdered_df_genes_interest <- resOrdered[ row.names(resOrdered) %in% c('Pdac_HC_chr14G0022900','Pdac_HC_chr14G0023100','Pdac_HC_chr14G0028200'), ]
res05Ordered_df_genes_interest <- res05Ordered[ row.names(res05Ordered) %in% c('Pdac_HC_chr14G0022900','Pdac_HC_chr14G0023100','Pdac_HC_chr14G0028200'), ]
getwd()
write.table(resOrdered_df_genes_interest, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/3 gene results table.csv")
write.table(res05Ordered_df_genes_interest, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/3 gene results table 0.05.csv")
read.table("3 gene results table.csv", header = T, sep = " ")
read.table("3 gene results table 0.05.csv", header = T, sep = " ")

Pdac_HC_chr14G0028200_counts_table <- plotCounts(dds, gene='Pdac_HC_chr14G0028200', intgroup="condition", 
           returnData=TRUE)
write.table(Pdac_HC_chr14G0028200_counts_table, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/Pdac_HC_chr14G0028200_counts_table.csv")
read.table("Pdac_HC_chr14G0028200_counts_table.csv", header = T, sep = " ")

png("Pdac_HC_chr14G0028200_counts_table.png")
plotCounts(dds, gene='Pdac_HC_chr14G0028200', intgroup="condition")
dev.off()





Pdac_HC_chr14G0022900_counts_table <- plotCounts(dds, gene='Pdac_HC_chr14G0022900', intgroup="condition", 
           returnData=TRUE)
write.table(Pdac_HC_chr14G0022900_counts_table, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/Pdac_HC_chr14G0022900_counts_table.csv")
read.table("Pdac_HC_chr14G0022900_counts_table.csv", header = T, sep = " ")

png("Pdac_HC_chr14G0022900_counts_table.png")
plotCounts(dds, gene='Pdac_HC_chr14G0022900', intgroup="condition")
dev.off()



Pdac_HC_chr14G0023100_counts_table <- plotCounts(dds, gene='Pdac_HC_chr14G0023100', intgroup="condition", 
                                                 returnData=TRUE)
write.table(Pdac_HC_chr14G0023100_counts_table, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/Pdac_HC_chr14G0023100_counts_table.csv")
read.table("Pdac_HC_chr14G0023100_counts_table.csv", header = T, sep = " ")

png("Pdac_HC_chr14G0023100_counts_table.png")
plotCounts(dds, gene='Pdac_HC_chr14G0023100', intgroup="condition")
dev.off()

resultsNames(dds)
res.shrunk <- lfcShrink(dds, coef = "condition_lowSucrose_vs_highSucrose" )
res.shrunkOrdered <- res.shrunk[order(res.shrunk$pvalue),]
res.shrunk.ashr <- lfcShrink(dds, contrast = c('condition','lowSucrose','highSucrose'), type = "ashr" )
res.shrunk.normal <- lfcShrink(dds, contrast = c('condition','lowSucrose','highSucrose'), type = "normal" )
res.shrunk.ashrOrdered <- res.shrunk.ashr[order(res.shrunk.ashr$pvalue),]
res.shrunk.normalOrdered <- res.shrunk.normal[order(res.shrunk.normal$pvalue),]
resOrdered
res.shrunkOrdered
res.shrunk.ashrOrdered
res.shrunk.normalOrdered

png("MA_unshrunken.png")
plotMA(res, ylim=c(-4,4))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()

png("MA_shrunken.png")
plotMA(res.shrunk,ylim=c(-4,4))
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()


res.shrunkOrdered_3_genes.table <- res.shrunkOrdered[ row.names(res.shrunkOrdered) %in% c('Pdac_HC_chr14G0022900','Pdac_HC_chr14G0023100','Pdac_HC_chr14G0028200'), ]
res.shrunk.ashrOrdered_3_genes.table <- res.shrunk.ashrOrdered[ row.names(res.shrunk.ashrOrdered) %in% c('Pdac_HC_chr14G0022900','Pdac_HC_chr14G0023100','Pdac_HC_chr14G0028200'), ]
res.shrunk.normalOrdered_3_genes.table <- res.shrunk.normalOrdered[row.names(res.shrunk.normalOrdered) %in% c('Pdac_HC_chr14G0022900','Pdac_HC_chr14G0023100','Pdac_HC_chr14G0028200'), ]
write.table(res.shrunkOrdered_3_genes.table, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/res.shrunkOrdered_3_genes.table.csv")
write.table(res.shrunk.ashrOrdered_3_genes.table, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/res.shrunk.ashrOrdered_3_genes.table.csv")
write.table(res.shrunk.normalOrdered_3_genes.table, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/res.shrunk.normalOrdered_3_genes.table.csv")
read.table("res.shrunkOrdered_3_genes.table.csv", header = T, sep = " ")
read.table("res.shrunk.ashrOrdered_3_genes.table.csv", header = T, sep = " ")
read.table("res.shrunk.normalOrdered_3_genes.table.csv", header = T, sep = " ")

summary(res.shrunkOrdered_3_genes.table) #set FDR as 0.1
res.shrunkOrdered_3_genes.table[, "log2FoldChange">= 1]
keepabs <- abs(res.shrunkOrdered_3_genes.table$log2FoldChange) >=1
keepabs.df <- as.data.frame(res.shrunkOrdered_3_genes.table[keepabs,])
write.table(keepabs.df, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/keepabs.df.csv")
read.table("keepabs.df.csv", header = T, sep = " ")



write.table(res.shrunkOrdered, file = "/Users/zunqiuwang/Dropbox/NYU_Bioinfo/NGS/WK9/res.shrunkOrdered.csv")  
res.shrunkOrdered <- read.table('res.shrunkOrdered.csv', header = T, sep = " ")
count.of.lfc_larger_zero <- count(res.shrunkOrdered %>% filter(log2FoldChange >0))
count.of.lfc_smaller_zero <- count(res.shrunkOrdered %>% filter(log2FoldChange <0))

count.df <- data.frame(cts = c(as.numeric(count.of.lfc_larger_zero), as.numeric(count.of.lfc_smaller_zero)), lfc= c('larger than 0', 'smaller than 0'))
count.df
ggplot(count.df, aes(x= lfc, y = cts, fill = lfc)) +
  geom_col()




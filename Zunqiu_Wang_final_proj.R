# create metadata files for tximport()

library(tximport)

sample_names <- c('SRR7819990','SRR7819991','SRR7819992','SRR7819993','SRR7819994','SRR7819995')
sample_condition <- c(rep('control',3),rep('treated',3))

files <- file.path("Salmon_quant",sample_names,paste(sample_names,".transcripts_quant",sep=""),'quant.sf')
names(files) <- sample_names

files %>% View()

tx2gene <- read.table("tx2gene.csv",header=F,sep=",")

tx2gene

all(file.exists(files))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

samples <- data.frame(sample_names=sample_names,condition=sample_condition)
row.names(samples) <- sample_names

glimpse(samples)

# Create DESeqDataSet obj

library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
class(ddsTxi)
ddsTxi

# remove genes with 10 or fewer reads

head(counts(ddsTxi))
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep, ]
ddsTxi

# rlog() and plotPCA()
rld <- rlog(ddsTxi)
library(ggplot2)
library(ggrepel)

                                  ### one outlier in treated? PC1 and PC2 both large +
                                           ### all seperates by condition on PC2

pca <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
png("final_proj_PCA.png")
ggplot(pca, aes(x=PC1, y=PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_text_repel(aes(label=condition)) +
  labs(title="PCA plot") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed()
dev.off()
                                           ### 2 seperates on PC1 1 outlier
# hcluster
dists <- dist(t(assay(rld)))
hclust(dists)
png("final_proj_hcluster.png")
plot(hclust(dists))
dev.off()   
                                          ### seems no batch effect, cluster according to condition      

### heatmap of first 20 genes    particular for samples whose experimental treatment suffered from an anormality that renders the data points obtained from these particular samples detrimental to our purpose.

library("pheatmap")
select <- order(rowMeans(counts(ddsTxi,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- data.frame(condition=colData(ddsTxi)[, "condition"])
rownames(df) <- colnames(assay(rld)[select,])
heatmap_pair <- pheatmap(cor(assay(rld)), annotation=df, main="heatmap of sample pair") 
heatmap <- pheatmap(assay(rld)[select,], annotation_col=df, main="heatmap of first 20 genes")
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()   
}                  # ref from https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
save_pheatmap_pdf(heatmap, "heatmap.pdf")
save_pheatmap_pdf(heatmap_pair, "heatmap_pair.pdf")

# run DeSeq
ddsTxi <- DESeq(ddsTxi)
class(ddsTxi) # Determine the type of object
ddsTxi

### save work in rds
saveRDS(ddsTxi, "ddsTxi.rds")

### read RDS
ddsTxi <- readRDS("ddsTxi.rds")

# dispersion-by-mean plot

estimateDispersions(ddsTxi)
png("disp_para.png")
plotDispEsts(ddsTxi)
dev.off()          #analysis on results

# Results
res <- results(ddsTxi, contrast = c('condition','control','treated')) ### control/treated
class(res) 
res
write.table(res, file = paste0(pwd, "/All_gene_res.csv"), sep = "\t", col.names = NA)

# lfcshrink, order by padj
res.shrunk <- lfcShrink(ddsTxi, contrast = c('condition','control','treated'), type = "ashr")
res.shrunkOrdered <- res.shrunk[order(res.shrunk$padj),]
head(res.shrunkOrdered, 10)

# plot(MA)
png("MA_shrunken.png")
plotMA(res.shrunk, ylim=c(-2,2)) #remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()                    ##lots of non biological significant lfc, 

# remove NA 
no.na.res.shrunk <- res.shrunk[!is.na(res.shrunk$padj),]
class(no.na.res.shrunk)
no.na.res.shrunk.df <- as.data.frame(no.na.res.shrunk)
class(no.na.res.shrunk.df)
str(no.na.res.shrunk.df)

# the number of biologically relevant differentially expressed genes
keeplog2fold <- abs(no.na.res.shrunk$log2FoldChange) >=1 #(2-fold)
keeplog2fold.df <- as.data.frame(no.na.res.shrunk[keeplog2fold,])
keeplog2fold.df 
nrow(keeplog2fold.df) #81 rows after removing NA in padj/pval
pwd <- getwd()
write.table(keeplog2fold.df, file = paste0(pwd, "/keeplog2fold.df.csv"), sep = "\t", col.names = NA)
read.table("keeplog2fold.df.csv", header = T, sep = "\t")
  
  
# the number of statistically significant genes at your chosen FDR(0.05)
keeppadj <- no.na.res.shrunk$padj < 0.05
res.shrunk_sig.df <- as.data.frame(no.na.res.shrunk[keeppadj, ])
res.shrunk_sig.df
nrow(res.shrunk_sig.df) #3062 rows
write.table(res.shrunk_sig.df, file = paste0(pwd, "/res.shrunk_sig.df.csv"), sep = "\t", col.names = NA)
read.table("res.shrunk_sig.df.csv", header = T, sep = "\t")

# a table of the number of biologically relevant differentially expressed genes and the number of statistically significant genes at your chosen FDR(0.05)
num_bio_sig <- nrow(keeplog2fold.df)
num_stats_sig <- nrow(res.shrunk_sig.df)
num.sig.df <- data.frame(num_bio_sig = num_bio_sig,
           num_stats_sig = num_stats_sig)
write.table(num.sig.df, file = paste0(pwd, "/num.sig.df.csv"), sep = "\t", row.names = FALSE)
read.table("num.sig.df.csv", header = T, sep = "\t")

# a table with the 10 most highly significant differentially expressed genes
library(dplyr)
no.na.res.shrunkOrdered <- no.na.res.shrunk[order(no.na.res.shrunk$padj),]
no.na.res.shrunkOrdered.df <- as.data.frame(no.na.res.shrunkOrdered)
no.na.res.shrunkOrdered.top.10.df <- no.na.res.shrunkOrdered.df %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1) %>% head(10)
write.table(no.na.res.shrunkOrdered.top.10.df, file = paste0(pwd, "/no.na.res.shrunkOrdered.top.10.df.csv"), sep = "\t", col.names = NA)
read.table("no.na.res.shrunkOrdered.top.10.df.csv", header = T, sep = "\t")
no.na.res.shrunkOrdered.filter.df <- no.na.res.shrunkOrdered.df %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) > 1) 
write.table(no.na.res.shrunkOrdered.filter.df, file = paste0(pwd, "/no.na.res.shrunkOrdered.filter.df.csv"), sep = "\t", col.names = NA)

#### the number of sig differentially expressed genes with lfc < -1 and lfc > 1 
no.na.res.shrunkOrdered.neg1.lfc.df <- no.na.res.shrunkOrdered.df %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1)
nrow(no.na.res.shrunkOrdered.neg1.lfc.df)  #68 rows/genes upregulated under treatment
no.na.res.shrunkOrdered.pos1.lfc.df <- no.na.res.shrunkOrdered.df %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1)
nrow(no.na.res.shrunkOrdered.pos1.lfc.df)    ##13 rows/genes downregulated under treatment

count.df <- data.frame(count = c(as.numeric(nrow(no.na.res.shrunkOrdered.pos1.lfc.df)), as.numeric(nrow(no.na.res.shrunkOrdered.neg1.lfc.df))), log2FoldChange= c('larger than 1', 'smaller than -1'))
write.table(count.df, file = paste0(pwd, "/count.biosig.genes.df.csv"), sep = "\t", row.names = FALSE)
read.table("count.biosig.genes.df.csv", header = T, sep = "\t")

png("Significant_differentially_expressed_genes.png")
ggplot(count.df, aes(x= log2FoldChange, y = count, fill = log2FoldChange)) +
  geom_col() +
  labs(title = "Significant differentially expressed genes") +
  geom_text(aes(label = count, vjust = -0.3)) +
  theme_bw()
dev.off()  

# DESeq2 results for ALL Genes as a tab-delimited file. ####including NA in padj/pval
res.shrunkOrdered.df <- as.data.frame(res.shrunkOrdered)
write.table(res.shrunkOrdered.df, file = paste0(pwd, "/res.shrunkOrdered.df.csv"), sep = "\t", col.names = NA)
read.table("res.shrunkOrdered.df.csv", header = T, sep = "\t")


# raw p value histogram: 
png("p_value_histogram_after_shrunk.png")
ggplot(as.data.frame(res.shrunk),aes(pvalue)) + geom_histogram(fill="light blue",color='black')
dev.off()                    ## an enrichment of low p-values. This is the expected result if there is a large class of differentially expressed genes between treatment and control.



# a table with the total number of reads and the mapping rate for each sample

fastqc_sequence_counts <- read.table("fastqc_sequence_counts_plot.tsv", header = T, sep = "\t")
fastqc_sequence_counts <- as.data.frame(fastqc_sequence_counts)
unique_reads <- as.numeric(fastqc_sequence_counts$Unique.Reads)
dup_reads <- as.numeric(fastqc_sequence_counts$Duplicate.Reads)
total_reads <- unique_reads + dup_reads
sample <- fastqc_sequence_counts$Category
mapping_rate <- c('91.2417%', '91.8589%', '92.8621%', '92.0663%', '92.3211%', '92.394%')
reads_map_rate.table <- data.frame(sample = sample,
           unique_reads = unique_reads,
           dup_reads = dup_reads,
           total_reads = total_reads,
           mapping_rate = mapping_rate)
reads_map_rate.table
write.table(reads_map_rate.table, file = paste0(pwd, "/reads_map_rate.table.csv"), sep = "\t", row.names = FALSE)
read.table("reads_map_rate.table.csv", header = T, sep = "\t")

# the number of statistically significant genes at your chosen FDR(0.05), the number of biologically relevant differentially expressed genes (defined using some criterion, such as a change in gene expression of two-fold or greater)

stat_sig_num <- as.numeric(nrow(res.shrunk_sig.df))
bio_sig_num <- as.numeric(nrow(keeplog2fold.df))
stat_bio_sig.table <- data.frame(stat_sig_num = stat_sig_num,
                                 bio_sig_num = bio_sig_num)
stat_bio_sig.table
write.table(stat_bio_sig.table, file = paste0(pwd, "/stat_bio_sig.table.csv"), sep = "\t", row.names = FALSE)
read.table("stat_bio_sig.table.csv", header = T, sep = "\t")


###metadata with gene symbol and gene_id
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

yyd <- toTable(org.Hs.egENSEMBL2EG)
xxd <- toTable(org.Hs.egSYMBOL)

gene_id_ensemble_symbol.df <- merge(yyd, xxd, by='gene_id')

no.na.res.shrunkOrdered.table <- no.na.res.shrunkOrdered.df %>% filter(padj < 0.05) %>% filter(abs(log2FoldChange) >=1)
no.na.res.shrunkOrdered.table$ensembl_id <- unlist(lapply(strsplit(rownames(no.na.res.shrunkOrdered.table), '\\.'),function(x){x[1]}))


no.na.res.shrunkOrdered.table <- merge(no.na.res.shrunkOrdered.table, gene_id_ensemble_symbol.df, by='ensembl_id')
metadata_symbol <- no.na.res.shrunkOrdered.table[order(no.na.res.shrunkOrdered.table$log2FoldChange),] 
write.table(metadata_symbol, file = paste0(pwd, "/metadata_symbol.csv"),sep = "\t", row.names = FALSE)
read.table("metadata_symbol.csv", header = T, sep = "\t") 

### volcano plot and add a new column indicating up/down regulation of DEG, named top 20 genes with padj< 0.05 and |lfc| > 1

no.na.res.shrunk.df$ensembl_id <- unlist(lapply(strsplit(rownames(no.na.res.shrunk.df), '\\.'),function(x){x[1]}))
no.na.res.shrunk.df <- merge(no.na.res.shrunk.df, gene_id_ensemble_symbol.df, by='ensembl_id')
no.na.res.shrunk.df$gene_id <- NULL
no.na.res.shrunk.df$diffexpressed <- "NO"
no.na.res.shrunk.df$diffexpressed[no.na.res.shrunk.df$log2FoldChange > 1 & no.na.res.shrunk.df$pvalue < 0.05] <- "UP"
no.na.res.shrunk.df$diffexpressed[no.na.res.shrunk.df$log2FoldChange < -1 & no.na.res.shrunk.df$pvalue < 0.05] <- "DOWN"

mycolors <- c("red", "green", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
no.na.res.shrunk.df$delabel <- NA
no.na.res.shrunk.df$delabel[no.na.res.shrunk.df$diffexpressed != "NO"] <- no.na.res.shrunk.df$symbol[no.na.res.shrunk.df$diffexpressed != "NO"]
ordered.shrunk.df <- no.na.res.shrunk.df[order(no.na.res.shrunk.df$padj),]
write.table(ordered.shrunk.df, file = paste0(pwd, "/ordered.shrunk.df.csv"),sep = "\t", row.names = FALSE)
ordered.filter.shrunk.df <- ordered.shrunk.df %>% filter(abs(log2FoldChange) > 1) %>% mutate(top20_symbol = "")
write.table(ordered.filter.shrunk.df, file = paste0(pwd, "/ordered.filter.shrunk.df.csv"),sep = "\t", row.names = FALSE)
ordered.filter.shrunk.df$top20_symbol[1:20] <- as.character(ordered.filter.shrunk.df$delabel[1:20])
options(ggrepel.max.overlaps = Inf)
png("volcano_plot_only_DE.png")
ggplot(data=ordered.filter.shrunk.df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-1, 1), col="orange") +
  geom_hline(yintercept=-log10(0.05), col="orange") +
  scale_color_manual(values = mycolors) +
  labs(col = "DE genes?") +
  geom_text_repel(aes(label=top20_symbol))
dev.off()

png("volcano_plot.png")
ggplot(data=ordered.shrunk.df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept=c(-1, 1), col="orange") +
  geom_hline(yintercept=-log10(0.05), col="orange") +
  scale_color_manual(values = mycolors) +
  labs(col = "DE genes?")
dev.off()


install.packages("gt")
library(gt)
BiocManager::install("tximport")
library(tximport)

sample_names <- c('PDAC253','PDAC282','PDAC286','PDAC316','PDAC266','PDAC273','PDAC306','PDAC318')
sample_condition <- c(rep('highSucrose',4),rep('lowSucrose',4))

files <- file.path(sample_names,paste(sample_names,".transcripts_quant",sep=""),'quant.sf')
names(files) <- sample_names

files %>% View()

tx2gene <- read.table("Pdac_Barhee_chr_unan_180126_maker_HC.tx2gene",header=F,sep=",")

tx2gene

all(file.exists(files))

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

samples <- data.frame(sample_names=sample_names,condition=sample_condition)
row.names(samples) <- sample_names

# create DastaSet
library("DESeq2")
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
class(ddsTxi)

ddsTxi

assays(ddsTxi)
ddsTxi

# remove genes with 10 or fewer reads
head(counts(ddsTxi))
keep <- rowSums(counts(ddsTxi)) >= 10
ddsTxi <- ddsTxi[keep, ]
ddsTxi #check genes left

# run DeSeq
ddsTxi <- DESeq(ddsTxi)
class(ddsTxi) # Determine the type of object
ddsTxi

# results()
res <- results(ddsTxi, contrast = c('condition','lowSucrose','highSucrose') )
class(res) 
head(res)
resOrdered <- res[order(res$padj), ]
head(resOrdered, 10)

# lfcshrink
res.shrunk <- lfcShrink(ddsTxi, contrast = c('condition','lowSucrose','highSucrose'), type = "ashr")
res.shrunkOrdered <- res.shrunk[order(res.shrunk$padj),]
head(res.shrunkOrdered, 10)

# visualize p value after lfcshrink

png("p_value_histogram_after_shrunk.png")
ggplot(as.data.frame(res.shrunk),aes(pvalue)) + geom_histogram(fill="light blue",color='black')
dev.off()

# report top 10 genes according to padj
res.shrunkOrdered <- res.shrunk[order(res.shrunk$padj),]
res.shrunkOrdered_top_10 <- head(res.shrunkOrdered, 10)
pwd <- getwd()
write.table(res.shrunkOrdered_top_10, file = paste0(pwd, "/res.shrunkOrdered_top_10.csv"))
read.table("res.shrunkOrdered_top_10.csv", header = T, sep = " ")

res.shrunkOrdered_df_genes_interest <- res.shrunkOrdered[ row.names(res.shrunkOrdered) %in% c('Pdac_HC_chr14G0022900','Pdac_HC_chr14G0023100','Pdac_HC_chr14G0028200'), ]
res.shrunkOrdered_df_genes_interest
write.table(res.shrunkOrdered_df_genes_interest, file = paste0(pwd, "/res.shrunkOrdered_df_genes_interest.csv"))
read.table("res.shrunkOrdered_df_genes_interest.csv", header = T, sep = " ")


#top10 genes
library(gt)
library(dplyr)

tb1 <- read.table("res.shrunkOrdered_top_10.csv", header = T, sep = " ")
tb1 %>% gt() %>% tab_header(title="top 10 differentially expressed genes")

#3 genes of interests
tb2 <- read.table("res.shrunkOrdered_df_genes_interest.csv", header = T, sep = " ") 
tb2 %>% gt() %>% tab_header(title="3 DEG")


# estimateDispersions() and plotDispEsts()
?estimateDispersions
?plotDispEsts
# para type
estimateDispersions(ddsTxi, fitType= "parametric")
png("disp_para.png")
plotDispEsts(ddsTxi)
dev.off()

#local type
estimateDispersions(ddsTxi, fitType= "local")
png("disp_local.png")
plotDispEsts(ddsTxi)
dev.off()

#mean type
estimateDispersions(ddsTxi, fitType= "mean")
png("disp_mean.png")
plotDispEsts(ddsTxi)
dev.off()







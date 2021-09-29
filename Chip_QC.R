SRR7207017.df <- read.table("SRR7207017_peaks.narrowPeak", header = F, sep = "\t", stringsAsFactors = F)
SRR7207017.df$width <- SRR7207017.df$V3-SRR7207017.df$V2
SRR7207017.mean <- mean(SRR7207017.df$width)
SRR7207017.mean

SRR7207011.df <- read.table("SRR7207011_peaks.narrowPeak", header = F, sep = "\t", stringsAsFactors = F)
SRR7207011.df$width <- SRR7207011.df$V3-SRR7207011.df$V2
SRR7207011.mean <- mean(SRR7207011.df$width)
SRR7207011.mean

69002/22517898
877004/20702991



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIC")

system("wget https://www.encodeproject.org/files/ENCFF000BDQ/@@download/ENCFF000BDQ.bam")
system("wget https://www.encodeproject.org/files/ENCFF000BFX/@@download/ENCFF000BFX.bam")

library(ChIC)

#  Set the working directory to the location where the ENCFF000BDQ.bam and ENCFF000BFX.bam are located:
setwd(getwd())

# Now define two variables with the prefix for the two BAM filenames (required by ChIC) and read the BAM files.
chipName <- "ENCFF000BFX"
inputName <- "ENCFF000BDQ"

chipBam <- readBamFile(chipName)
inputBam <- readBamFile(inputName)

# Define variables which we will use to run the analysis with multiple threads below.
mc <- 4
cluster <- parallel::makeCluster( mc )

# Calculate QC-metrics from CrossCorrelation analysis: 
# cross-correlation profiling with ChIC is to determine binding peak separation distance 
# and approximate window size that should be used for binding detection

chip_binding.characteristics <- get.binding.characteristics(
  chipBam,
  srange=c(0,500),
  bin = 5,
  accept.all.tags = TRUE,
  cluster = cluster)
input_binding.characteristics <- get.binding.characteristics(
  inputBam,
  srange=c(0,500),
  bin = 5,
  accept.all.tags = TRUE,
  cluster = cluster)
parallel::stopCluster( cluster )

class(chip_binding.characteristics)
summary(chip_binding.characteristics$cross.correlation)
str(chip_binding.characteristics$cross.correlation)
chip_binding.characteristics

# Now using the output binding characteristics from above, we can calculate the cross correlation QC-metrics (e.g., NSC and RSC) 
# that are used to measure signal-to-noise for the ChIP sample and generate the cross-correlation plot. Execution of the command below should create a pdf in your working directory called CrossCorrelation.pdf.

crossvalues_Chip <- getCrossCorrelationScores( chipBam ,
                                               chip_binding.characteristics,
                                               read_length = 36,
                                               annotationID="hg19",
                                               savePlotPath = ".",
                                               mc = mc)
class(crossvalues_Chip)

crossvalues_Chip$CC_NSC
crossvalues_Chip$CC_RSC

# getCrossCorrelationScores function for the control (“input”) sample and generate a cross-correlation profile plot.

crossvalues_Input <- getCrossCorrelationScores( inputBam ,
                                               input_binding.characteristics,
                                               read_length = 36,
                                               annotationID="hg19",
                                               savePlotPath = ".",
                                               mc = mc)
class(crossvalues_Input)

crossvalues_Input$CC_NSC
crossvalues_Input$CC_RSC


# We first remove anomalies in the BAMs as follows (see vignette)

selectedTags <- removeLocalTagAnomalies(chipBam,
                                        inputBam,
                                        chip_binding.characteristics,
                                        input_binding.characteristics)
inputBamSelected <- selectedTags$input.dataSelected
chipBamSelected <- selectedTags$chip.dataSelected

# extract the “tag shift” value from the cross-correlation output.

finalTagShift <- crossvalues_Chip$tag.shift

# calculate smoothed tag (read) density distributions. 

smoothedChip <- tagDensity(chipBamSelected,
                           annotationID = "hg19",
                           tag.shift = finalTagShift, mc = mc)

smoothedInput <- tagDensity(inputBamSelected,
                            annotationID = "hg19",
                            tag.shift = finalTagShift, mc = mc)

# fingerprint plot and calculate global entrichment metrics

Ch_Results <- qualityScores_GM(densityChip = smoothedChip,
                               densityInput = smoothedInput,
                               savePlotPath = ".")

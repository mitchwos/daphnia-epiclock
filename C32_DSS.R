library(DSS)
library(readr)

setwd("~/Documents/IRP/R/clock/C32")

# Read in all files
file.list = list.files(("./"),pattern="*_positions.txt")
samples <- lapply(file.list, function(f){
  read_delim(f, "\t", escape_double = F, col_names = T, trim_ws = T)
})

# make BSobj 
BSobj = makeBSseqData(samples,
                      c("C32_HF_old_1", "C32_HF_old_2", "C32_HF_old_3",
                        "C32_HF_old_4", "C32_HF_old_5", "C32_HF_old_6",
                        "C32_HF_young_1", "C32_HF_young_2", "C32_HF_young_3",
                        "C32_HF_young_4", "C32_HF_young_5", "C32_HF_young_6",
                        "C32_LF_old_1", "C32_LF_old_2", "C32_LF_old_3",
                        "C32_LF_old_4", "C32_LF_old_5", "C32_LF_old_6",
                        "C32_LF_young_1", "C32_LF_young_2", "C32_LF_young_3",
                        "C32_LF_young_4", "C32_LF_young_5", "C32_LF_young_6"))

# Make the metadata
Age = base::as.factor(c(rep("OLD",6),rep("YOUNG",6),rep("OLD",6),rep("YOUNG",6)))
Diet <- c(rep("HF",12),rep("LF",12))
design = data.frame(Age, Diet)
design$Age <- as.numeric(design$Age)

# Define model
mod = model.matrix(~Age+ Diet, design)

# Run the linear model
DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Age+Diet + Age:Diet)

# See the names of the coefficients
colnames(DMLfit$X)

#Turn off Scientific Notation
options(scipen = 999) 

#Age Differences

# Pull out the result of interest - this gives all the individual CpGs with pvals
Age = DMLtest.multiFactor(DMLfit, coef="Age")

Significant_Age<- subset(Age , Age$fdrs < 0.05)
nrow(Significant_Age)
write.table(Significant_Age, "C32_Significant_Age_CpG.txt", quote = F, sep = "\t")

# Filter this result for regions with many significant CpGs (this is actually quite nice)
significant_age_pos <- callDMR(Age, p.threshold=0.05)
nrow(significant_age_pos)
write.table(significant_age_pos, "C32_Significant_Age_Regions.txt", quote = F, sep = "\t")

#---------------------------------------------
# Extract methylated counts and total coverage
methylated_counts <- getCoverage(BSobj, type = "M")
total_counts <- getCoverage(BSobj, type = "Cov")

# Calculate percentage methylation
percent_meth <- (methylated_counts / total_counts) * 100
percent_methylation <- as.data.frame(percent_meth)

#Extract seqnames and start positions
seqnames_values <- as.character(seqnames(SummarizedExperiment::rowRanges(BSobj)))
start_positions <- start(SummarizedExperiment::rowRanges(BSobj))

percent_methylation <-cbind(percent_methylation,scaf = seqnames_values, start=start_positions)

# Combine 'scaf' and 'start' into a new column called 'Row_ID'
percent_methylation$Row_ID <- paste(percent_methylation$scaf, percent_methylation$start, sep = ".")

# Remove the 'scaf' and 'start' columns
percent_methylation <- percent_methylation[, !(colnames(percent_methylation) %in% c("scaf", "start"))]

fwrite(percent_methylation, file = "~/Documents/IRP/R/clock/C32/C32_daphnia_percent")

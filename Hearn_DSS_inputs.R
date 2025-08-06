# ------------------------------------------------------------------------
# Making DSS input files
# ------------------------------------------------------------------------

# NOTE: run this on ALICE - 5hrs, 80GB memory.
library(methylKit)
setwd("/scratch/evo-epi/mjw85/methyl_data")

## -------------------------------------------------------------------------
sample.list <- list("Hearn/ERR3535037_C32_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535038_C32_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535039_C32_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535040_C32_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535041_C32_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535042_C32_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535043_C32_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535044_C32_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535045_C32_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535046_C32_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535047_C32_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535048_C32_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535049_C32_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535050_C32_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535051_C32_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535052_C32_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535053_C32_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535054_C32_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535055_C32_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535056_C32_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535057_C32_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535058_C32_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535059_C32_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535060_C32_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535061_KA53_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535062_KA53_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535063_KA53_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535064_KA53_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535065_KA53_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535066_KA53_HF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535067_KA53_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535068_KA53_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535069_KA53_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535070_KA53_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535071_KA53_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535072_KA53_LF_OLD_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535073_KA53_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535074_KA53_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535075_KA53_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535076_KA53_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535077_KA53_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535078_KA53_HF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535079_KA53_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535080_KA53_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535081_KA53_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535082_KA53_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535083_KA53_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz",
                    "Hearn/ERR3535084_KA53_LF_YOUNG_trim.CpG_report.merged_CpG_evidence.cov.gz")

# Read in data
CPGRaw <- methRead(sample.list, 
                   sample.id = list("C32_HF_old_1", "C32_HF_old_2", "C32_HF_old_3",
                                    "C32_HF_old_4", "C32_HF_old_5", "C32_HF_old_6",
                                    "C32_LF_old_1", "C32_LF_old_2", "C32_LF_old_3",
                                    "C32_LF_old_4", "C32_LF_old_5", "C32_LF_old_6",
                                    "C32_HF_young_1", "C32_HF_young_2", "C32_HF_young_3",
                                    "C32_HF_young_4", "C32_HF_young_5", "C32_HF_young_6",
                                    "C32_LF_young_1", "C32_LF_young_2", "C32_LF_young_3",
                                    "C32_LF_young_4", "C32_LF_young_5", "C32_LF_young_6",
                                    "KA53_HF_old_1", "KA53_HF_old_2", "KA53_HF_old_3",
                                    "KA53_HF_old_4", "KA53_HF_old_5", "KA53_HF_old_6",
                                    "KA53_LF_old_1", "KA53_LF_old_2", "KA53_LF_old_3",
                                    "KA53_LF_old_4", "KA53_LF_old_5", "KA53_LF_old_6",
                                    "KA53_HF_young_1", "KA53_HF_young_2", "KA53_HF_young_3",
                                    "KA53_HF_young_4", "KA53_HF_young_5", "KA53_HF_young_6",
                                    "KA53_LF_young_1", "KA53_LF_young_2", "KA53_LF_young_3",
                                    "KA53_LF_young_4", "KA53_LF_young_5", "KA53_LF_young_6"),
                   assembly="daphnia",
                   treatment=c(0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,
                               5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7),
                   context="CpG",
                   dbtype = NA,
                   pipeline = "bismarkCoverage",
                   header = T, sep = "\t", mincov=1,
                   dbdir = '/scratch/evo-epi/mjw85/methyl_data')

# Filter on at least 10 coverage per cpg
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
rm(CPGRaw)

# Keep only CpGs which have coverage in at least two replicates per condition
meth_all_data <- unite(filtered_data, destrand=F, min.per.group = 2L)
rm(filtered_data)

# How many CpGs have enough coverage?
nrow(meth_all_data) # 2224847

# Get data from methylobject
df_meth_all <- getData(meth_all_data)

# Remove first few uninformative columns
df_meth_all_1 <- df_meth_all[,-c(1:4)]

# Subset data into individual dataframes per sample
all_data <- lapply(seq(1, ncol(df_meth_all_1), by=3), function(i) 
  df_meth_all_1[i: pmin((i+1), ncol(df_meth_all_1))])

# Define binomial test
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}

# Set up output dataframe
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

# Run binomial 
for (i in seq_along(all_data)) {
  all_data[[i]] <- all_data[[i]][complete.cases(all_data[[i]]),]
  colnames(all_data[[i]]) <- c("CT", "Ccount")
  all_data[[i]]$pVal <- mapply(bt, all_data[[i]]$Ccount, all_data[[i]]$CT)
  all_data[[i]]$FDR <- p.adjust(all_data[[i]]$pVal, method = "BH", n = length(all_data[[i]]$pVal))
  all_data[[i]]$row <- as.numeric(rownames(all_data[[i]]))
  dfmeth <- subset(all_data[[i]], all_data[[i]]$FDR < 0.05)
  allrows <- rbind(allrows, dfmeth)
}

# Get positions which are methylated in at least one sample
meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions)

# Pull out dataframes again and subset to keep only these positions
all_data <- lapply(seq(1, ncol(df_meth_all_1), by=3), function(i) 
  df_meth_all_1[i: pmin((i+1), ncol(df_meth_all_1))])

# Add sample names to dataframes in the list
sample.id = list("C32_HF_old_1", "C32_HF_old_2", "C32_HF_old_3",
                 "C32_HF_old_4", "C32_HF_old_5", "C32_HF_old_6",
                 "C32_LF_old_1", "C32_LF_old_2", "C32_LF_old_3",
                 "C32_LF_old_4", "C32_LF_old_5", "C32_LF_old_6",
                 "C32_HF_young_1", "C32_HF_young_2", "C32_HF_young_3",
                 "C32_HF_young_4", "C32_HF_young_5", "C32_HF_young_6",
                 "C32_LF_young_1", "C32_LF_young_2", "C32_LF_young_3",
                 "C32_LF_young_4", "C32_LF_young_5", "C32_LF_young_6",
                 "KA53_HF_old_1", "KA53_HF_old_2", "KA53_HF_old_3",
                 "KA53_HF_old_4", "KA53_HF_old_5", "KA53_HF_old_6",
                 "KA53_LF_old_1", "KA53_LF_old_2", "KA53_LF_old_3",
                 "KA53_LF_old_4", "KA53_LF_old_5", "KA53_LF_old_6",
                 "KA53_HF_young_1", "KA53_HF_young_2", "KA53_HF_young_3",
                 "KA53_HF_young_4", "KA53_HF_young_5", "KA53_HF_young_6",
                 "KA53_LF_young_1", "KA53_LF_young_2", "KA53_LF_young_3",
                 "KA53_LF_young_4", "KA53_LF_young_5", "KA53_LF_young_6")
names(all_data) <- sample.id

# Add back in chr and position and filter
chr_pos <- df_meth_all[,c(1,2)]

for( i in seq_along(all_data)){
  all_data[[i]] <- cbind(chr_pos, all_data[[i]])
  all_data[[i]] <- subset(all_data[[i]], row.names(all_data[[i]]) %in% meth_positions)
  colnames(all_data[[i]]) <- c("chr","pos", "N", "X")
  myfile <- file.path("./", paste0(names(all_data[i]),"_","positions.txt"))
  write.table(all_data[[i]], file=myfile, quote=F, sep="\t", row.names=F)
}


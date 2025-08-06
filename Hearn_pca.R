## -------------------------------------------------------------------------
## Making PCA plot for all samples together
## -------------------------------------------------------------------------

# NOTE: run this on a cluster - 5hrs, 80GB memory.
# If using less samples it will be fine on your own computer
library(grid)
library(readr)
library(ggplot2)
library(methylKit)
library(ggdendro)
library(ggfortify)
library(ggrepel)

## -------------------------------------------------------------------------
# Define samples
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
                   dbdir= path)

# Filter on at least 10 coverage per cpg
filtered_data <- filterByCoverage(CPGRaw,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)
rm(CPGRaw)

# Normalise by library size
normalized <- normalizeCoverage(filtered_data)
rm(filtered_data)

# Keep only CpGs which have coverage in at least two replicates per condition
meth_all_data <- unite(normalized, destrand=F, min.per.group = 2L)
rm(normalized)

# How many CpGs have enough coverage?
nrow(meth_all_data) # 2224847

# get data from methylobject
df_meth_all <- getData(meth_all_data)

# remove first few uninformative columns
df_meth_all <- df_meth_all[,-c(1:4)]

# subset data into individual dataframes per sample
all_data <- lapply(seq(1, ncol(df_meth_all), by=3), function(i) 
  df_meth_all[i: pmin((i+1), ncol(df_meth_all))])

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
length(meth_positions) # 81282

# Keep only these positions for diff meth test
subset_methBase <- methylKit::select(meth_all_data, meth_positions)
head(subset_methBase)

# Write out data for later use
subset_methBase_out <- getData(subset_methBase)
write.table(subset_methBase_out, file="sites_for_PCA.txt", sep="\t", quote=F, row.names = F, col.names = T)

## -------------------------------------------------------------------------
# Make a nice PCA

# Run PCA analysis
PCA_data <- PCASamples(subset_methBase, screeplot=F, obj.return = T)

# Pull out dataframe
PCA_data1 <- as.data.frame(PCA_data$x)

# Add sample metadata
PCA_data1$sample <- row.names(PCA_data1)
PCA_data1$Age <- c(rep("50",12),rep("10",12),rep("50",12),rep("10",12))
PCA_data1$genotype <- c(rep('C32',24),rep('KA53',24))
PCA_data1$age_clonal_line <- c(rep("C32 - 50",12),rep("C32 - 10",12),rep("KA53 - 50",12),rep("KA53 - 10",12))
PCA_data1$Diet <- c(rep('HF',6),rep('LF',6),rep('HF',6),rep('LF',6),rep('HF',6),rep('LF',6),rep('HF',6),rep('LF',6))

# Write out PCA data to mess with later
write.table(PCA_data1, file="PCA_data.txt", sep ="\t", quote=F, row.names = F, col.names = T)

# Get percentages for PCA axis
percentage <- round(PCA_data$sdev / sum(PCA_data$sdev) * 100, 0)
percentage <- paste(colnames(PCA_data), "(", paste( as.character(percentage), "%", ")", sep="") )

# Write out percentages to mess with later
write.table(percentage, file="percentages_all_PCAs.txt", sep ="\t", quote=F, row.names = F, col.names = T)

# Plot both C32 and KA53
ggplot(PCA_data1, aes(PC1, PC2, colour=age_clonal_line))+
  geom_point(size=10)+
  scale_color_manual(values=c('#94DEFF','#00ADFF','#FFAF8A','#F56F31')) +
  xlab(paste0("PC1:",percentage[1]," variance")) +
  ylab(paste0("PC2:",percentage[2]," variance")) +
  labs(
    title = 'C32 / KA53 PCA',
    subtitle = 'All Samples',
    x = 'PC1 : 9% Variance',
    y = 'PC2 : 5% Variance',
    colour = 'Genotype - Age'
    ) +
  theme_bw() +
  theme(
    plot.title=element_text(face='bold',size=16,hjust=0.5),
    plot.subtitle=element_text(size=12,hjust=0.5),
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12,hjust=0.5))

# Create C32 only df
C32 <- PCA_data1[PCA_data1$genotype == 'C32', ]

# Plot C32
ggplot(C32, aes(PC1, PC2, colour=Age, shape=Diet, stroke=Diet, size=Diet))+
  geom_point() +
  scale_size_manual(values=c(10,8)) +
  scale_color_manual(values=c('#94DEFF','#00ADFF')) +
  scale_shape_manual(values=c('HF'=16,'LF'=1)) +
  scale_discrete_manual(
    aesthetics = 'stroke',
    values=c('HF'=0,'LF'=2)) +
  labs(
    title = 'C32 PCA',
    subtitle = 'All Samples',
    x = 'PC1 : 9% Variance',
    y = 'PC2 : 5% Variance',
    colour = 'Age',
  ) +
  theme_bw() +
  theme(
    plot.title=element_text(face='bold',size=16,hjust=0.5),
    plot.subtitle=element_text(size=12,hjust=0.5),
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12,hjust=0.5)) +
  guides(colour=guide_legend(override.aes = list(size=10)))

# Create KA53 only df
KA53 <- PCA_data1[PCA_data1$genotype == 'KA53', ]

# Plot KA53
ggplot(KA53, aes(PC1, PC2, colour=Age, shape=Diet, stroke=Diet, size=Diet))+
  geom_point() +
  scale_size_manual(values=c(10,8)) +
  scale_color_manual(values=c('#FFAF8A','#F56F31')) +
  scale_shape_manual(values=c('HF'=16,'LF'=1)) +
  scale_discrete_manual(
    aesthetics = 'stroke',
    values=c('HF'=0,'LF'=2)) +
  labs(
    title = 'KA53 PCA',
    subtitle = 'All Samples',
    x = 'PC1 : 9% Variance',
    y = 'PC2 : 5% Variance',
  ) +
  theme_bw() +
  theme(
    plot.title=element_text(face='bold',size=16,hjust=0.5),
    plot.subtitle=element_text(size=12,hjust=0.5),
    axis.text=element_text(size=12),
    axis.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.title=element_text(size=12,hjust=0.5)) +
  guides(colour=guide_legend(override.aes = list(size=10)))

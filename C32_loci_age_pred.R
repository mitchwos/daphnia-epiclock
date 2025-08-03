library(dplyr)
library(ggplot2)

setwd('~/Documents/IRP/R/clock/Hearn_meth/C32')

# Load clock loci 
clock_loci <- read.csv('~/Documents/IRP/R/clock/clock_loci.csv', header = FALSE)
intercept <- as.numeric(clock_loci[2,2]) # Extract intercept
clock_loci <- clock_loci[3:44,1:2] 
colnames(clock_loci) <- c('cpg','m') 

# Load coverage files
methyl_files = list.files(("./"),pattern="*.cov")

# Create empty ages df
ages <- data.frame(sample = character(0), age = numeric(0), stringsAsFactors = FALSE)

# Predict ages for each sample
for (f in methyl_files) {
  cov <- read.delim(f, header = FALSE)
  colnames(cov) <- c('chr','pos','end','meth_per','meth_count','unmeth_count')
  cov$cpg <- paste0(cov$chr, ".", cov$pos) # Format to match clock loci
  cov <- cov %>% select(cpg, meth_per) # Keep only cpg and meth_per cols
  matched <- left_join(clock_loci, cov, by = 'cpg') # Merge dfs matched by common loci
  matched$meth_per[is.na(matched$meth_per)] <- 0 # Convert na values to zero
  matched$meth_per <- matched$meth_per /100 # Convert methylation percentage to decimal
  
  matched$meth_per <- as.numeric(matched$meth_per)
  matched$m <- as.numeric(matched$m)
  
  matched$age <- matched$meth_per * matched$m # Multiply methylation percentages by loci values
  age <- sum(matched$age) + intercept # Sum up age col and intercept to predict age
  ages <- rbind(ages, data.frame(sample = f, age = age, stringsAsFactors = FALSE)) # Add sample and predicted age to ages df
}

# Clean up sample names
ages$sample <- gsub("_trim.*$", "", ages$sample)

write.csv(ages, "~/Documents/IRP/R/clock/C32_ages.csv", row.names = FALSE)

# Create age factor
ages$group <- ifelse(grepl("OLD$", ages$sample), "OLD", "YOUNG")

# Boxplot comparison of OLD and YOUNG predicted ages
ggplot(ages, aes(x = group, y = age, fill = group)) +
  geom_boxplot() +
  geom_jitter(width=0.2, size=2, alpha=0.5) +
  labs(title = "C32 - Predicted Age with Clock Loci",
       x = "Group",
       y = "Predicted Age") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title = element_text(size = 14),  # Increase axis title font size
    axis.text = element_text(size = 12),   # Increase axis text font size
    legend.title = element_text(size = 14),  # Adjust legend title font size
    legend.text = element_text(size = 12),   # Adjust legend text font size
    legend.position = "top"  # Move legend to the top
  )

ggsave("~/Documents/IRP/R/clock/Hearn_meth/C32/C32_clock_loci.pdf")


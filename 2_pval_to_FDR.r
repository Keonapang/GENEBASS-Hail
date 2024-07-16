
#!/usr/bin/env Rscript

# Script:/UKBIOBANK/pangk/Genebass_results/2_pval_to_FDR.r

# Input (argument): one /2_pval_filtered/500k_trait_pval.txt.gz file containing PLINK_SNP_NAME and trait_pval
# Output: to directory /UKBIOBANK/pangk/Genebass_results/3_FDR
        # PLINK_SNP_NAME height_FDR
        # 1:925920:A:C 0.771294151937272

# tool: qvalue package in R
###########################################################################################

args <- commandArgs(trailingOnly = TRUE)

library(qvalue)

# Get the trait from the filename
filename <- args[1] # 500k_trait_pval.txt.gz file 
print(paste("Loading file:", filename))

# get trait name 
trait <- sub("^500k_", "", basename(filename)) # Remove prefix
trait <- sub("_pval\\.txt\\.gz$", "", trait) # Remove suffix
trait <- paste0(trait, "_FDR")
print(trait)

# Read the input file to get the p-value column

x <- read.table(gzfile(filename), header=TRUE, sep=" ")
pval_col <- as.vector(x[, 2])

# Calculate FDR using Benjamini-Hochberg procedure
print("Calculating FDR using Benjamini-Hochberg procedure")
FDR_col <- p.adjust(pval_col, method = "BH")
new_df <- data.frame(x[, 1], FDR_col)

# Create a new dataframe with the first column of 'x' and 'adjusted_pvalues'
colnames(new_df) <- c("PLINK_SNP_NAME", trait)

# Print the dimensions of new_df
print(paste("Output file consists of columns PLINK_SNP_NAME trait_FDR:", dim(new_df)))

# Write to directory 3 within genebass
output_file <- paste("/UKBIOBANK/pangk/Genebass_results/3_FDR/", trait, "_7m.txt", sep = "")
write.table(new_df, output_file, sep = " ", row.names = FALSE, quote = FALSE)
print(paste("Output written to:", output_file))
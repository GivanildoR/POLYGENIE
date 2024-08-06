library(vcfR)
library(adegenet)
library(ggplot2)

# Read the VCF file
vcf <- read.vcfR("yourfile.vcf")

# Function to calculate nucleotide diversity (π) for polyploids using vcfR
calc_nucleotide_diversity_vcfR <- function(vcf, ploidy) {
  geno <- extract.gt(vcf, return.alleles = TRUE)
  total_diff <- 0
  total_comparisons <- 0
  
  for (i in 1:ncol(geno)) {
    locus_genotypes <- geno[, i]
    alleles <- unlist(strsplit(locus_genotypes, split = "[|/]"))
    
    # Calculate pairwise differences
    for (j in 1:(length(alleles) - 1)) {
      for (k in (j + 1):length(alleles)) {
        total_diff <- total_diff + (alleles[j] != alleles[k])
      }
    }
    
    # Calculate the number of comparisons for this locus
    total_comparisons <- total_comparisons + (length(alleles) * (length(alleles) - 1) / 2)
  }
  
  # Calculate nucleotide diversity (π)
  pi <- total_diff / total_comparisons
  return(pi)
}

# Calculate π for the dataset
ploidy_level <- 4  # Adjust this based on your polyploid level
pi_value <- calc_nucleotide_diversity_vcfR(vcf, ploidy_level)
print(pi_value)


# Plot nucleotide diversity (π) across regions
ggplot(pi_data, aes(x = Region, y = Pi, fill = Region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nucleotide Diversity (π) Across Regions",
       x = "Region",
       y = "Nucleotide Diversity (π)") +
  theme_minimal()

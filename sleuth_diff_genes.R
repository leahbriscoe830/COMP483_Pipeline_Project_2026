# Find the differentially expressed genes between samples 2 and 6 dpi using sleuth

library(sleuth)
library(dplyr)

# Read in sample, condition, path table
stab = read.table(snakemake@input[[1]], header=TRUE)

# Fit the model
so = sleuth_prep(stab)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

# Extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 

# Filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval < 0.05) |> dplyr::arrange(qval) 

# Subset only the target id, test statistic, and the p & q values
sleuth_stats = subset(sleuth_significant, select = c(target_id, test_stat, pval, qval))

# Write FDR < 0.05 transcripts to temp txt file
write.table(sleuth_stats, file=snakemake@output[[1]], append=TRUE, quote = FALSE, row.names = FALSE, sep = "\t")
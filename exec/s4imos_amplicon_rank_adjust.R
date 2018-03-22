##==================================================================================================================
#
# When simulate data, primers used for amplicon simulation may not cover all the genomes used in the simulation.
# Through this script, you can keep the current relative abundances from shotgun simulation and the script is going
# to correct the abundances in a proportional way for the covered genomes.
# 
# Before, use s4imos_amplicon_search.pl script and select the covered genomes by primers.
# 
# This script can be used by an Rscript call as:
#   Rscript s4imos_amplicon_rank_adjust.R --input [INPUT] --out [OUTPUT]
# 
# The input file should be a TSV (without columns names):
# 
#   taxonomyId | relativeAbundance
# 
# The output file follows the same structure, differing the values with corrected relative abundances.
#
##==================================================================================================================


##==================================================================================================================
## Required libraries
##==================================================================================================================
require(optparse)


##==================================================================================================================
## Parameters settings
##==================================================================================================================
option_list = list(
  make_option(c("-i", "--input"), type="character", default = "amplicon_ranks.tsv",
              help="Input file with amplicon ranks - first column as ID, second column as RELATIVE ABUNDANCE", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="adjusted_amplicon_ranks.tsv",
              help="Output file with the amplicons adjusted ranks", metavar="character")
); 


##==================================================================================================================
## Parameters parsing
##==================================================================================================================
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


##==================================================================================================================
## Data loading
##==================================================================================================================
ranks.filename <- opt$input
ranks.df <- read.delim(ranks.filename, header = F, col.names = c("id", "relative_abundance"))


##==================================================================================================================
## Data normalization
##==================================================================================================================
normalization.factor <- 100/sum(ranks.df$relative_abundance)
ranks.df$relative_abundance <- ranks.df$relative_abundance * normalization.factor


##==================================================================================================================
## Data exporting
##==================================================================================================================
write.table(ranks.df, opt$out, quote = F, row.names = F, col.names = F, sep = "\t")

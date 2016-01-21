library(dplyr)
source("../R/mutation_null_model.R")

maps_fit = function(synonymous_vars, mu_snp){
  # synonymous vars should have chr, pos, ref, alt, allele_count, tri_nuc_context
  # mu snps is three columns: from, to, mu_snp

  synonymous_tri = ddply(neutral_vars, c('tri_nuc_context', 'ref', 'alt'), function(x) {
    data.frame(n=nrow(x), singletons=sum(x$allele_count == 1))})

  synonymous_tri = merge(synonymous_tri, mu_snp, by.x = "tri_nuc_context", by.y = "from")

  expected_proportion_singleton_lm = lm(singletons/n ~ mu_snp, synonymous_tri, weights = synonymous_tri$n)

  return(rates)
}


MAPS = function(variants, split_factor, maps_model){
  # take a dataframe of variants that are separated by some identifier called split_factor
  # maps_model is a linear model learned on presumed synonymous mutations using lm(singletons/n ~ mu_snp)
  # with mu_snp from the tri-nucleotide mutation model

  # calculate the singleton_ratio_raw for each split_factor

  # calculate the mu_snp

}

# get all gencode seqs
#gencode = read.table("~/phd/code/CNE/data/gencode_protein_coding_genes_v19_+strand.txt", header = TRUE, sep = "\t")
#library(BSgenome.Hsapiens.UCSC.hg19)
#takes ~30 minutes on 16Gb RAM mac - but saves lots of time later.
#seqs = mapply(function(c,s,e) getSeq(Hsapiens, c, s, e), as.character(gencode$chr), gencode$start, gencode$stop)
#s = sapply(seqs, as.character)
#gencode$seq = s
#write.table(gencode, file = "../data/gencode_protein_coding_genes_v19_+strand.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

gencode = read.table("../data/gencode_protein_coding_genes_v19_+strand.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
parental_gencode = read.table("../data/gencode_parental_alleles_FULL.txt", header = TRUE, sep = "\t")
synonymous_parental_vars = subset(parental_gencode, vep_consequence == "synonymous_variant")

context = get_context(paste0("chr", synonymous_parental_vars$chr), synonymous_parental_vars$pos, gencode)




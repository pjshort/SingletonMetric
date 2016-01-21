# Null model for de novo regulatory mutations

# Probability of mutation in middle base of each trinucleotide is calculated from empirical data.
# The probabilities are assumed poisson, independent and summed across each sequence
library(GenomicRanges)

mut_rates <- read.table("../data/forSanger_1KG_mutation_rate_table.txt", header=TRUE)

### get the trinucleotide context

get_tri = function(CNE_idx, pos, CNEs) {
  start = pos - CNEs$start[CNE_idx] - 1
  end = pos - CNEs$start[CNE_idx] + 1
  chr = CNEs$chr[CNE_idx]
  seq = CNEs$seq[CNE_idx]

  if (start < 0 | end > nchar(seq)) {
    # retrieve using getSeq
    tri = getSeq(Hsapiens, chr, pos - 1, pos + 1)  # uses absolute start/stop
  } else {
    tri = substr(seq, start, end)  # use relative position
  }

}

get_context = function(chromosomes, positions, intervals = NULL) {
  # faster version using intervals with $seq column. must be sure that all chr, pos are contained in this set of intervals

  p = GRanges(seqnames=Rle(chromosomes), ranges = IRanges(start = positions, end = positions + 1))

  if (!is.null(intervals)) {  # the fast way
    i = GRanges(seqnames=Rle(intervals$chr), ranges = IRanges(start = intervals$start, end = intervals$stop))
  } else {
    # the SLOW way using hg19 sequence retrieval
  }

  # find overlap between denovos and annotated CNEs
  interval_hits_idx = findOverlaps(p, i, select = "first")

  context = mapply(get_tri, interval_hits_idx, positions, MoreArgs = list("CNEs" = intervals))
  context = sapply(context, as.character)

  return(context)

}



### SNP MODEL - based on sequence context ####

# the indel null model is inferred directly from global snp mutation rate and data - see rupit_core.R for generate_indel_null function

p_base <- function(from){

  # probability of mutation from base in position 2 to any of three other possible bases

  p = mut_rates$mu_snp[c(mut_rates$from == from)]
  return(sum(p))
}

p_position <- function(sequence, normalize = FALSE){

  # return vector length nchar(sequence) - 2 with probability of mutation at each base

  sequence = as.character(sequence)
  p = sapply(seq(1:nchar(sequence)-2), function(z) p_base(substr(sequence, z, z+2)))
  if (normalize == FALSE){
    return(p)
  } else {
    return(p/sum(p))
  }
}

p_sequence <- function(sequence){

  # sum p_all across each trinucleotide sliding along sequence

  p = p_position(sequence)
  return(sum(as.numeric(p)))

}

sample_alt <- function(ref_tri) {

  # given a ref sequence and mutated position, sample a new alt based on null model (middle base changes)

  m = mut_rates[which(mut_rates$from == ref_tri), c("to", "mu_snp")]
  tri_alt = sample(m$to, size=1, prob=m$mu_snp)
  return(as.character(substr(tri_alt,2,2)))
}


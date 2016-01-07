# retrieve CADD scores for a set of alleles using the pysam library
# the entire CADD genome is stored in tabix file (on the farm)

import pysam
import argparse


def get_options():
    """ get the command line options
    """

    parser = argparse.ArgumentParser(description="Get CADD scores for set of alleles.")
    parser.add_argument("--tabix", default=sys.stdout,
        help="location of genome-wide CADD score tabix file")
    parser.add_argument("--variants", default=sys.stdout,
        help="location of variants in chr\tpos\tref\talt (minimal vcf) format")
    parser.add_argument("--variants_out", default=sys.stdout,
        help="path to send the list of variants with CADD scores to.")
    parser.add_argument("--score", default=sys.stdout,
        help="Which score to compute. Accepts 'CADD' or 'Genomiser'.")
    args = parser.parse_args()

    return args

def get_regions(regions_path):
    """ input: args.regions - create a dictionary from chr to tuple of positions for regions
    """
    regions = {}

    with open(regions_path, "r") as myregions:
        myregions.readline()  # skip the header
        for line in myregions:
            line = line.split("\t")
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            if chr in regions:
                regions[chr].append((start, end))
            else:
                regions[chr] = [(start, end)]

    return regions


def get_variants(variants_path):
    chr = []
    pos = []
    ref = []
    alt = []
    with open(variants_path, "r") as myvariants:
        myvariants.readline()  # skip the header
        for line in myvariants:
            line = line.split("\t")
            chr.append(line[0])
            pos.append(int(line[1]))
            ref.append(line[2])
            alt.append(line[3].rstrip())

    return chr, pos, ref, alt

def parse_cadd_tabix(args):
  tabixfile = pysam.Tabixfile(args.tabix)
  #regions = get_regions(args$regions)
  chr, pos, ref, alt = get_variants(args.tabix)

  myfile = open(args$variants_out)
  myfile.write("\t".join(["chr", "pos", "ref", "alt", "unscaled_CADD", "scaled_CADD\n"]))
  for c, p, r, a in zip(chr, pos, ref, alt):

    t = tabixfile.fetch(c, p-1, p)
    for line in t:
      line = line.split("\t") # chr, pos, ref, alt, unscaled CADD, scaled CADD
      if (r == line[2]) & (a == line[3]):  # matches ref and alt
        myfile.write("\t".join([c,str(p),r,a,str(line[4]), str(line[5]) + "\n"]))
        break
      else:
        pass
  myfile.close()

  return

def parse_genomiser_tabix(args):
  tabixfile = pysam.Tabixfile(args.tabix)
  #regions = get_regions(args$regions)
  chr, pos, ref, alt = get_variants(args.variants)

  myfile = open(args$variants_out)
  myfile.write("\t".join(["chr", "pos", "ref", "alt", "genomiser_ReMM"]))
for c, p, r, a in zip(chr, pos, ref, alt):

  t = tabixfile.fetch(c, p-1, p)
  for line in t: # only a single line, not allele specific
    line = line.split("\t") # chr, pos, ref, alt, unscaled CADD, scaled CADD
    myfile.write("\t".join([c,str(p),r,a,str(line[2]) + "\n"]))

myfile.close()

if __name__ == "main":
  args = get_options()

  if args.score == "CADD":
    parse_cadd_tabix(args)
  elif args.score == "Genomiser":
    parse_genomiser_tabix(args)

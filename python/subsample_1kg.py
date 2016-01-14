# subsample 1000 genomes regions

thousand_genomes_dir = "/lustre/scratch113/projects/ddd/resources/1000G_MAF_20141114"

# sub sample non-coding variants from 1000 genomes

import pysam
import argparse
import sys
import random

chromosome_sizes = {"1":	249250621, "2":	243199373, "3":	198022430, "4": 191154276, "5": 180915260,
"6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516,
"12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210,
"18": 78077248, "19": 59128983, "20": 63025520, "21": 48129895, "22": 51304566}

def get_options():
    """ get the command line options
    """

    parser = argparse.ArgumentParser(description="Get CADD scores for set of alleles.")
    parser.add_argument("--chr", default="ALL",
        help="Select a chromosome to sample from (default is ALL).")
    parser.add_argument("--variants_per_chr", default=10000,
        help="Number of variants to sample per chromosome (default is 10,000)")
    parser.add_argument("--variants_out", default=sys.stdout,
        help="Where to save the output variants.")
    args = parser.parse_args()

    return args


def max_af(var, allele_number):
  entries = var.split(";")
  af_list = []
  for e in entries:
    af = float(e.split("=")[1].split(",")[allele_number])
    af_list.append(af)

  return max(af_list)

def subsample_1kg(args):

  if args.chr == "ALL":
    chromosomes = range(1,23)
  else:
    chromosomes = [args.chr]

  myfile = open(args.variants_out, 'w')
  myfile.write("\t".join(["chr", "pos", "ref", "alt", "MAX_1kg_MAF\n"]))

  for chr in chromosomes:
    print 'Sampling 10,000 variants from 1000 Genomes (chromosome %i)' % chr
    tabixfile = pysam.Tabixfile("/lustre/scratch113/projects/ddd/resources/1000G_MAF_20141114/ALL.chr%i.vcf.gz" % chr)
    count = 0
    chr_length = chromosome_sizes[str(chr)]
    while count < args.variants_per_chr:
      pos = random.randint(1, chr_length + 1)
      t = tabixfile.fetch(chr, pos-1, pos)
      vars = list(t)
      if len(vars) == 0: # no 1kg var at this position
        continue
      else:
        v = vars[0].split("\t") # chr, pos, id, ref, alt, NA, NA, INFO

        if len(v[4].split(",")) > 1:  # more than one alt
          idx = random.randint(0, len(v[4].split(",")) - 1)
          v[4] = v[4].split(",")[idx]  # take nth alt allele
          max_AF = max_af(v[7], idx)
        else:
          max_AF = max_af(v[7], 0)  # only one alt allele at ref pos

        myfile.write("\t".join([str(chr), str(pos), v[3], v[4], str(max_AF) + "\n"]))
        count += 1

  myfile.close()

if __name__ == '__main__':

  thousand_genomes_dir = "/lustre/scratch113/projects/ddd/resources/1000G_MAF_20141114"
  chromosome_sizes = {"1":	249250621, "2":	243199373, "3":	198022430, "4": 191154276, "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516,"12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210,"18": 78077248, "19": 59128983, "20": 63025520, "21": 48129895, "22": 51304566}

  args = get_options()
  subsample_1kg(args)

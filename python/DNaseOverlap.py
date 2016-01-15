import pysam
import argparse
import sys
import os

def get_options():
    """ get the command line options
    """

    parser = argparse.ArgumentParser(description="Get CADD scores for set of alleles.")
    parser.add_argument("--variants", default="/path/to/input/variants",
        help="location of variants in chr\tpos\tref\talt (minimal vcf) format")
    parser.add_argument("--variants_out", default=sys.stdout,
        help="path to send the list of variants with DNASE overlap binary values.")
    parser.add_argument("--roadmap_epigenome_ids", default = ["E002", "E010", "E053", "E072", "E081", "E082", "E083"],
    help = "list of roadmap epigenome project ideas in form E###.")
    args = parser.parse_args()

    return args


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


def check_dnase_overlap(chr, pos, ref, alt, id, TABIX_DIR):

  tabixfile = pysam.Tabixfile(TABIX_DIR + "regions_enh_%s.bed.gz" % id) # e.g. E081 is male fetal brain

  overlap = []
  for c, p, r, a in zip(chr, pos, ref, alt):

    t = tabixfile.fetch("chr" + c, p-1, p)

    if len(list(t)) > 0:  # overlaps entry
      overlap.append(1)
    else:  # no overlap
      overlap.append(0)

  return overlap


if __name__ == "__main__":
  args = get_options()
  chr, pos, ref, alt = get_variants(args.variants)
  TABIX_DIR = "/lustre/scratch113/projects/ddd/users/ps14/REP/"

  try:
    id_list = [line.rstrip() for line in open(args.roadmap_epigenome_ids)]
  except NameError:
    pass  # take default (brain and brain developmental tissues)

  overlap_list = []
  for id in id_list:
    if !os.path.isfile("regions_enh_%s.bed.gz" % id):
      print "No DNase Peaks for %s. Skipping and moving to the next tissue." % id
      continue
    print "Intersecting parental alleles with %s DNase peaks." % id
    overlap = check_dnase_overlap(chr, pos, ref, alt, id, TABIX_DIR)
    overlap_list.append(overlap)

  myvariants = open(variants_path, "r")
  variant_header =  myvariants.readline().rstrip()
  variants = myvariants.readlines()

  #lines = [variants] + overlap_list
  myfile = open(args.variants_out, 'w')

  # write header
  header = variant_header + "\t" + "\t".join(id_list) + "\n"
  myfile.write(header)

  # write lines
  i = 0
  for overlaps in zip(*overlap_list): # the * unpacks the list of lists
    var = variants[i].rstrip()
    myfile.write(var + "\t".join(str(o) for o in overlaps) + "\n")
    i += 1

  print 'Finished!'




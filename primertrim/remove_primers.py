import sys
import argparse
import os

class FastqRead(object):
   def __init__(self, read):
      self.desc, self.seq, self.qual = read
   def __repr__(self):
      return self.desc + "\n" + self.seq + "\n+\n" + self.qual + "\n"

def _grouper(iterable, n):
   "Collect data into fixed-length chunks or blocks"
   # grouper('ABCDEFG', 3) --> ABC DEF
   args = [iter(iterable)] * n
   return zip(*args)

def parse_fastq(f):
   for desc, seq, _, qual in _grouper(f, 4):
      desc = desc.rstrip()[1:]
      seq = seq.rstrip()
      qual = qual.rstrip()
      yield desc, seq, qual

def main(argv=None):
   p = argparse.ArgumentParser()
   p.add_argument("-r", "--read", help="remove_fwd_primer_from_R2 OR remove_rev_primer_from_R1")
   p.add_argument("-i", "--in_fastq_FP", help="PATH to input fastq to be trimmed")
   p.add_argument("-o", "--out_fastq_FP", help="PATH to output fastq after trimming")
   args = p.parse_args(argv)

   which_primer = args.read
   in_fastq_FP = args.in_fastq_FP
   out_fastq_FP = args.out_fastq_FP

   if which_primer == "remove_fwd_primer_from_R2":
      ### if input or output fastq file name ends with R1.fastq, then something is wrong
      if args.in_fastq_FP.endswith("R1.fastq") or out_fastq_FP.endswith("R1.fastq"):
         print("Error: R2 sequences expected in both input and output")
         sys.exit(0)
      elif in_fastq_FP.endswith("R2.fastq") and out_fastq_FP.endswith("R2.fastq"):
         primer = "TTACTTCCTCTAAATGACCAAG" ### rev complement of ITS1F primer	
      else: 
         print("Error: invalid input/output fastq file paths")
         sys.exit(0)
   elif which_primer == "remove_rev_primer_from_R1":
	### if input or output fastq file name ends with R2.fastq, then something is wrong
      if in_fastq_FP.endswith("R2.fastq") or out_fastq_FP.endswith("R2.fastq"):
         print("Error: R1 sequences expected in both input and output")
         sys.exit(0)
      elif in_fastq_FP.endswith("R1.fastq") and out_fastq_FP.endswith("R1.fastq"):
         primer = "GCATCGATGAAGAACGCAGC" ### rev complement of ITS2R primer    
      else: 
         print("Error: invalid input/output fastq file paths")
         sys.exit(0)
   else: 
      print("Error: invalid argument.\nProvide either 'remove_fwd_primer_from_R2' or 'remove_rev_primer_from_R1'")

   log_file_base = os.path.splitext(out_fastq_FP)[0]
   log_file_FP = log_file_base + ".log" 
   in_fastq = open(in_fastq_FP)
   reads = (FastqRead(x) for x in parse_fastq(in_fastq))

   with open(out_fastq_FP, "w") as out, open(log_file_FP, "w") as log:
      for read in reads:
         if primer in read.seq:
            primer_loc = read.seq.index(primer) ### returns first position of primer sequences match
            trimmed_seq = read.seq[:primer_loc]
            trimmed_qual = read.qual[:primer_loc]
         else:
            primer_loc = "NA"
            trimmed_seq = read.seq
            trimmed_qual = read.qual
         out.write("@%s\n%s\n+\n%s\n" % (read.desc, trimmed_seq, trimmed_qual))
         log.write("%s\t%s\n" % (read.desc, primer_loc))
   in_fastq.close()


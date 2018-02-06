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
   p.add_argument(
      "-i", "--in_fastq", type=argparse.FileType('r'),
      help="PATH to input fastq to be trimmed")
   p.add_argument(
      "-o", "--out_fastq", type=argparse.FileType('w'),
      help="PATH to output fastq after trimming")
   p.add_argument(
      "--log", type=argparse.FileType('w'),
      help="Log file to record location of primers detected")
   args = p.parse_args(argv)

   which_primer = args.read
   in_fastq = args.in_fastq
   out_fastq = args.out_fastq
   log = args.log

   if which_primer == "remove_fwd_primer_from_R2":
         primer = "TTACTTCCTCTAAATGACCAAG" ### rev complement of ITS1F primer	
   elif which_primer == "remove_rev_primer_from_R1":
         primer = "GCATCGATGAAGAACGCAGC" ### rev complement of ITS2R primer    

   reads = (FastqRead(x) for x in parse_fastq(in_fastq))

   for read in reads:
       if primer in read.seq:
          primer_loc = read.seq.index(primer) ### returns first position of primer sequences match
          trimmed_seq = read.seq[:primer_loc]
          trimmed_qual = read.qual[:primer_loc]
       else:
          primer_loc = "NA"
          trimmed_seq = read.seq
          trimmed_qual = read.qual
       out_fastq.write("@%s\n%s\n+\n%s\n" % (read.desc, trimmed_seq, trimmed_qual))
       log.write("%s\t%s\n" % (read.desc, primer_loc))


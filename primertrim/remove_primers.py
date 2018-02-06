import sys
import argparse
import os

class FastqRead(object):
   def __init__(self, read):
      self.desc, self.seq, self.qual = read

   def trim(self, idx):
      if idx is None:
         return self
      else:
         return FastqRead((self.desc, self.seq[:idx], self.qual[idx]))

   def format_fastq(self):
      return "@{0}\n{1}\n+\n{2}\n".format(self.desc, self.seq, self.qual)

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

class ExactMatcher(object):
    def __init__(self, queryset):
        self.queryset = queryset

    def find_match(self, seq):
        return find_queryset(seq, self.queryset)

def find_queryset(s, queryset):
    idx = -1
    for query in queryset:
        idx = s.find(query)
        if idx != -1:
            break
    if idx == -1:
        idx = None
    return idx

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

   m = ExactMatcher([primer])
   reads = (FastqRead(x) for x in parse_fastq(in_fastq))

   for read in reads:
       idx = m.find_match(read.seq)
       trimmed_read = read.trim(idx)
       out_fastq.write(trimmed_read.format_fastq())

       if idx is None:
           primer_loc = "NA"
       else:
           primer_loc = str(idx)
       log.write("%s\t%s\n" % (trimmed_read.desc, primer_loc))

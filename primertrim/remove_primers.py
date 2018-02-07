import sys
import argparse
import os

class FastqRead(object):
   def __init__(self, desc, seq, qual):
      self.desc = desc
      self.seq = seq
      self.qual = qual

   def trim(self, idx):
      if idx is None:
         return self
      else:
         return self.__class__(self.desc, self.seq[:idx], self.qual[:idx])

   def format_fastq(self):
      return "@{0}\n{1}\n+\n{2}\n".format(self.desc, self.seq, self.qual)

   @classmethod
   def parse(cls, f):
      for desc, seq, _, qual in _grouper(f, 4):
         desc = desc.rstrip()[1:]
         seq = seq.rstrip()
         qual = qual.rstrip()
         yield cls(desc, seq, qual)

def _grouper(iterable, n):
   "Collect data into fixed-length chunks or blocks"
   # grouper('ABCDEFG', 3) --> ABC DEF
   args = [iter(iterable)] * n
   return zip(*args)

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
   p.add_argument(
      "primer",
      help="Primer sequence to be trimmed")
   p.add_argument(
      "-i", "--input_fastq", type=argparse.FileType('r'),
      help="Input FASTQ file to be trimmed (default: standard input)")
   p.add_argument(
      "-o", "--output_fastq", type=argparse.FileType('w'),
      help="Output fastq after trimming (default: standard output)")
   p.add_argument(
      "--log", type=argparse.FileType('w'),
      help="Log file to record location of primers detected (default: not written)")
   args = p.parse_args(argv)

   if args.input_fastq is None:
      args.input_fastq = sys.stdin

   if args.output_fastq is None:
      args.output_fastq = sys.stdout

   m = ExactMatcher([args.primer])
   reads = FastqRead.parse(args.input_fastq)

   for read in reads:
       idx = m.find_match(read.seq)
       trimmed_read = read.trim(idx)
       args.output_fastq.write(trimmed_read.format_fastq())

       if args.log:
          if idx is None:
             primer_loc = "NA"
          else:
             primer_loc = str(idx)
          args.log.write("{0}\t{1}\n".format(trimmed_read.desc, primer_loc))

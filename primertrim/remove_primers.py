import itertools
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

class Matcher(object):
    def __init__(self, queryset):
        self.queryset = queryset.copy()

    def find_match(self, seq):
        idx = -1
        for query in self.queryset:
            idx = seq.find(query)
            if idx != -1:
                break
        if idx == -1:
            idx = None
        return idx

class CompleteMatcher(Matcher):
    def __init__(self, queryset, max_mismatch):
        super().__init__(queryset)

        if max_mismatch > 0:
           n = 1
           while n <= max_mismatch:
              self.queryset.extend(self._mismatched_queries(queryset, n))
              n += 1

    def _mismatched_queries(self, queryset, n_mismatch):
        # This algorithm is terrible unless the number of mismatches is very small
        assert(n_mismatch in [1, 2, 3])
        for query in queryset:
            idx_sets = itertools.combinations(range(len(query)), n_mismatch)
            for idx_set in idx_sets:
                # Replace base at each position with a literal "N", to match
                # ambiguous bases in the reference
                yield replace_with_n(query, idx_set)
                # Change to list because strings are immutable
                qchars = list(query)
                # Replace the base at each mismatch position with an
                # ambiguous base specifying all possibilities BUT the one
                # we see.
                for idx in idx_set:
                    qchars[idx] = AMBIGUOUS_BASES_COMPLEMENT[qchars[idx]]
                    # Expand to all possibilities for mismatching at this
                    # particular set of positions
                for query_with_mismatches in deambiguate(qchars):
                    yield query_with_mismatches

def replace_with_n(seq, idxs):
    chars = list(seq)
    for idx in idxs:
        chars[idx] = "N"
    return "".join(chars)

AMBIGUOUS_BASES = {
    "T": "T",
    "C": "C",
    "A": "A",
    "G": "G",
    "R": "AG",
    "Y": "TC",
    "M": "CA",
    "K": "TG",
    "S": "CG",
    "W": "TA",
    "H": "TCA",
    "B": "TCG",
    "V": "CAG",
    "D": "TAG",
    "N": "TCAG",
    }


# Ambiguous base codes for all bases EXCEPT the key
AMBIGUOUS_BASES_COMPLEMENT = {
    "T": "V",
    "C": "D",
    "A": "B",
    "G": "H",
    }


def deambiguate(seq):
    nt_choices = [AMBIGUOUS_BASES[x] for x in seq]
    return ["".join(c) for c in itertools.product(*nt_choices)]


COMPLEMENT_BASES = {
    "T": "A",
    "C": "G",
    "A": "T",
    "G": "C",
    }


def reverse_complement(seq):
    rc = [COMPLEMENT_BASES[x] for x in seq]
    rc.reverse()
    return ''.join(rc)

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
   p.add_argument(
      "-n", "--num_mismatches", type=int, default=0,
      help="Number of mismatches to primer allowed (default: %(default)s)")
   args = p.parse_args(argv)

   if args.input_fastq is None:
      args.input_fastq = sys.stdin

   if args.output_fastq is None:
      args.output_fastq = sys.stdout

   queryset = deambiguate(args.primer)
   m = CompleteMatcher(queryset, args.num_mismatches)

   for read in FastqRead.parse(args.input_fastq):
       idx = m.find_match(read.seq)
       trimmed_read = read.trim(idx)
       args.output_fastq.write(trimmed_read.format_fastq())

       if args.log:
          if idx is None:
             primer_loc = "NA"
          else:
             primer_loc = str(idx)
          args.log.write("{0}\t{1}\n".format(trimmed_read.desc, primer_loc))

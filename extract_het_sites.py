import sys
from vcf import Reader
import argparse

parser = argparse.ArgumentParser(description='''
Extracts all het-sites for bi-allelic markers for one sample.
''')

parser.add_argument('vcf', nargs=1,
                    help='VCF-file of called variants.')
parser.add_argument('sample', nargs=1,
                    help='Sample to extract het-sites for.')

parser.add_argument('--min-read-count', 
                    default=10, action='store', type=int,
                    help='Minimum number of reads covering each allele.')

args = parser.parse_args()

# FIXME: handle missing files more gracefully
if args.vcf[0] == '-':
    vcffile = sys.stdin
else:
    vcffile = open(args.vcf[0])

sample = args.sample[0]

for record in Reader(vcffile):
    if record.FILTER:
        continue
    if len(record.ALT) != 1:
        continue
    
    call = record.genotype(sample)
    if call.is_het and len(call.data.AD) == 2:
        count1, count2 = call.data.AD
        if min(count1,count2) > args.min_read_count:
            print "%s\t%d\t%s\t%d\t%s\t%d" % (
                record.CHROM, record.POS, 
                record.REF, count1, 
                record.ALT[0], count2)
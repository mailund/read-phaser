import sys
from vcf import Reader

if len(sys.argv) != 3:
    print 'Usage:', sys.argv[0], 'vcf-file sample'
    sys.exit(1)
    
vcffile = sys.argv[1]
sample = sys.argv[2]

for record in Reader(open(vcffile)):
    if record.FILTER:
        continue
    if len(record.ALT) != 1:
        continue # FIXME: this isn't necessary but makes the output 
                 # slightly easier to work with
    
    call = record.genotype(sample)
    if call.is_het and len(call.data.AD) == 2:
        count1, count2 = call.data.AD
        if min(count1,count2) > 5:
            print record.CHROM, record.POS, 
            print record.REF, count1, 
            print record.ALT[0], count2
import sys
import pysam
import argparse

parser = argparse.ArgumentParser(description='Phase pairs of het-sites.')

parser.add_argument('hets', nargs=1,
                    help='List of het-sites (use - for stdin).')
parser.add_argument('bam', nargs=1,
                    help='Bam-file of reads (assumed covering a single '
                         'individual).')

parser.add_argument('--min-read-count', 
                    default=10, action='store', type=int,
                    help='Minimum number of reads covering each phase of a '
                    'pair of het sites.')
parser.add_argument('--max-marker-distance', 
                    default=1000, action='store', type=int,
                    help='Maximal distance between two het-markers that '
                    'will be considered.')

args = parser.parse_args()

# Threshold constants
MIN_COUNT_THRESHOLD = args.min_read_count
MAX_MARKER_DISTANCE = args.max_marker_distance

# FIXME: handle missing files more gracefully
bamfile = pysam.Samfile(args.bam[0], 'rb')
if args.hets[0] == '-':
    hets = sys.stdin
else:
    hets = open(args.hets[0])

# FIXME: THIS SHOULD PROBABLY BE STREAMED SO WE DON'T HAVE TO READ IN ALL
# THE HET SITES BEFORE WE START OUTPUTTING PHASED PAIRS
var_sites = []
for line in hets:
    x = line.split()
    chrom,pos = x[0],int(x[1])-1
    var_sites.append( (chrom,pos) )
hets.close()

def get_pileup_column(site):
    '''Extract a specific column of the alignment and get 
    all pileup reads there. Extract the query name and the called
    nucleotide'''
    pileup = bamfile.pileup(site[0], site[1], site[1]+1)
    for col in pileup:
        if col.pos == site[1]:
            return [(read.alignment.qname, read.alignment.seq[read.qpos])
                    for read in col.pileups]
    return []


for i in xrange(len(var_sites)):

    i_reads = dict( (name,allele) 
                    for name,allele in get_pileup_column(var_sites[i]))

    # Collect the j indices relevant for this i index.
    valid_j_indices = []
    for j in xrange(i+1,len(var_sites)):

        if var_sites[j][0] != var_sites[i][0]:
            break # different chromosome, no need to continue here
    
        if var_sites[j][1] - var_sites[i][1] > MAX_MARKER_DISTANCE:
            break
            
        valid_j_indices.append(j)
    
    # now phase the i/j sites.
    for j in valid_j_indices:
        j_reads = dict( (name,allele) 
                        for name,allele in get_pileup_column(var_sites[j]))

        names_overlap = set(i_reads.keys()).intersection(j_reads.keys())
        if len(names_overlap) == 0:
            continue

        hap_type_count = {}
        for qname in names_overlap:
            i_allele = i_reads[qname]
            j_allele = j_reads[qname]
            key = (i_allele,j_allele)
            try:
                hap_type_count[key] += 1
            except:
                hap_type_count[key] = 1
      
        likely_calls = dict( (k,v) for k,v in hap_type_count.items()
                                   if v >= MIN_COUNT_THRESHOLD )
        if len(likely_calls) != 2:
            continue

        (gtype1,count1), (gtype2,count2) = likely_calls.items()

        # check that we actually have two alleles at both sites
        alleles1 = set([gtype1[0],gtype2[0]])
        alleles2 = set([gtype1[1],gtype2[1]])
        if len(alleles1) != 2 or len(alleles2) != 2:
            continue
        
        print "%s\t%d\t%d\t%s%s\t%d\t%s%s\t%d" % (
            var_sites[i][0], var_sites[i][1], var_sites[j][1],
            gtype1[0],gtype1[1], count1,
            gtype2[0],gtype2[1], count2)
